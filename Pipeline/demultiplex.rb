## split fastq files based on barcode,
## barcode matching using hamming distance ( <=2 to real and >2 to the other)
##  alternatively, we can generate a set of possible codes with hd <=2 and then use hash

#require "narray"

#class String   ## hamming distance and xor function from Kirk Haines
#  def xor(other)
#    if other.empty?
#      self
#    else
    #  left = self
#      right = other

#      if left.length < right.length
#        n,r = right.length.divmod(left.length)
#        left = left * n + left[0, r]
#      elsif right.length < left.length
#        n,r = left.length.divmod(right.length)
#        right = right * n + right[0, r]
#      end

#    (NArray.to_na(self, "byte") ^ NArray.to_na(other, "byte")).to_s
#    end
#  end
  
#  def hamming_distance(other)
#    self.xor(other).tr("\x00",'').length
#  end
# end

# Output fastq file naming convention RUNNAME_laneX_SAMPLEID_ENDS.fastq  (where RUNNAME,laneX,SAMPLEID have no _ and ENDS = 1 (SE) OR 3(PE)
require 'zlib'

def main
  inputdir = ARGV[0]
  outdir = ARGV[1]#.tr("_", "-")         #replace all _ with - in RUNNAME
  prefix = ARGV[2]
  barcode = ARGV[3] # comma-delimited csv file
  
  nt = ARGV[4]  # max threads

  ## targetfastq = ARGV[1]  # barcode is the first 6bp of the reads
  
  if nt == nil
    nt = 4
  else
    nt = nt.to_i
  end



  barcodesize = 6  # barcode size

  coding = readBar(barcode) # return a hash {lane => { code => sample  }}
 
### assume the original fastq files are: s_1_1.fastq for lane 1

  multiplex = {}
  $reads_lane ={}

  coding.each do |lane, ch|
    multiplex[lane] = {}
    ch.each do |str,sampleID|
      mutate1(str).each do |strmut|   # only allow hamming distance of 1
        multiplex[lane][strmut] = sampleID.tr("_", "-") + "_" + str		#replace all _ with - in sampleID
      end
    end
  end
  ##  $stderr.puts "allowed multiplexing codes: #{multiplex.keys.size}"
  $stderr.puts "multiplex mapping: #{multiplex}"

  outprefix = outdir + "/" + prefix 
  stat = sanity_check(inputdir, outprefix, coding) 
  if stat == 0
	$stderr.puts "ERROR: Run Failed Sanity Chcheck "
  else
	$stderr.puts "Sanity Check Sucessful "
  end

 ##  outputio["discarded"] = File.new(targetfastq + "_discarded.fastq", "w" )
  assignment = decode(inputdir, multiplex, outprefix, barcodesize, nt)
  
  ## print summary on reads - mean and stdev per lane
  $reads_lane.each do |lane,arr_samples| 
         $stderr.puts "Lane #{lane} \t" + arr_samples.join("\t") 
  end
  $stderr.puts "\n"
  $reads_lane.each do |lane,arr_samples|
         mean = arr_samples.inject{ |sum, el| sum + el }.to_f / arr_samples.size
         variance = arr_samples.inject(0) { |variance, x| variance += (x - mean) ** 2 }
         stddev = Math.sqrt(variance/(arr_samples.size-1))
         $stderr.puts "Lane #{lane}\tMean: #{mean.to_i}\tStdev: #{stddev.to_i}"
  end


end

def decode(inputdir, multiplex, outprefix, barcodesize, nt )
  
  # get the list of files in the dir
  #   puts inputdir
  barcodefq = Dir.new(inputdir).select {|a| a.match(/s\_\d+\_2\.fastq$/) }
  
  $stderr.puts "barcode files: \n#{barcodefq.join("\n")}"

  nprocess = 0

  barcodefq.sort.each do |bfq|
    if bfq.match(/s\_(\d+)\_2\.(\S+)/)
      lane = "#{$1}"
      targetfq1 = "#{inputdir}/s_#{lane}_1.#{$2}"
      targetfq3 = "#{inputdir}/s_#{lane}_3.#{$2}"
      
      
      
      bfql = "#{inputdir}/#{bfq}"
      next unless multiplex.key?(lane)
      $stderr.puts "working on lane #{lane}"
      Process.fork do
        doSplit(bfql, targetfq1, lane, multiplex, outprefix, barcodesize) if File.exist?(targetfq1)
        doSplit(bfql, targetfq3, lane, multiplex, outprefix, barcodesize) if File.exist?(targetfq3)
      end

      nprocess += 1

      if nprocess >= nt # need wait
        Process.waitall
        nprocess = 0 
      end
    end
  end
  Process.waitall
end


def doSplit(bfq, targetfq, lane, multiplex, outprefix, barcodesize)
  
  ndecode = 0
  ndiscard = 0
  
  ends = "1"
  if targetfq =~ /s\_\d+\_(\d+)\.f(\S+)$/
    ends = $1
  end


  outio = {}
  reads_sample ={}
  multiplex[lane].values.sort.each do |sampleID|
    outName = outprefix + "_lane#{lane}_#{sampleID}_#{ends}.fastq"
    outio[sampleID] = File.new(outName,'w')
    reads_sample[sampleID] =0;
 end
  
  outio[:discard] = File.new(outprefix + "_lane#{lane}_#{ends}.unknown.fastq", 'w')
  
  if bfq.match(/.gz$/)
    bio = Zlib::GzipReader.new(File.open(bfq))
  else
    bio = File.new(bfq, 'r')
  end

  if targetfq.match(/.gz$/)
    tio = Zlib::GzipReader.new(File.open(targetfq))
  else
    tio = File.new(targetfq, 'r')
  end
  
#   Zlib::GzipReader.new
 
     # chunk = 1000000  # 1 million lines each time
  while !bio.eof? # not end
    bunit = bio.take(4)
    tunit = tio.take(4)
    bc = bunit[1][0..barcodesize-1]
    
    if multiplex[lane].key?(bc) ## 
      
      sampleID = multiplex[lane][bc]
      
      outio[sampleID].puts tunit
      reads_sample[sampleID] += 1
      ndecode += 1
    else # discarded
      outio[:discard].puts tunit
      ndiscard += 1
    end
  end

 ### Output summary of reads per sample and reads per lane etc.

  $stderr.puts "#{targetfq}:\t#decoded=#{ndecode}\t#unknown=#{ndiscard}"
  read_summary_file =  File.open(outprefix + "_summary_lane#{lane}.stats",'a+')
  lane_out="lane"
  temp_out="#{lane}"
  $reads_lane[lane] = []

  reads_sample.each do |sample,reads|
	lane_out << "\t#{sample}"
	temp_out << "\t#{reads}"
	$reads_lane[lane].push(reads)	
  end
  lane_out << "\tDecoded"
  temp_out << "\t#{ndecode}"
  lane_out << "\tUnknown"
  temp_out << "\t#{ndiscard}"
  read_summary_file.puts lane_out
  read_summary_file.puts temp_out
  read_summary_file.close

  outio.values.each {|oio| oio.close}
  bio.close
  tio.close

end

def mutate2(string)
  ## mutation the string within hammingDist of 2
  hd1 = mutate1(string)
  hd2=[]
  hd1.each do |s|
    hd2.concat(mutate1(s))
  end
  return hd1.concat(hd2)
end

def mutate1(string)
  l = string.length
  hd1 = []
  0.upto(l-1) do |i|
    cp = String.new(string)
    ["A", "T", "G", "C", "."].each do |nt|
      cp[i] = nt
      hd1 << String.new(cp)
    end
  end
  return hd1
end


def readBar(b)
  # example:
  # Ac039eabxx, 1 ,15,hg18,TGACCA,Target,N,SE100,Erin Bush,Richard Mayeux
  # Ac039eabxx, 1 ,16,hg18,ACAGTG,Target,N,SE100,Erin Bush,Richard Mayeux
  # Ac039eabxx, 2 ,R1d492,bacteria,CGATGT,DNA Single End 50bp multiplexing=6,N,SE100,Erin Bush,Largus Angenent
  # Ac039eabxx, 2 ,R1d520,bacteria,TGACCA,DNA Single End 50bp multiplexing=6,N,SE100,Erin Bush,Largus Angenent

  #111012_SN828_0089_BD0CWFACXX    1       YELLOW  YELLOW  1       Bacteria        ATCACG  DNA     --      SE100   Erin Bush       PI      110926_DEREK_ASHLEY_24_BACTERIA_DNA_SE100_HISEQ aefranks@microbio.umass.edu, jward@microbio.umass.edu, la249@Cornell.edu

  coding = {}
  File.new(b, 'r').each do |line|
    next if line.match(/^#/) # header line
    

    cols = line.chomp.split(/\t/)

    if cols.size < 5 
      cols  = line.chomp.split(/\t/)
    end

    run, lane, sampleID, code = cols[0].strip, cols[1].strip, cols[4].strip, cols[6].strip
    
    next if code == "--"
    next if code == "-"
    
    if !coding.key?(lane)
      coding[lane] = {}
    end

	if coding[lane].key?(code)
                $stderr.puts "Duplicate barcodes in lane : #{lane}  for #{code} from sample #{coding[lane][code]}"
		exit
	end

## replace / with _, and space with _ in sampleID	
    coding[lane][code] = sampleID.tr("/", "_").tr(" ","_")
  end
  return coding
end

def sanity_check(inputdir, outprefix, coding)
 barcodefq = Dir.new(inputdir).select {|a| a.match(/s\_\d+\_2\.fastq.barcode-stats$/) }
 
 flag=1
 barcodefq.sort.each do |bfq|
	if bfq.match(/s\_(\d+)\_2\.(\S+)/)
		lane = "#{$1}"
		num_samples_this_lane = coding[lane].length
		count = 0
		all_codes=""
		f = File.new("#{inputdir}/#{bfq}","r")
		begin
		    while (line = f.readline)
		        line.chomp
			break if line =~ /\.\.\./
			cols = line.chomp.split(/\t/)
			all_codes <<  cols[0].strip 
			all_codes << "\t"
				count += 1	
		    end
		rescue EOFError
		    f.close
		end
#		print "\nlane = #{lane}\t all_codes = #{all_codes}"
		if count == num_samples_this_lane
			sanity_flag=1
			coding[lane].each_key do |code|
#	                        print  "\n\t#{code}"
				if all_codes.scan( "#{code}A" ).empty?		##Assuming one polyA base for HiSeq runs
#					print "\tfailed\t"
					sanity_flag=0
#                                       print  all_codes.scan("#{code}" )
					break
				else
#					print "\tpassed\t"
#					print  all_codes.scan("#{code}" )
				end
			end
			if sanity_flag == 1
				$stderr.puts "Sanity Check - Lane #{lane} : Successful"
			else
                        	$stderr.puts "Sanity Check - Lane #{lane} : Failed"
	                        flag=0
                	end
		else
			$stderr.puts "Sanity Check - Lane #{lane} : Failed"
			flag=0
		end
	end
 end
 return flag
end

main()


## NOte: optimal file reading:
# File.open(ARGV.first, 'r') do |fi|
#  fsize = fi.size
#  fi.read(fsize).lines { |l| 
#  }
# end
