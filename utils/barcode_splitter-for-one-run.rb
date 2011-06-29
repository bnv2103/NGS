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

require 'zlib'

def main
  barcode = ARGV[0] # comma-delimited csv file
  ## targetfastq = ARGV[1]  # barcode is the first 6bp of the reads
  outprefix = ARGV[1]
  inputdir = ARGV[2]

  barcodesize = 6  # barcode size

  coding = readBar(barcode) # return a hash {lane => { code => sample  }}
  

### assume the original fastq files are: s_1_1.fastq for lane 1

  multiplex = {}

  coding.each do |lane, ch|
    multiplex[lane] = {}
    ch.each do |str,sampleID|
      mutate1(str).each do |strmut|   # only allow hamming distance of 1
        multiplex[lane][strmut] = sampleID
      end
    end
  end
##  $stderr.puts "allowed multiplexing codes: #{multiplex.keys.size}"
  $stderr.puts "multiplex mapping: #{multiplex}"

##  outputio["discarded"] = File.new(targetfastq + "_discarded.fastq", "w" )
  assignment = decode(inputdir, multiplex, outprefix, barcodesize)
  
end

def decode(inputdir, multiplex, outprefix, barcodesize)
  
  # get the list of files in the dir
  #   puts inputdir
  barcodefq = Dir.new(inputdir).select {|a| a.match(/s\_\d+\_2\.\S+/) }
  
  $stderr.puts "barcode files: \n#{barcodefq.join("\n")}"

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
  multiplex[lane].values.sort.each do |sampleID|
    outName = outprefix + "_lane_#{lane}_#{sampleID}_#{ends}.fastq"
    outio[sampleID] = File.new(outName,'w')
  end
  
  outio[:discard] = File.new(outprefix + "_lane_#{lane}_#{ends}.unknown.fastq", 'w')
  
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
      ndecode += 1
    else # discarded
      outio[:discard].puts tunit
      ndiscard += 1
    end
  end

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


  coding = {}
  File.new(b, 'r').each do |line|
    next if line.match(/^#/) # header line

    cols = line.chomp.split(',')
    run, lane, sampleID, code = cols[0].strip, cols[1].strip, cols[2].strip, cols[4].strip

    if !coding.key?(lane)
      coding[lane] = {}
    end
    coding[lane][code] = sampleID
  end
  return coding
end

main()


## NOte: optimal file reading:
# File.open(ARGV.first, 'r') do |fi|
#  fsize = fi.size
#  fi.read(fsize).lines { |l| 
#  }
# end
