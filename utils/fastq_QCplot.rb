require 'csv'

$VERBOSE = nil

## Plots the Quality Distribution, ACTG per cycle from raw fastq reads. Can take paired end file as 3rd argument. 
def main
  output = ARGV[0] #pdf file
  fq = ARGV[1] # sampl.fastq
  fq3 = ARGV[2] # sampl.fastq3 if there is 

  qground = 64 ## For illumina 1.3+, 1.5+ use -64, for illumina 1.8+ use -33 

  outfile = File.new("#{output}_ACTG.txt", "w")
  outfile1 = File.new("#{output}_Qual.txt", "w")
#  qsc_f = File.new("#{output}_QScore_f.txt", "w")
#  qsc_r = File.new("#{output}_QScore_r.txt", "w")
  qsc_f_hist = File.new("#{output}_QScore_f_hist.txt", "w")
  qsc_r_hist = File.new("#{output}_QScore_r_hist.txt", "w")
  utils = "/ifs/data/c2b2/ngs_lab/ngs/code/NGS/utils/"

  numReads = 0
  slen=0
  if fq3.nil? 
	pair = 1
  else
	if File.exists?(fq3)
		pair =2
	else
		pair =1
	end
  end
 
  qscore_bucket_fwd = Array.new(45){ |i| 0*i }
  qscore_bucket_rvs = Array.new(45){ |i| 0*i }
  infile = File.new(fq, "r")

while !infile.eof?
    s = infile.take(1)
    readName = s
    s = infile.take(1)
    s = s.to_s
    s = s[2..-5]
    readSeq = s

   if numReads ==0
     slen = s.size
     qual = Array.new(slen*pair){ |i| 0*i }
     valueA = Array.new(slen*pair) { |i| 0*i }
     valueT = Array.new(slen*pair) { |i| 0*i }
     valueC = Array.new(slen*pair) { |i| 0*i }
     valueG = Array.new(slen*pair) { |i| 0*i }
   end

    s = infile.take(1)
    s = infile.take(1)
    s = s.to_s
    s = s[2..-5]
    readQual = s
    numReads = numReads + 1
    for i in 0..(slen-1)		##Extract ACTG info
      bp = readSeq[i]
      if bp != nil
        if bp == "A"
          valueA[i] = valueA[i] + 1
        elsif bp =="T"
          valueT[i] = valueT[i] + 1
        elsif bp =="C"
          valueC[i] = valueC[i] + 1
        elsif bp =="G"
          valueG[i] = valueG[i] + 1
        end
      end
    end 
    for i in 0..(slen-1)		##Extract Quality distribution info
      bp = readQual[i]
      if bp != nil
        q = readQual[i].ord - qground
        qual[i] = qual[i] + q        
      end
    end 
    aqscore=0   		#Extract Qscore histogram info
    readQual.each_char do | item |
	aqscore += (item.ord - qground)
    end
        avgqscore = aqscore/slen
#        qsc_f.puts avgqscore    #Extract Qscore histogram info
        qscore_bucket_fwd[avgqscore]+=1
  end

if pair == 2
infile.close
infile = File.new(fq3, "r")
  while !infile.eof?
    s = infile.take(1)
    readName = s
    s = infile.take(1)
    s = s.to_s
    s = s[2..-5]
    readSeq = s
    s = infile.take(1)
    s = infile.take(1)
    s = s.to_s
    s = s[2..-5]
    readQual = s
    numReads = numReads + 1
    for j in 0..(slen -1)	##Extract ACTG info
      bp = readSeq[j]
	i = (j + slen)
      if bp != nil
        if bp == "A"
          valueA[i] = valueA[i] + 1
        elsif bp =="T"
          valueT[i] = valueT[i] + 1
        elsif bp =="C"
          valueC[i] = valueC[i] + 1
        elsif bp =="G"
          valueG[i] = valueG[i] + 1
        end
      end
    end 
    for j in 0..(slen-1)		##Extract Quality distribution info
      bp = readQual[j]
        i = (j + slen)
      if bp != nil
        q = readQual[j].ord - qground
        qual[i] = qual[i] + q         
      end
    end    
    aqscore=0
    readQual.each_char do | item |
        aqscore += (item.ord - qground)
    end
        avgqscore = aqscore/slen
#	qsc_r.puts avgqscore	#Extract Qscore histogram info
	qscore_bucket_rvs[avgqscore]+=1
  end
end

  numReads = numReads -1
  outfile.print "A \t"
  for i in 0..(slen*pair)
    outfile.print "#{valueA[i]} \t"
  end
  outfile.print "\n T \t"
  for i in 0..(slen*pair)
    outfile.print "#{valueT[i]} \t"
  end
  outfile.print "\n C \t"
  for i in 0..(slen*pair)
    outfile.print "#{valueC[i]} \t"
  end
  outfile.print "\n G \t"
  for i in 0..(slen*pair)
    outfile.print "#{valueG[i]} \t"
  end
  outfile.print "\n"
  # outfile.print "\n qual \t"
  for i in 0..(slen*pair)
     outfile1.print "#{qual[i]} \t"
  end

  qscore_bucket_fwd.each do | item |
	qsc_f_hist.puts item
  end
  qscore_bucket_rvs.each do | item |
        qsc_r_hist.puts item
  end

  outfile1.print "\n"
  outfile1.close
  infile.close
  outfile.close
  qsc_r_hist.close
  qsc_f_hist.close
#  qsc_r.close
#  qsc_f.close


puts "Rscript #{utils}/QCplot_qscore.R  #{numReads} #{output}.pdf #{output}_ACTG.txt #{output}_Qual.txt #{output}_QScore_f_hist.txt #{output}_QScore_r_hist.txt "

`Rscript #{utils}/QCplot_qscore.R #{numReads} #{output}.pdf   #{output}_ACTG.txt #{output}_Qual.txt #{output}_QScore_f_hist.txt #{output}_QScore_r_hist.txt `

end

main()
