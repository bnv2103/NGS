require 'csv'
$VERBOSE = nil

def main
  fq = ARGV[0] # sampl.fastq
  output = ARGV[1] #pdf file

  # `samtools view #{fq} > accepted_hits.sam`
  infile = File.new(fq, "r")
  outfile = File.new("ACTG.txt", "w")
  outfile1 = File.new("Qual.txt", "w")
  numReads = 0
  

  qual = Array.new(101){ |i| 0*i }
  valueA = Array.new(101) { |i| 0*i }
  valueT = Array.new(101) { |i| 0*i }
  valueC = Array.new(101) { |i| 0*i }
  valueG = Array.new(101) { |i| 0*i }

  while !infile.eof?
    
    # puts line
    # numReads = numReads + 1
    s = infile.take(1)
    # cols = line.chomp.split(' ')
    
    # readName = cols[0]
    # readSeq = cols[9]
    # readQual = cols[10]
    
    readName = s
    s = infile.take(1)
    s = s.to_s
    s = s[2..102]
    # puts s
    readSeq = s
    len = readSeq.length()
    # puts len
    s = infile.take(1)
    s = infile.take(1)
    s = s.to_s
    s = s[2..102]
    readQual = s
    numReads = numReads + 1
   # puts readQual
    for i in 0..readSeq.length()
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

    
   
    for i in 0..readQual.length()
      bp = readQual[i]
      if bp != nil
        q = readQual[i].ord - 33
        qual[i] = qual[i] + q
        # puts "Value of local variable is #{readQual[i]} #{q}"
      end
    end
   
   # puts qual[1]
   # outfile.puts "#{readName} \t\t #{qsum}"
  end

  numReads = numReads -1
  
  puts numReads 
  
  outfile.print "A \t"
  for i in 0..101
    outfile.print "#{valueA[i]} \t"
  end
  outfile.print "\n T \t"
  for i in 0..101
    outfile.print "#{valueT[i]} \t"
  end
  outfile.print "\n C \t"
  for i in 0..101
    outfile.print "#{valueC[i]} \t"
  end
  outfile.print "\n G \t"
  for i in 0..101
    outfile.print "#{valueG[i]} \t"
  end
  outfile.print "\n"
  # outfile.print "\n qual \t"
  for i in 0..101
     outfile1.print "#{qual[i]} \t"
   end
  outfile1.print "\n"
  outfile1.close





  infile.close

  outfile.close

  `Rscript /ifs/scratch/c2b2/ngs_lab/xs2182/code/QCplot.R #{numReads} #{output}`

end

main()
