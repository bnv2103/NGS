require 'csv'
$VERBOSE = nil

def main
  fq = ARGV[0] # bam file dir
  # output = ARGV[1] # QC.pdf
  puts fq
  `samtools view #{fq}/accepted_hits.bam > #{fq}/accepted_hits_QC.sam`
  infile = File.new("#{fq}/accepted_hits_QC.sam", "r")
  outfile = File.new("#{fq}/ACTG.txt", "w")
  outfile1 = File.new("#{fq}/Qual.txt", "w")
  
  numReads = 0
  arrLen = 101
  qual = Array.new(arrLen){ |i| 0*i }
  valueA = Array.new(arrLen) { |i| 0*i }
  valueT = Array.new(arrLen) { |i| 0*i }
  valueC = Array.new(arrLen) { |i| 0*i }
  valueG = Array.new(arrLen) { |i| 0*i }


  infile.each {
    |line|
    # puts line

   
    numReads = numReads + 1
    cols = line.chomp.split(' ')
    
    readName = cols[0]
    readSeq = cols[9]
    readQual = cols[10]
       
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
  }

  
  
  puts numReads 
  
  outfile.print "A \t"
  for i in 0..arrLen
    outfile.print "#{valueA[i]} \t"
  end
  outfile.print "\n T \t"
  for i in 0..arrLen
    outfile.print "#{valueT[i]} \t"
  end
  outfile.print "\n C \t"
  for i in 0..arrLen
    outfile.print "#{valueC[i]} \t"
  end
  outfile.print "\n G \t"
  for i in 0..arrLen
    outfile.print "#{valueG[i]} \t"
  end
  outfile.print "\n"

  for i in 0..arrLen
     outfile1.print "#{qual[i]} \t"
   end
  outfile1.print "\n"
  outfile1.close


  infile.close

  outfile.close

  `Rscript /ifs/scratch/c2b2/ngs_lab/xs2182/code/QCplot.R #{numReads} "#{fq}/ACTG.txt" "#{fq}/Qual.txt" "#{fq}/QC.pdf"`
  `rm "#{fq}/*.sam"`
end

main()
