
def main
  genome = ARGV[0]
  group1 = ARGV[1]
  group2 = ARGV[2]
  dir = ARGV[3]
  output = ARGV[4]

  group1_samples = group1.split(',')
  group2_samples = group2.split(',')
  
  aFile = File.new("#{output}transcripts.txt", "w")
  bamFile=File.new("#{output}bams.txt","w")
  i = 0
  group1Len = group1_samples.length
  
  while i < group1Len  do
    file1 = `echo #{dir}/*_#{group1_samples[i]}_*`.strip!
    temp = "#{file1}/accepted_hits.bam_cufflinks_ref/transcripts.gtf"
    aFile.puts(temp)
    bams = "#{file1}/accepted_hits.bam"
    if i == 0
      bamFile.write(bams)
    else
      bamFile.write(","+bams)
    end
    puts(file1)
   i = i + 1
  end
  bamFile.write(" ")
  i = 0
  group2Len = group2_samples.length
  while i < group2Len  do
    file1 = `echo #{dir}/*_#{group2_samples[i]}_*`.strip!
    temp = "#{file1}/accepted_hits.bam_cufflinks_ref/transcripts.gtf"
    bams = "#{file1}/accepted_hits.bam"
    aFile.puts(temp)
    if i == 0
      bamFile.write(bams)
    else
      bamFile.write(","+bams)
    end

    puts(file1)
   i = i + 1
  end
  bamFile.close
  aFile.close

end
  

main()
