require 'csv'

def main
  # 1. calculate the number of raw reads, mapped reads from the bam file.
 
  dir = ARGV[0]
  output = ARGV[1]

  sampleName = dir
  nreads = 0
  dirb = File.basename(dir)
  #  puts dirb
  if dirb =~ /^(\S+)\_random\_(\d+)/
    sampleName = $1
    nreads = $2
  end

  if nreads == 0
    # get number of reads
    nreads = `grep reads_in #{dir}/left_kept_reads.info | cut -f2 -d '='`.to_i
  end

   mappedline = `samtools flagstat #{dir}/accepted_hits.bam | grep mapped | head -1`
   if mappedline =~ /^(\d+)\s+/
     mapped = $1
   end
  # `samtools view #{dir}/accepted_hits.bam | cut -f3 | sort | uniq -c > #{dir}/reads`
  
    # get the sample name 
    sampleName_short = sampleName.split('_')
    
    # print sample name with statistical result
    writer = CSV.open(output, 'a') do |csv|
      csv << [sampleName_short[5], nreads, mapped]
      end
  # `Rscript /ifs/scratch/c2b2/ngs_lab/xs2182/code/printSpikeIn.R #{dir}/reads #{sampleName_short[5]}`
end


main()
