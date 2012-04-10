require 'csv'

def main
  # 1. calculate the number of raw reads, unique reads from the bam file.
  # 2. call summaryR.R to get the basic RPKM information of isoforms/genes 
  dir = ARGV[0]
  output = ARGV[1]
  flag = ARGV[2]

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

  if flag == nil
    genes = "#{dir}/accepted_hits.bam_cufflinks_ref/genes.fpkm_tracking"
    isoforms = "#{dir}/accepted_hits.bam_cufflinks_ref/isoforms.fpkm_tracking"
  else
    isoforms = "#{dir}/accepted_hits.bam_cufflinks_ref/genes.fpkm_tracking"
  end

  # get number of unique reads
  uqMapped = 0
  uqMapped = `samtools view #{dir}/accepted_hits.bam | cut -f1 | sort -u -S 10G | wc -l`.to_i
  continueMap = `samtools view #{dir}/accepted_hits.bam |cut -f6 |grep "101M" | wc -l`.to_i

  if File.exist?(isoforms)
    a = `Rscript /ifs/scratch/c2b2/ngs_lab/xs2182/code/summaryR.R #{isoforms} #{genes} #{dirb}`
    b = a.split(' ')
    # get the sample name 
    sampleName_short = sampleName.split('_')
    
    # print sample name with statistical result
    writer = CSV.open(output, 'a') do |csv|
      csv << [sampleName_short[5], nreads, uqMapped, continueMap, b[1], b[2], b[3], b[4], b[5], b[6]]
      end
  end
end


main()
