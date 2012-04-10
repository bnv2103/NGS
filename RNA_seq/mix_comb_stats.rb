require 'csv'

def main
  dir = ARGV[0] # Human data
  dir2 = ARGV[1] # Mouse data
  output = ARGV[2]
  flag = ARGV[3]

  sampleName = dir
  sampleName2 = dir2
  nreads = 0
  dirb = File.basename(dir)
  dirb2 = File.basename(dir2)
  #  puts dirb
  if dirb =~ /^(\S+)\_random\_(\d+)/
    sampleName = $1
    nreads = $2
  end

  if nreads == 0
    nreads = `grep reads_in #{dir}/left_kept_reads.info | cut -f2 -d '='`.to_i
    nreads2 = `grep reads_in #{dir2}/left_kept_reads.info | cut -f2 -d '='`.to_i
  end

  if flag == nil
    genes = "#{dir}/accepted_hits.bam_cufflinks_ref/genes.fpkm_tracking"
    isoforms = "#{dir}/accepted_hits.bam_cufflinks_ref/isoforms.fpkm_tracking"
    genes2 = "#{dir2}/accepted_hits.bam_cufflinks_ref/genes.fpkm_tracking"
    isoforms2 = "#{dir2}/accepted_hits.bam_cufflinks_ref/isoforms.fpkm_tracking"
 
  end
  overlaps = 0
  uqMapped = 0
  uqMapped2 = 0
  #uqMapped = `samtools view #{dir}/accepted_hits.bam | cut -f1 | sort -u -S 4000M | wc -l`.to_i
  #uqMapped2 = `samtools view #{dir2}/accepted_hits.bam | cut -f1 | sort -u -S 4000M | wc -l`.to_i  
  
  o = File.new("readsName3.txt", 'w')
  # rawLines = `samtools view #{dir}/accepted_hits.bam | cut -f1 | sort -u -S 2G`
  o.puts `samtools view #{dir}/accepted_hits.bam | cut -f1 | sort -u -S 4G`
  o.close
  
  o = File.new("readsName4.txt", 'w')
  #rawLines = `samtools view #{dir2}/accepted_hits.bam | cut -f1 | sort -u -S 2G`
  o.puts `samtools view #{dir2}/accepted_hits.bam | cut -f1 | sort -u -S 4G`
  o.close

  # overlaps = `grep -f readsName1.txt readsName2.txt | wc -l`.to_i #does not work if txt file is large
  overlaps = `comm -12 readsName3.txt readsName4.txt| wc -l`.to_i
  
  mapped = 0;
  mappedline = `samtools flagstat #{dir}/accepted_hits.bam | grep mapped | head -1`
  if mappedline =~ /^(\d+)\s+/
    mapped = $1
  end
  mapped2 = 0;
  mappedline = `samtools flagstat #{dir2}/accepted_hits.bam | grep mapped | head -1`
  if mappedline =~ /^(\d+)\s+/
    mapped2 = $1
  end



  if File.exist?(isoforms)
    a = `Rscript /ifs/scratch/c2b2/ngs_lab/xs2182/code/summaryR.R #{isoforms} #{genes} #{dirb}`
    b = a.split(' ')
    
    a2 = `Rscript /ifs/scratch/c2b2/ngs_lab/xs2182/code/summaryR.R #{isoforms2} #{genes2} #{dirb2}`
    b2 = a2.split(' ')

    sampleName_short = sampleName.split('_')

    writer = CSV.open(output, 'a') do |csv|
      csv << [sampleName_short[3], overlaps, nreads, uqMapped, mapped, b[1], b[2], b[3], b[4], b[5], b[6]]
      csv << [sampleName_short[3], overlaps, nreads2, uqMapped2, mapped2, b2[1], b2[2], b2[3], b2[4], b2[5], b2[6]]
      end
  end
end


main()
