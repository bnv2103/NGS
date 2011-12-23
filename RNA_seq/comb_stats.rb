lrequire 'csv'

def main
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
    nreads = `grep out #{dir}/left_kept_reads.info | cut -f2 -d '='`.to_i
  end

  if flag == nil
    genes = "#{dir}/accepted_hits.bam_cufflinks_ref/genes.fpkm_tracking"
    isoforms = "#{dir}/accepted_hits.bam_cufflinks_ref/isoforms.fpkm_tracking"
  else
    isoforms = "#{dir}/accepted_hits.bam_cufflinks_ref/genes.fpkm_tracking"
  end

  mapped = 0
  mappedline = `samtools flagstat #{dir}/accepted_hits.bam | grep mapped | head -1`
  if mappedline =~ /^(\d+)\s+/
    mapped = $1
  end


  if File.exist?(isoforms)
    # a = `Rscript summaryR.R #{isoforms} #{genes} #{dirb}`.sub("[1]","").strip.split(/\s+/).join("\n")
    a = `Rscript /ifs/scratch/c2b2/ngs_lab/xs2182/code/summaryR.R #{isoforms} #{genes} #{dirb}`
    b = a.split(' ')
    sampleName_short = sampleName.split('_')
    # puts "#{dir}\n#{sampleName}\n#{nreads}\n#{mapped}\n#{a}"
    
    writer = CSV.open(output, 'a') do |csv|
      # csv << [#{dir},#{sampleName},#{nreads},#{mapped},#{a}]
      csv << [sampleName_short[3], mapped, b[1], b[2], b[3], b[4], b[5], b[6]]
      end
  end
end


main()