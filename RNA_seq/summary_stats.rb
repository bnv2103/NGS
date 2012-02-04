#


def main
  dir = ARGV[0]
  flag = ARGV[1]
 
  sampleName = dir
  nreads = 0
  dirb = File.basename(dir)
  if dirb =~ /^(\S+)\_random\_(\d+)/
    sampleName = $1
    nreads = $2
  end

  if nreads == 0
    nreads = `grep out #{dir}/left_kept_reads.info | cut -f2 -d '='`.to_i
  end

  if flag == nil

    isoforms = "#{dir}/accepted_hits.bam_cufflinks_ref/isoforms.fpkm_tracking"
  else
    isoforms = "#{dir}/accepted_hits.bam_cufflinks_ref/genes.fpkm_tracking"
  end
  
  mapped = 0
  mappedline = `samtools flagstat #{dir}/accepted_hits.bam | grep mapped | grep -v different`
  if mappedline =~ /^(\d+)\s+/
    mapped = $1
  end

  puts isoforms
  if File.exist?(isoforms)
    a = `Rscript ~/code/NGS/RNA_seq/cufflinks.summary.R #{isoforms}`.sub("[1]","").strip.split(/\s+/).join("\t")
    
    puts "#{dir}\t#{sampleName}\t#{nreads}\t#{mapped}\t#{a}"
  end
end


main()
