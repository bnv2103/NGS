
def main

dir = ARGV[0]
o = File.new("transcripts.txt", 'a+')
genes = "#{dir}/accepted_hits.bam_cufflinks_ref/transcripts.gtf"
puts genes
o.puts genes
o.close

end


main()
