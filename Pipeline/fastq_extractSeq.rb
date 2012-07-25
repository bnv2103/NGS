require 'csv'

$VERBOSE = nil

def main
  seq = ARGV[0] #pdf file
  fq = ARGV[1] # sampl.fastq
  fq3 = ARGV[2] # sampl.fastq3 if there is 

  seq = "CTGTAGGCACCATCAATAGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG"

  extract = File.new("#{fq}_extracted.fastq", "w")
  remaining = File.new("#{fq}_remaining.fastq", "w")

  numReads = 0
  numReadsExtract = 0
  if fq3.nil?		##TODO Also deal with the paired-end read 
	pair = 1
  else
	if File.exists?(fq3)
		pair =2
	else
		pair =1
	end
  end
 
  infile = File.new(fq, "r")

while !infile.eof?
	s = infile.take(1)
	s1 = infile.take(1)
	s2 = infile.take(1)
	s3 = infile.take(1)
	if s1.to_s.include? seq 	# "CTGTAGGCACCATCAATAGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG"
		extract.puts s
		extract.puts s1
		extract.puts s2
		extract.puts s3
		numReadsExtract+=1
	else
		remaining.puts s
		remaining.puts s1
		remaining.puts s2
		remaining.puts s3
	end
	numReads+=1
end

infile.close
extract.close
remaining.close

puts "Total # of reads = #{numReads}"
puts "Total # of reads extracted = #{numReadsExtract}"
puts "Sequence of interest = #{seq}"

end

main()
