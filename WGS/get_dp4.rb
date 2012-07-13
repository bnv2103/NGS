#!/usr/bin/env ruby
#$ -cwd

require 'getoptlong'

# At the read base column, a dot stands for a match to the reference base on the forward strand, a comma for a match on the reverse strand, `ACGTN' for a mismatch on the forward strand and `acgtn' for a mismatch on the reverse strand. A pattern `\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion between this reference position and the next reference position. The length of the insertion is given by the integer in the pattern, followed by the inserted sequence. Here is an example of 2bp insertions on three reads:


#default values
def main 
  puts calculate("C",13,"A..Tgc,g")

exit
 puts read.count "."
 puts read.count ","
exit
ref="C"
  alt="dc"
 if ref == alt
 	puts "equal"
elsif ref == alt.upcase
	puts "finally equal"
else
	puts "not equal"
end
exit

  settings = {}
  optHash = getopt()
  vcf = optHash["--vcf"]
  sample1  = optHash["--s1"]
  sample2  = optHash["--s2"]

  settings.keys.sort.each do |s|
    if optHash.key?(s)
      settings[s] = optHash[s]
    end
  end
  
  calc_dp4(vcf,sample1,sample2,settings) 
end

def  calc_dp4(vcf,sample1,sample2,settings)
TODO  dphash = {}
  o = File.new(vcf + ".dp4.vcf", 'w')
  File.new(vcf, 'r').each do |line|
    if line.match("^#") 
	o.puts line
    else
	cols=line.chomp.split(/\t/) ##remove the word chr and rais a flag if there is the word chr. so that it can be put while printing"
	if not dphash[cols[0]]  #if hash element does not exist then create it
		dphash[cols[0]] = {}
	else
		dphash[cols[0]][cols[1]] = [cols[0]..cols[-1]]
	end
    end
  end

  cols[-3] << ":DP4"
  File.new(sample1, 'r').each do |line|
	cols=line.chomp.split(/\t/) 
	if dphash[cols[0]][cols[1]]
		cols[-2] << ":"
		cols[-2] << calculate(cols[2],cols[3],cols[4])
	end
  end
  File.new(sample2, 'r').each do |line|
        cols=line.chomp.split(/\t/)
        if dphash[cols[0]][cols[1]]
                cols[-1] << ":"
		cols[-1] << calculate(cols[2],cols[3],cols[4])
        end
  end

  ##Sort and print , add "chr" if chrFlag is set
  chrarray = dphash.sort_by { |chr, pos| chr }
  chrarray.each do |chr,pos|
	posarray = pos.sort_by { |pos, line | pos }
	posarray.each do | item |
		puts o item.join("\t") ##here add the word chr if needed
	end
  end
end
def calculate(ref, depth, reads)
 arr = reads.split(//)
 fr,rr,fa,ra = 0,0,0,0
# fA, fC,fT,fG
# fa fc ft fg

 i = -1;
 while i < reads.length  do
	i+=1
	base = reads[i]
	i+=1 next if base == "^"
	next if base == "$"
	if base == "+" or base == "-"
		j=i+1
		lengthIndel=""
		while reads[j].match(/[0-9]/)
			lengthIndel << reads[j]
			j+=1
		end
		indelString=reads[j..(j + lengthIndel.to_i - 1)].join
		i=j+lengthIndel.to_i-1
		next
	end		
	fr += 1 if base == "."
	rr += 1 if base == ","
	fa += 1 if base.match(/[ACTG]/)
	ra += 1 if base.match(/[actg]/)
#	puts base
 end
 puts fr+rr+fa+ra
 puts depth
 dp4="#{fr},#{rr},#{fa},#{ra}"
# puts o "FreqA= #{(fA+fa)/depth} ; FreqC= #{(fc+fC)/depth} ;  fG+fg/depth , \n"
  return dp4
end

def getopt
  opts = GetoptLong.new(
                        ["--vcf", "-v", GetoptLong::REQUIRED_ARGUMENT],
                        ["--s1", "-1", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--s2", "-2", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--help", "-h", GetoptLong::NO_ARGUMENT]
                        )
  optHash = {}
  opts.each do |opt, arg|
    optHash[opt] = arg
  end
  if optHash.key?("--help") or !optHash.key?("--vcf") or !optHash.key?("--s1")  or !optHash.key?("--s2")
    $stderr.puts "Usage: ruby __.rb --vcf VCF --s1 file1.pileup --s2 file2.pileup"
    $stderr.puts "     options: "
    exit
  end
  return optHash
end

main()


