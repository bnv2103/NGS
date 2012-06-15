#!/usr/bin/env ruby
#$ -cwd

require 'getoptlong'

# At the read base column, a dot stands for a match to the reference base on the forward strand, a comma for a match on the reverse strand, `ACGTN' for a mismatch on the forward strand and `acgtn' for a mismatch on the reverse strand. A pattern `\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion between this reference position and the next reference position. The length of the insertion is given by the integer in the pattern, followed by the inserted sequence. Here is an example of 2bp insertions on three reads:

#joint.ac1-ac2.var.flt.6.vcf:22  16051968        .       C       A    182        .       DP=44;AF1=1;CLR=1;AC1=4;DP4=4,6,6,14;MQ=18;FQ=-30.6;PV4=0.69,1,1,0.45   GT:PL:DP:SP:GQ  1/1:82,0,1:15:5:6       1/1:134,18,0:15:0:24
#small.22.AC1.pileup:22  16051968        C       24      aAa,,AAAaaAaA,A,a.,.a,A.        D#JFI###DC@@-##3CF;1A#!6
#small.22.AC2.pileup:22  16051968        C       20      aAA,.aA,a.a.AaAaaaaA    !#BJ#C;GI#B;C#CB=B<F

#default values
def main
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
  dphash = {}
  cols=[]
  o = File.new(vcf + ".dp4.vcf", 'w')
  File.new(vcf, 'r').each do |line|
    if line.match("^#") 
	o.puts line
    else
	cols=line.chomp.split(/\t/) ##remove the word chr and rais a flag if there is the word chr. so that it can be put while printing"
	if not dphash[cols[0]]  #if hash element does not exist then create it
		dphash[cols[0]] = {}
	else
		cols[-3] << ":DP4"
		dphash[cols[0]][cols[1]] = cols[0..-1]
	end
    end
  end

#  puts dphash

  File.new(sample1, 'r').each do |line|
	cols=line.chomp.split(/\t/) 
	if not dphash[cols[0]][cols[1]].nil?

# puts  "#{cols[1]}	#{cols[2]}	#{dphash[cols[0]][cols[1]][4]}	#{cols[4]}	#{cols[5]} \t " << calculate(dphash[cols[0]][cols[1]][3],cols[4],cols[5])
	
		dphash[cols[0]][cols[1]][-2] << ":" << calculate(dphash[cols[0]][cols[1]][4],cols[4],cols[5])
	end
  end
  File.new(sample2, 'r').each do |line|
        cols=line.chomp.split(/\t/)
        if not dphash[cols[0]][cols[1]].nil?
#  puts  "#{cols[1]}       #{cols[2]}      #{dphash[cols[0]][cols[1]][4]}  #{cols[4]}      #{cols[5]} \t " << calculate(dphash[cols[0]][cols[1]][3],cols[4],cols[5])
                dphash[cols[0]][cols[1]][-1] << ":" << calculate(dphash[cols[0]][cols[1]][4],cols[4],cols[5])
        end
  end

  ##Sort and print , add "chr" if chrFlag is set
  chrarray = dphash.sort_by { |chr, pos| chr }
  chrarray.each do |chr,pos|
	posarray = pos.sort_by { |pos, line | pos }
	posarray.each do | item |
		o.puts item.join("\t") ##here add the word chr if needed
	end
  end
 o.close
end
def calculate(alt, reads, qual)		## Ignore bases with base quality [ascii interger-33] < 13
 arr = reads.split(//)
 qarr = qual.split(//)
 fr,rr,fa,ra = 0,0,0,0
 altarr = alt.split(/,/)
 i = -1;
 count=0
 while i < (reads.length-1)  do
	i+=1
	base = reads[i]
#	puts "#{base}	#{qarr[count]}	#{i}	#{count}\n"
	if base == "^"
		i+=1
		next
	end
	next if base == "$"
	if base == "+" or base == "-"
		j=i+1
		lengthIndel=""
		while reads[j].match(/[0-9]/)
			lengthIndel << reads[j]
			j+=1
		end
#		indelString=reads[j..(j + lengthIndel.to_i - 1)].join
		i=j+lengthIndel.to_i-1
		next
	end		
	fr += 1 if base == "."  and (qarr[count].ord - 33) >= 13
	rr += 1 if base == ","  and (qarr[count].ord - 33) >= 13
	if altarr.include?(base)  and (qarr[count].ord - 33) >= 13
		fa += 1
	elsif altarr.include?(base.upcase)  and (qarr[count].ord - 33) >= 13
		ra += 1
	end
	count+=1
 end
 dp4="#{fr},#{rr},#{fa},#{ra}"
#  puts "#{alt}\t#{reads}\t#{qual}\t#{dp4}\n"
 return dp4
end

def getopt
  opts = GetoptLong.new(
                        ["--vcf", "-v", GetoptLong::REQUIRED_ARGUMENT],
                        ["--s1", "-1", GetoptLong::REQUIRED_ARGUMENT],
                        ["--s2", "-2", GetoptLong::REQUIRED_ARGUMENT],
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


