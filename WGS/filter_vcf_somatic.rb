#!/usr/bin/env ruby

require 'getoptlong'


### filter VCF based on:
# 1.DP4 >= 2,2 
# 2. PV4>= 0.001 for het 
# 3. CLR>20 for somatic 

#default values
def main 
  settings = {}
  settings["--dp4"]=2
  settings["--pv41"]=0.0001
  settings["--pv42"]=0.0001
  settings["--pv43"]=0.0001
  settings["--pv44"]=0.0001
  settings["--clr"]=20
  settings["--normal"]=1
 
  optHash = getopt()
  vcf = optHash["--vcf"]
  
  settings.keys.sort.each do |s|
    if optHash.key?(s)
      settings[s] = optHash[s].to_f
    end
  end
  
  nsample=countSamples(vcf)
  
  filterVCF(vcf,settings,nsample)  # gt: gene -> pos -> sample -> genotype, 

end

def filterVCF(vcf, settings, nsample)
  o = File.new(vcf + ".somaticfiltered.vcf", 'w')
  e = File.new(vcf + ".dropped.vcf", 'w')
  File.new(vcf, 'r').each do |line|
    if line.match("^#") 
      o.puts line
      e.puts line
    else
      cols=line.chomp.split(/\t/)
      qual, info = cols[5].to_f, cols[7].split(';')
#10      61852   .       T       C       999     .       DP=76;AF1=0.5;AC1=2;DP4=9,17,11,15;MQ=40;FQ=999;PV4=0.78,0.17,0.42,0.043        GT:PL:DP:SP:GQ  0/1:215,0,221:28:4:99   0/1:219,0,208:24:0:99
#10      68575   .       C       T       999     .       DP=72;AF1=1;AC1=4;DP4=4,6,24,25;MQ=29;FQ=-60.5;PV4=0.73,1,1,1   GT:PL:DP:SP:GQ  1/1:235,29,0:31:4:63    1/1:239,41,0:28:2:75

	format = cols[8].split(':')
	i = 0
	format.each do |gt|
		if gt == "GT"
			break
		end
		i+=1
	end
	if settings["--normal"] == 1
		gt_norm,gt_tumor =  cols[9].split(':')[i] ,cols[10].split(':')[i]
	else
                gt_tumor,gt_norm =  cols[9].split(':')[i] ,cols[10].split(':')[i]
	end
	
	clrflag = 0     
      	flag = 0 
	if gt_norm == "0/0" and  gt_tumor == "0/1"
		flag = 0
	else
		flag = 1
        end
	if flag == 0
        info.each do |item|
          k,v=item.split('=')[0..1]
          if k == "DP4" 
	    v1=v.split(',')
            if v1[2].to_f < settings["--dp4"] or v1[3].to_f < settings["--dp4"] 
              flag = 1
            end
          elsif k == "PV4" 
	        v1=v.split(',')
            	if v1[0].to_f < settings["--pv41"] or v1[1].to_f < settings["--pv42"] or v1[2].to_f < settings["--pv43"] or v1[3].to_f < settings["--pv44"]
			flag = 1
		end
          elsif k == "CLR" 
		clrflag = 1
		if v.to_f < settings["--clr"]
	        	flag = 1
		end
          end
	end
	end
      if flag == 0 and clrflag == 1
        o.puts line
      else
        e.puts line 
      end
    end
  end
  o.close
  e.close
end

def countSamples(vcf)
  n = 0
  File.new(vcf, 'r').each do |line|
    n += 1
    if line.match("^#CHROM") 
      cols = line.chomp.split(/\s+/)
      return cols.size - 9 
    elsif n > 1000
      return 0
    end
  end
end

def getopt
  
  opts = GetoptLong.new(
                        ["--vcf", "-v", GetoptLong::REQUIRED_ARGUMENT],
                        ["--dp4", "-d", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--pv41", "-p", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--pv42", "-q", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--pv43", "-r", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--pv44", "-s", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--clr", "-c", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--normal", "-n", GetoptLong::REQUIRED_ARGUMENT],
                        ["--help", "-h", GetoptLong::NO_ARGUMENT]
                        )
  optHash = {}
  opts.each do |opt, arg|
    optHash[opt] = arg
  end
  if optHash.key?("--help") or !optHash.key?("--vcf") or !optHash.key?("--normal")
    $stderr.puts "Usage: ruby __.rb --vcf VCF --normal [1/2] [--dp4[2] --pv41[0.0001] --pv42[0.0001] --pv43[0.0001] --pv44[0.0001] --clr[20] ]"
    $stderr.puts "     options: "
    exit
  end
  return optHash
  
end


main()


