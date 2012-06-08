#!/usr/bin/env ruby

require 'getoptlong'

### filter annovar based on:
# 1.DP4 >= 2,2 
# 2. PV4>= 0.001 for het 

#default values
def main 
  settings = {}
  settings["--dp4"]=2
  settings["--pv41"]=0.0001
  settings["--pv42"]=0.0001
  settings["--pv43"]=0.0001
  settings["--pv44"]=0.0001
 
  optHash = getopt()
  vcf = optHash["--annovar"]
  
  settings.keys.sort.each do |s|
    if optHash.key?(s)
      settings[s] = optHash[s].to_f
    end
  end
  
#  nsample=countSamples(vcf)
  
  filterVCF(vcf,settings)  # gt: gene -> pos -> sample -> genotype, 

end

def filterVCF(vcf, settings)
  o = File.new(vcf + ".filtered", 'w')	#variants passing filter
  e = File.new(vcf + ".dropped", 'w')   #variants NOT passing filter criteria
  File.new(vcf, 'r').each do |line|
    if line.match("^Func") 
      o.puts line
      e.puts line
    else
      cols=line.chomp.split(/\t/)
      func, fclass, info = cols[0].gsub(/"/,''), cols[2].gsub(/"/,''), cols[33].gsub(/"/,'').split(';') 
      
	flag = 0 
        if func.include? "splic"	##exclude splicing variants
		flag = 1
	else
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
          end
         end
	 end
      if flag == 0 
        o.puts line
      else
        e.puts line 
      end
    end
  end
  o.close
end
def getopt
  
  opts = GetoptLong.new(
                        ["--annovar", "-v", GetoptLong::REQUIRED_ARGUMENT],
                        ["--dp4", "-d", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--pv41", "-p", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--pv42", "-q", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--pv43", "-r", GetoptLong::OPTIONAL_ARGUMENT],
                        ["--pv44", "-s", GetoptLong::OPTIONAL_ARGUMENT],
                         ["--help", "-h", GetoptLong::NO_ARGUMENT]
                        )
  optHash = {}
  opts.each do |opt, arg|
    optHash[opt] = arg
  end
  if optHash.key?("--help") or !optHash.key?("--annovar") 
    $stderr.puts "Usage: ruby __.rb -v ANNOVAR.tsv [options]"
    $stderr.puts "     options: "
    exit
  end
  return optHash
  
end


main()

exit
