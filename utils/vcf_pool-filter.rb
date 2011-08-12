## filter VCF from pooled sequencing
# heuristics: 
# qual >= 20
# DP4: non-ref >=  ratio  > (singleton freq = 0.5 / N * 3), N = number of samples in each pool
# DP4: non-ref reads >= 3 in each direction

require 'getoptlong'
def main
  
###   nsample = ARGV.shift
  optHash = getopt()
  rcut, ncut = 0.2, 3
  tbc, bbc, mbc, sbc = -20, -50, -40, -20
  if optHash.key?("--ratio") 
    rcut = optHash["--ratio"].to_f
  end

  if optHash.key?("--tail")
    tbc = optHash["--tail"].to_f
  end

  if optHash.key?("--baseQ")
    bbc = optHash["--baseQ"].to_f
  end
  
  if optHash.key?("--mapQ")
    mbc = optHash["--mapQ"].to_f
  end

  if optHash.key?("--strand")
    sbc = optHash["-strand"].to_f
  end

  if optHash.key?("--nreads")
    ncut = optHash["--nreads"].to_i
  end

  File.new(optHash["--vcf"], 'r').each do |line|
    if line.match("^#")
      puts line
    else
      cols = line.chomp.split(/\s+/)
      info = cols[7]
      rflag, nflag = 1, 1
      bflag = 1

      info.split(';').each do |f|
        if f=~ /DP4\=(\d+)\,(\d+)\,(\d+),(\d+)/
          ref1, ref2, alt1, alt2 = $1.to_i, $2.to_i, $3.to_i, $4.to_i
          ratio = (alt1 + alt2 ).to_f / (ref1 + ref2 + alt1 + alt2)
          if ratio < rcut
            rflag = 0
          end
            
          if alt1 < ncut or alt2 < ncut
            nflag = 0
          end

        elsif f=~ /PV4\=(\S+)\,(\S+),(\S+),(\S+)/
          sb , bb, mb, tb = $1.to_f, $2.to_f, $3.to_f, $4.to_f
          
          if  bb == 0.0 or tb == 0.0
            bflag = 0
          elsif Math.log10(sb) < sbc or Math.log10(bb) < bbc or Math.log10(mb) < mbc or Math.log10(tb) < tbc
            bflag = 0
          end
        end
      end

      flag = []
      if rflag == 0
        flag << "RATIO"
      end
      
      if nflag == 0
        flag << "NUM-NONREF"
      end

      if bflag == 0
        flag << "BIAS"
      end

      if flag.size > 0
        flagstring = "FILTER:" + flag.join(";")
      else
        flagstring = "PASS"
      end  
      puts "#{cols[0..5].join("\t")}\t#{flagstring}\t#{cols[7..-1].join("\t")}"
    end
    
  end
end

def getopt
  
  opts = GetoptLong.new(
                        ["--vcf", "-v", GetoptLong::REQUIRED_ARGUMENT],
                        ["--ratio", "-r", GetoptLong::REQUIRED_ARGUMENT],
                        ["--tail", "-t", GetoptLong::REQUIRED_ARGUMENT],
                        ["--baseQ", "-b", GetoptLong::REQUIRED_ARGUMENT],
                        ["--mapQ", "-m", GetoptLong::REQUIRED_ARGUMENT],
                        ["--strand", "-s", GetoptLong::REQUIRED_ARGUMENT],
                        ["--nreads", "-n", GetoptLong::REQUIRED_ARGUMENT],
                        ["--help", "-h", GetoptLong::NO_ARGUMENT]
                        )
  optHash = {}
  opts.each do |opt, arg|
    optHash[opt] = arg
  end
  if optHash.key?("--help") || !optHash.key?("--vcf")
    $stderr.puts "Usage: ruby __.rb [options] -v FOO.VCF"
    
    helpinfo = <<-EOS
    options
      -v  vcf file
      -r  min ratio of non-ref reads (default: 0.2)
      -t  max p-value of tail distance bias in log10 scale (default: -20)
      -b  max p-value of baseQ bias in log10 scale (default: -50)
      -m  max p-value of mapQ bias in log10 scale (default: -40)
      -s  max p-value of strand bias in log10 scale (default: -20)
      -n  min number of reads for non-ref allele on each strand (default: 3)
    EOS
    
    $stderr.puts helpinfo    
    exit
  end
  return optHash
  
end


main()
