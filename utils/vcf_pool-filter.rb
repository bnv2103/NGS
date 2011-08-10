## filter VCF from pooled sequencing
# heuristics: 
# qual >= 20
# DP4: non-ref >=  ratio  > (singleton freq = 0.5 / N * 3), N = number of samples in each pool
# DP4: non-ref reads >= 3 in each direction


def main
  
  nsample = ARGV.shift
  
  rcut = 0.25 / nsample.to_f
  ccut = 3

  while line = ARGF.gets do 
    if line.match("^#")
      puts line
    else
    
      cols = line.chomp.split(/\s+/)
      
      info = cols[7]
      
      info.split(';').each do |f|
        if f=~ /DP4\=(\d+)\,(\d+)\,(\d+),(\d+)/
          ref1, ref2, alt1, alt2 = $1.to_i, $2.to_i, $3.to_i, $4.to_i
          ratio = (alt1 + alt2 + 0.0001) / (ref1 + ref2 + alt1 + alt2 + 0.0001)
          if ratio >= rcut and alt1 >= ccut and alt2 >= ccut
            puts "#{cols[0..5].join("\t")}\tPASS\t#{cols[7..-1].join("\t")}"
          else
            puts "#{cols[0..5].join("\t")}\tFILTER\t#{cols[7..-1].join("\t")}"
          end
      end
        
      end
    end
  end
end
main()
