## change  gene-based positions to chr-based positions

# input: chr, shift, vcf

# SORL1: 121322961-121504471

def main
  chr = ARGV.shift
  s = ARGV.shift.to_i
  
  while line = ARGF.gets do 
    if line.match("^#")
      puts line
    else
      cols = line.chomp.split(/\t/)
      pos = cols[1].to_i + s - 1 
      puts "#{chr}\t#{pos}\t#{cols[2..-1].join("\t")}"
      
    end
  end
end

main()
