## prepare K Bio file

def main
  candf = ARGV[0]
  ref = ARGV[1]

  candidates = parse(candf)

  findContext(candidates, ref)
  
  output(candidates)
end


def output(c)
  c.each_key do |chr|
    c[chr].keys.sort.each do |pos|
      puts "#{c[chr][pos][:name]}\t#{c[chr][pos][:context]}\t#{chr}\t#{pos}\t#{c[chr][pos][:qual]}"
    end
  end
end

def process(chr, seq, c)
  c[chr].keys.sort.each do |pos|
    ref = seq[pos-1]
    c[chr][pos][:context] = seq[(pos - 51)..(pos-2)] + "[#{ref}/#{c[chr][pos][:alt]}]" + seq[(pos+1)..(pos+50)]
    if ref != c[chr][pos][:ref]
      $stderr.puts "Warning: reference allele does not match: #{chr} #{pos} #{c[chr][pos][:name]}"
    end
  end
  return c
end


def findContext(c, ref)
  flag = 0
  seq = ''
  chr = ''
  File.new(ref, 'r').each do |line|
    if line=~ /^\>(\S+)/
      if seq.size > 0 
        process(chr, seq, c)
      end
      chr = $1
      seq = ''
      $stderr.puts "load chr #{chr}"
      if c.key?(chr)
        flag = 1
      else
        flag = 0
      end
      
    elsif flag == 1
      seq <<  line.chomp
    end

  end
  if seq.size > 0
    process(chr, seq, c)
  end
  
end


def parse(f)
  c = {}
  File.new(f,'r').each do |line|
    cols = line.chomp.split(/\s+/)
    chr,pos,name,ref,alt,qual = cols[0..-1]
    c[chr] = {} if !c.key?(chr)
    c[chr][pos.to_i] = {:name => name, :ref => ref, :alt => alt, :qual => qual.to_f}
  end
  return c

end

main()
