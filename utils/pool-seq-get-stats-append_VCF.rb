# example :

#chr11   121361015       .       G       A       -0      PASS    DP=13335;AF1=0;AC1=0;DP4=5724,6804,356,451;MQ=37;FQ=-212;PV4=0.4,0.029,0.00054,0.098    GT:PL:GQ        0/0:0,255,255:99        0/0:0,193,255:99

def main
#  puts "pos\ttitv\tqual\tratio\tref1\tref2\talt1\talt2\tsb\tqb\tmb\ttb"

  while line = ARGF.gets do
    if line.match("^#CHR") # header
      puts line.chomp + "\tfunctionalClass\tTiTv\tRatio\tref1\tref2\talt1\talt2\tstrand-bias\tbaseQ-bias\tmapQ-bias\ttail-bias"
      next
    end

    next if line.match(/^#/)
    cols = line.split(/\s+/)
    chr, pos , ref, alt, qual, info = cols[0], cols[1].to_i, cols[3], cols[4], cols[5].to_f, cols[7]

    ti = judgeTiTv(ref, alt)
    ref1, ref2, alt1, alt2 = 0, 0, 0, 0
    sb , qb, mb, tb  = 0.0, 0.0, 0.0, 0.0
    ratio = 0
    fclass = ""

    info.split(';').each do |f|
      if f=~ /DP4\=(\d+)\,(\d+)\,(\d+),(\d+)/
        ref1, ref2, alt1, alt2 = $1.to_i, $2.to_i, $3.to_i, $4.to_i
        ratio = (alt1 + alt2).to_f / (ref1 + ref2 + alt1 + alt2)
      elsif f=~ /PV4\=(\S+)\,(\S+),(\S+),(\S+)/
        sb , qb, mb, tb = $1.to_f, $2.to_f, $3.to_f, $4.to_f
      elsif f=~ /functionalClass\=(\S+)/
        fclass = $1
      end
    end

    if sb < 1e-300 
      sb = 1e-300 
    end
    
    if qb < 1e-300
      qb = 1e-300
    end
    
    if mb < 1e-300
      mb = 1e-300
    end

    if tb < 1e-300
      tb = 1e-300
    end

    puts "#{line.chomp}\t#{fclass}\t#{ti}\t#{ratio.round(5)}\t#{ref1}\t#{ref2}\t#{alt1}\t#{alt2}\t#{sb}\t#{qb}\t#{mb}\t#{tb}"
    
  end
end


def judgeTiTv(a1, a2)
  ti = 0
  if a1 == 'A'
    if a2 == 'G'  # ti
      ti = 1
    end
  elsif a1 == 'C'
    if a2 == 'T'
      ti = 1
    end
  elsif a1 == 'T'
    if a2 == 'C'
      ti = 1
    end
  elsif a1 == 'G'
    if a2 == 'A'
      ti = 1
    end
  end
  return ti
  
end



main()
