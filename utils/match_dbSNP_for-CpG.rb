## purpose:
# 1. assign dbSNP ID based on chr and pos
# 2. find additional SNPs within 10bp window of current SNP


def main
  list = ARGV[0]
  dbSNP = ARGV[1]

  $snpH = {}
  $snpA = []
  readlist(list)
  
#  $stderr.puts $snpH.keys
  scan(dbSNP)
  output()
end


def output
  $snpA.each do |snp|
    puts "#{snp[:id]}\t#{snp[:chr]}\t#{snp[:pos]}\t#{snp[:dbsnp]}\t#{snp[:window11].join(";")}\t#{snp[:window21].join(";")}"
    
  end
end

def sortbychr(h, chr)
  
  a = []
  
  if h.key?(chr)
    a = h[chr].keys.sort
    
  end
  return a
end

def scan(db)
  
  
  currentChr = "0"
  
  posArr  = sortbychr($snpH, currentChr)
  # db is sorted by chr
  File.new(db, 'r').each do |line|
    cols = line.split(/\s+/)
    chr, pos, name = cols[1], cols[3].to_i, cols[4]
#    $stderr.puts pos
    if chr=~ /^chr(\S+)/
      chr= $1
      if currentChr != chr
        currentChr = chr
        posArr  = sortbychr($snpH, currentChr)
        $stderr.puts "chr#{chr}\t#{posArr.size}"
      end
    end
 
  
#    $stderr.puts "#{pos}\t\t#{posArr.join("; ")}"    
    while posArr.size > 0   and pos >=  posArr[0] + 10
      posArr.shift
    end
#    $stderr.puts "#{pos}\t\t#{posArr.join("; ")}\n\n"    

    next if posArr.size < 1

    
    
    j = 0
    while 1
      i = posArr[j]
      if pos == i # match
        $snpH[chr][i][:dbsnp] = name
        # $stderr.puts "match: #{name}\t#{chr}:#{pos}"
      elsif pos > i - 6  and pos < i + 6 # within 11bp window
        $snpH[chr][i][:window11] << name
        #     $stderr.puts "window: #{name}\t#{chr}:#{pos}"
      elsif pos > i - 10
        $snpH[chr][i][:window21] << name
      end
      j = j + 1
      break if j >= posArr.size or  pos < posArr[j] - 9
    end
  end
  return 
end



def readlist(list)
  

  File.new(list, 'r').each do |line|
    cols=line.split(/\s+/)
   # $stderr.puts "#{cols.size}"
    next if cols.size != 3
    
    id, chr, pos = cols[0], cols[1], cols[2].strip.to_i
    s = {}
    
    s[:id],s[:chr], s[:pos] = id, chr, pos
    s[:dbsnp] = ""
    s[:window11] = []
    s[:window21] = []
   # $stderr.puts s
    $snpA << s
    if !$snpH.key?(cols[1])
      $snpH[chr] = {}
    end

    $snpH[chr][pos] = s
    
  end
  return
end


main()
