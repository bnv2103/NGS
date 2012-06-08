## randomly take first N reads from fastq file

$VERBOSE = nil

def main
  fq = ARGV[0]
  num = ARGV[1]
  output = ARGV[2]

  
  num = num.to_i
  
  io = File.new(fq, 'r')
  
  o = File.new(output, 'w')
  
  array = []
  i = 0
  while i < num
    s = io.take(1)
    # o.puts s
    s = io.take(1)
    temp = s.to_s
    # puts temp.length
    temp = temp[2..101]
    o.puts temp
    s = io.take(2)
    i = i + 1
  end
  
  io.close
  
  o.close
  
end


main()
