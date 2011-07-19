fq = ARGV[0]

io = File.new(fq, 'r')

bcode = {}
while !io.eof?
  bc = (io.take(4))[1].chomp
  if !bcode.key?(bc)
    bcode[bc] = 0
  end
  bcode[bc] += 1
end

bcode.keys.sort {|a,b| bcode[b] <=> bcode[a] }.each  do |k|
  v = bcode[k]
  puts "#{k}\t#{v}"
end

io.close
