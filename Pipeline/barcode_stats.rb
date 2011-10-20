fq = ARGV[0]

io = File.new(fq, 'r')

bcode = {}
while !io.eof?
  bcl = (io.take(4))[1].chomp
  bc = bcl[0..6]
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
