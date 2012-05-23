input = ARGV[0]
mapf = ARGV[1]

map = {}
File.new(mapf, 'r').each do |line|
  if line=~/^\s*(\S+)\s+(\S+)/
    map[$1] = $2
  end
end

File.new(input, 'r').each do |line|
  line.chomp!
  if line=~ /^\s*(\S+)\s+(.*)/
    id, info = $1, $2
    if map.key?(id)
      id = map[id]
    end
    puts "#{id} #{info}"
  end
end

