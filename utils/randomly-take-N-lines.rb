list = ARGV[0]
num = ARGV[1].to_i

array = []
File.new(list, 'r').each do |line|
  if line=~ /^(\S+)/
    array << line.chomp
  end
end

newarray = array.sort_by {rand}

puts newarray[0,num].join("\n")

