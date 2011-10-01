## sort numerically (in the presence of 1e-05)

def main
  file=ARGV[0]
  column = ARGV[1]

  if column != nil
    c = column.to_i - 1
  else
    c = 0
  end

#  $stderr.puts c

  header = ""
  lines = {}
  File.new(file,'r').each do |line|
    cols = line.chomp.split(/\t/)
    lines[line.chomp] = cols[c].to_f
  end
  
  x = lines.keys.sort_by {|a,b| lines[a] <=> lines[b] }
  puts x.join("\n")
end

main()
