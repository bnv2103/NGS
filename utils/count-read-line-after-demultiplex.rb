
def main

  n = 0 
  while line=ARGF.gets do 
    cols = line.chomp.split(/\s+/)
    n += cols[1].to_i
    
  end
  puts n 
end

main()
