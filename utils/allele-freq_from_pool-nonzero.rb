## only print SNVs with MAF > 0

## example: 
##chr    pos     marker  15      16      9
#chr11   121299288       rs12576704HA-a  0;4;0;16;0      4;3;0;13;0      16;0;0;4;0

def main

  while line = ARGF.gets do 
    next if line.match("^#")
    
    cols = line.chomp.split(/\s+/)
    
    allelec = [0,0,0,0]
    pos = cols[1] 
    cols[3..-1].each do |gt|
      a = gt.split(';')
      1.upto(4) do |i|
        n = a[i].to_i
#      a[1..-1].each do |n|
        if n > 0
          allelec[i-1] += 1
        end
        #        allelec += n.to_i
      end
    end
    
    ga =  allelec.select {|g| g > 0}
    if ga.size > 1
      puts line
#      if pos == "121320556"
 #       $stderr.puts ga.join("\t")
 #     end
    end
  end
end

main()
