## only print SNVs with MAF > 0

## example: 
##chr    pos     marker  15      16      9
#chr11   121299288       rs12576704HA-a  0;4;0;16;0      4;3;0;13;0      16;0;0;4;0

def main

  while line = ARGF.gets do 
    next if line.match("^#")
    
    cols = line.chomp.split(/\s+/)
    
    allelec = 0
    
    cols[3..-1].each do |gt|
      a = gt.split(';')
      a[1..-1].each do |n|
        allelec += n.to_i
      end
    end
    if allelec > 0
      puts line
    end
  end
end

main()
