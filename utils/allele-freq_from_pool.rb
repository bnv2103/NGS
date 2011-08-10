## allele frequency from pool 

# input:  cvs file
# output:  tab-delimited file

def main
  csv = ARGV.shift
  posf = ARGV.shift

  gt = readCSV(csv)
  pos = readPos(posf)
  maf(gt, pos)
end

def readPos(posf)
  pos = []
  ## lifted over from hg18 to hg19
  File.new(posf, 'r').each do |line|
    cols = line.chomp.split(/\s+/)
    pos <<  cols[1].to_i
  end
  return pos
end

def maf(gt, pos)
  puts "\#chr\tpos\tmarker\t#{gt[0][:pools].keys.sort.join("\t")}"
  gt.each do |h| ## h: {:marker, :pools}
    m = h[:marker]
    p = pos.shift
    cmap = {'1' => 'A', '2' => 'C', '3' => 'G', '4' => 'T', '0' => 'N'}
    print "chr11\t#{p}\t#{m}"
    h[:pools].keys.sort.each do |p|
      # x = h[:pools][p].keys.sort_by {|a,b| h[:pools][p][b] <=> h[:pools][p][a] }
      print "\t#{h[:pools][p].join(";")}"
    end
    print "\n"
  end
end

def readCSV(csv)
  flag = 0
  gt = []
  markers = []
  num  = 0
  File.new(csv, 'r').each do |line|
    cols = line.chomp.split(",")
    if flag == 0 
      if cols[0].match("Pool")  # header line
        flag = 1
        cols[6..-1].each do |mk|
          markers << mk
        end
        num = markers.size  / 2
        0.upto(num-1) do |i|
          h = {:marker => markers[i*2], :pools => {} }
          gt << h
        end
      else
        next
      end
    elsif flag == 1
    
      
      pool, iid = cols[0], cols[2]
      # gt[pool] = {} unless gt.key?(pool)
      # gt[pool][iid] = []
      0.upto(num-1) do |i|
        # m = markers[i*2]
        a1 = cols[6 + i*2].to_i
        a2 = cols[7 + i*2].to_i
        if !gt[i][:pools].key?(pool)
        gt[i][:pools][pool] = [0,0,0,0,0] ## {'1' => 0, '2' => 0, '3' => 0, '4' => 0, '0' => 0}
          
        end
#        $stderr.puts "#{a1}\t#{a2}"
        gt[i][:pools][pool][a1] += 1
        gt[i][:pools][pool][a2] += 1
      end
      
    end
  end
  return gt
end

main()
