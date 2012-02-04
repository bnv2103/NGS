## check missing tiles

def main
  tilelist = ARGV[0]

  bclDir = ARGV[1]

  tiles = readList(tilelist)
  
  check(bclDir, tiles)

end

def check(dir, tiles)
  target = dir + "/Data/Intensities/BaseCalls/"
  
  tiles.each_key do |lane| 
    1.upto(209) do |i|
      circle = target + "/L00#{lane}/C#{i}.1/" 
      tiles[lane].each do |tile|
        bcl = circle + "s_#{lane}_#{tile}.bcl"
        stats = circle + "s_#{lane}_#{tile}.stats"
        
        if File.exist?(bcl) and File.exist?(stats)
          $stderr.puts "good: #{lane}\t#{tile}\t#{i}"
        else
          puts "#{lane}\t#{tile}\t#{i}"
        end
      end
    end
  end
end


def readList(l)
  tiles = {}
  File.new(l,'r').each do |line|
    cols = line.split(/\s+/)
    cols.each do |t|
      if t=~ /^s\_(\d)\_(\d+)/
        lane , tile = $1, $2

        if !tiles.key?(lane)
          tiles[lane] = []
        end

        tiles[lane] << tile

      end

      

    end
  end

  return tiles
end

main
