require 'csv'
$VERBOSE = nil

def main
  fq = ARGV[0]
  outfile = ARGV[1]
  infile = File.new(fq, "r")
  outfile = File.new(outfile, "w")

  reads = Hash.new
  infile.each {
    |line|
    cols = line.chomp.split(' ')
    readName = cols[0]
    chrNum = cols[2]
    startPos = cols[3]
    endPos = cols[7]
    if (endPos.to_i > 0 && startPos.to_i > 0)
      # if (startPos.to_i < endPos.to_i)
        for i in startPos.to_i .. endPos.to_i
          outfile.write(i.to_s + "\n")
        end
      # end
      # outfile.write(startPos + "\t" + endPos + "\n")
    end
  }

  infile.close()
  outfile.close()

end

main()
