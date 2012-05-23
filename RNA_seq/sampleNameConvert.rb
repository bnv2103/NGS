
def main
  # input = ARGV[0]
  # output = ARGV[1]
  infile = File.new("sampleName.txt", "r")
  outfile = File.new("junk.jnk", "w")
  infile.each {
    |i|
    outfile.write i
  }
outfile.close()
infile.close()


end
