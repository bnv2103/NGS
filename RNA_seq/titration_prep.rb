## 

# input: isoform or gene output from Cufflinks, under a series of number of total reads 

def main
  list = ARGV[0]

  files = {}
  File.new(list, 'r').each do |line|
    cols = line.chomp.split(/\s+/)
    files[cols[0]] = cols[1]
  end

  transcripts = {},

end

main()
