require 'csv'

def main
  output = ARGV[0]
  writer = CSV.open(output, 'w') do |csv|
    csv << ["Sample Name", "Number of Mapped Reads", "Median FPKM (isoforms)", "Mean FPKM (isoforms)", "Num of transcripts (FPKM > 1)", "Num of transcripts (FPKM > 0.1)", "Num of Genes (FPKM > 1)", "Num of Genes (FPKM > 0.1)"] 
  end

end

main()
