## count the number of reads from fastq file
require 'csv'

$VERBOSE = nil

def main
  fq = ARGV[0]
  output = ARGV[1]


  io = File.new(fq, 'r')

  o = File.new(output, 'a+')
  
  num = 0
  while !io.eof?
    s = io.take(4)
    num = num + 1
    
  end
  puts fq
  puts num

  # writer = CSV.open(output, 'a+') do |csv|
  #    csv << [fq, num]
  #     end

  o.puts fq + " #{num}" 
  # o.puts num
  io.close

   o.close

end


main()
