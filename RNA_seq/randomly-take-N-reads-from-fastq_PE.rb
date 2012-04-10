## randomly take N reads from fastq file

$VERBOSE = nil

def main
  fq = ARGV[0]
  fq2 = ARGV[1]
  num = ARGV[2]
  output = ARGV[3]
  output3 = ARGV[4]

  if fq == nil or fq2 == nil or num == nil
    $stderr.puts "Usage: ruby __.rb input.fastq number_of_reads [output.fastq]"
    exit
  end

  if output == nil or output3 == nil
    output = File.basename(fq, ".fastq") + "_rand_" + num + ".fastq"
    output3 = File.basename(fq2, ".fastq") + "_rand_" + num + ".fastq"

  end


  chunk = 2000000  # 2 million reads
  chunk = 500000

  if fq == nil or fq2 == nil or num == nil
    puts "Usage: ruby _.rb fastq NUM > random_NUM.fastq"
    exit
  end

  num = num.to_i
  total = `wc -l #{fq} | cut -f1 -d ' '`
  totalN = total.to_i / 4
  puts totalN
  io = File.new(fq, 'r')
  io3 = File.new(fq2, 'r')
  o = File.new(output, 'w')
  o3 = File.new(output3, 'w')
  array = []
  while !io.eof?
    s = io.take(4)
    s3 = io3.take(4)
    array << [s, s3]
    
    if array.size > chunk
      randomTake(array, num, totalN, o, o3)
      array = []
    end
  end
  puts array.size
  puts num
  puts totalN
  randomTake(array, num, totalN, o, o3)
  io.close
  io3.close
  o.close
  o3.close
end

def randomTake(array, num, t, o, o3)
  n = num * array.size / t
  return if  n < 1
  temp = array.sort_by {rand}
  arraylength = temp.size
  a1 = []
  a2 = []
  for ss in 0...n.to_i
    a1 << temp[ss][0]
    a2 << temp[ss][1]
  end
  o.puts a1
  o3.puts a2

end

main()
