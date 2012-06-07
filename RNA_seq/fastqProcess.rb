#!/usr/bin/ruby

def main()

name1 = ARGV[0]
infile1 = File.new(name1, "r")
outfile = File.new("chr.sizes", "w")
count = 0
outfile1 = File.new("test", "w")
name = " "
infile1.each {
    |line|
    
    if line.match(">")
        if count > 0          
          outfile.write name + "\t" + count.to_s() + "\n"
          # outfile1 = File.new(name, "w")
          # outfile1.write ">" + name
          
          outfile1.close()
        end
      
      cols = line.chomp.split(">")
      name = cols[1]
      puts count
      puts cols[1]
      fileName = name + ".fa"
      outfile1 = File.new(fileName, "w")
      outfile1.write ">" + name + "\n"


      count = 0
    else
      count = count + line.chomp.length
      outfile1.write line.chomp + "\n"
      # puts line      
    end
}
outfile.write name + "\t" + count.to_s() + "\n"


infile1.close()
outfile.close()


end


main()
