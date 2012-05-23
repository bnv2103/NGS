# change original ID to NX0000, and make the mapping file



def main
  ped = ARGV[0]

  newpedo = File.new(ped + ".anonymized.ped", "w")
  mapo = File.new(ped + ".ID-map", "w")

  iindex = 0
  aindex = "A"

  File.new(ped, 'r').each do |line|
    if line.match(/^IndividualID/)
      newpedo.puts line
    else 
      cols = line.chomp.split(/\s+/)
    
      oid = cols[0]
      nid = "N" + aindex + "#{iindex}"
      mapo.puts "#{nid}\t#{oid}"
      newpedo.puts "#{nid} #{cols[1..-1].join(" ")}"

      iindex += 1
      if iindex > 9999
        iindex = 0
        aindex = (aindex.ord + 1).chr 
      end
    end

  end



end

main()
