def main

  
  snm = Array.new(300, 0)
  ins = Array.new(300, 0)
  del = Array.new(300, 0)
  nr = Array.new(300,0)

  while line=ARGF.gets
    cols = line.split(/\s+/)
    flag, cigar, s, md = cols[1].to_i, cols[5], cols[9].size, cols[13]
    
    ## next if flag != 0

    nr[s] += 1

    next if md == nil or md == ""
    mda = md.scan(/([0-9]+)([A-Z]|\^[A-Z]+)/)
    index = -1
    mda.each do |error|
      index += error[0].to_i + 1 
      if error[1] =~ /\^/  # insertion
        ins[index] += 1
      else
        snm[index] += 1
      end
    end
    
    index = 0
    delsA = cigar.scan(/(\S+)D/)
    
    delsA.each do |dd|
      x = dd[0].scan(/\d+/)
      x.each do |y|
        index += y.to_i
      end
    end
    
    if index > 0
      del[index] += 1
    end
    
  end

  nreads = Array.new(300, 0)
  0.upto(299) do |i|
    i.upto(299) do |j|
      nreads[i] += nr[j]
    end
  end
    
  0.upto(299) do |i|
    puts "#{i}\t#{nreads[i]}\t#{snm[i]}\t#{del[i]}\t#{ins[i]}"
  end
 
end


main()
