#!/bin/env ruby

colN = ARGV.shift

if colN == nil
    num = 0
else
    num = colN.to_i
end

total = 0.0
while line = ARGF.gets do 
    cols = line.chomp.split(/\s+/)
    total += cols[num].to_f
end

puts total
