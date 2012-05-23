def main

name1 = ARGV[0]
name2 = ARGV[1]
dir = ARGV[2]

sampleName1 = name1.split('_')
sampleName2 = name2.split('_')

if name1!=name2 &&  sampleName1[5]==sampleName2[5]
    `cat #{name1} #{name2} > #{dir}"/combine"#{name2}`
# puts name1
# puts name2
end


end

main()
