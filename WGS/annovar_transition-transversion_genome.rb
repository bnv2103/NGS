#!/usr/bin/env ruby

## calcualte transition/transversion ratio from genome_summary annovar file
## per sample

def main
  infile = ARGV[0]
#Class @collumn0
#"downstream"
#"exonic"
#"exonic;splicing"
#"intergenic"
#"intronic"
#"ncRNA_exonic"
#"ncRNA_intronic"
#"ncRNA_UTR3"
#"ncRNA_UTR5"
#"splicing"
#"upstream"
#"upstream;downstream"
#"UTR3"
#"UTR5"
#
#
	lastpos = -1 

	bname= File.basename(infile)
	infilearr =  bname.chomp.split(/\./)
	sid = infilearr[0].sub(" ","_")        

      	intergenic, intronic, downstream,upstream,splicing, exonic, utr3, utr5,ncrna = 0,0,0,0,0,0,0,0,0,0
	knownindel, unknownindel, knownsnv, unknownsnv, knownti, knowntv ,novelti, noveltv,onek,dbsnp  = 0,0,0,0,0,0,0,0,0,0

File.new(infile, 'r').each do |line|
    cols=line.chomp.split(/\t/)
    next if cols[0].include? "Func"
    bclass, fclass,pos,name,ref,alt,onekname =  cols[0].gsub(/"/,''),cols[2].gsub(/"/,''),cols[22].to_i,cols[8],cols[24].gsub(/"/,''),cols[25].gsub(/"/,''), cols[7]
      info,gt = cols[-3].gsub(/"/,''),  cols[-1..-1]

      ## SYN AND NONSYN
	coding, indelFlag = 0,0

	if bclass.include? "intergenic"
		intergenic+=1
	elsif bclass == "intronic"
		intronic+=1
	elsif bclass.include? "downstream"
		downstream+=1
	elsif bclass.include? "upstream"
		upstream+=1
	elsif bclass == "exonic"
		exonic+=1
 		coding =1
	elsif bclass == "UTR3"
		utr3+=1
        elsif bclass == "UTR5"
                utr5+=1
	end
	if bclass.include? "splicing"
		splicing +=1
	end
	if bclass.include? "ncRNA"
		ncrna+=1
	end 

#    if coding == 1
#	if  fclass.include? "SNV"	
 #         coding = 1
#	elsif fclass.include? "frameshift" 
#	  indelFlag  =1
#	  coding = 1
#	else 
#	  coding = 0
#	end
 #   end
     
      next if pos == lastpos       #ignore same variant i.e. same position in tht chromose
#      next if coding == 0
      if onekname != '' 
	onek+=1
	end
      if name =~ /^rs/ 
	dbsnp+=1
	end
        
   if ref.size == alt.size and (not ref.include? "-")  and (not alt.include? "-")   # it is SNV
      ti = judge(ref,alt)
      if name =~ /^rs/  # known
        if ti == 1
          knownti +=1
        else
          knowntv +=1
        end
	knownsnv +=1
      else  # novel
        if ti == 1
          novelti +=1
        else
          noveltv +=1
        end
	unknownsnv+=1
      end
   else 
      if name =~ /^rs/  # known
	knownindel+=1
	else
	unknownindel+=1
	end
   end

  lastpos = pos
end

  puts "\Samples\t#{sid}"
  print "all_variants\t#{unknownsnv+knownsnv+knownindel+unknownindel}"
  print "\n"
  
  print "all_known(dbSNP132)"
  print "\t#{dbsnp}"
  print "\n"

  print "all_novel"
  print "\t#{unknownsnv+unknownindel}"
  print "\n"

  print "all_known%(dbSNP132)"
  print "\t#{(dbsnp/(unknownsnv+knownsnv+knownindel+unknownindel).to_f).round(2)}"
  print "\n"

  print "all_1KG"
print "\t#{onek}"
  print "\n"

  print "all_1KG%"
 print "\t#{(onek/(unknownsnv+knownsnv+knownindel+unknownindel).to_f).round(2)}"
  print "\n"

  print "all_SNV\t#{knownsnv+unknownsnv}\n"
  print "all_known_SNV(dbsnp132)\t#{knownsnv}\n"
  print "all_novel_SNV\t#{unknownsnv}\n"
  print "all_INDEL\t#{knownindel+unknownindel}\n"
  print "all_known_INDEL(dbsnp132)\t#{knownindel}\n"
  print "all_unknown_INDEL\t#{unknownindel}\n"

  print "known:ti/tv"
 print "\t#{knownti}/#{knowntv}"
  print "\n"

  print "known:ti/tv-ratio"
 print "\t#{(knownti/knowntv.to_f).round(3)}"
  print "\n"


  print "novel:ti/tv"
 print "\t#{novelti}/#{noveltv}"
  print "\n"

  print "novel:ti/tv-ratio"
 print "\t#{(novelti/noveltv.to_f).round(3)}"
  print "\n"

  print "Intronic\t#{intronic}\n"
  print "Intergenic\t#{intergenic}\n"
  print "Exonic\t#{exonic}\n"
  print "Splicing\t#{splicing}\n"
  print "Downstream\t#{downstream}\n"
  print "Upstream\t#{upstream}\n"
  print "5UTR\t#{utr5}\n"
  print "3UTR\t#{utr3}\n"
  print "ncRNA\t#{ncrna}\n"

end


def judge(a1, a2)
  ti = 0
  if a1 == "A"
    if a2 == "G"  # ti
      ti = 1
    end
  elsif a1 == "C"
    if a2 == "T"
      ti = 1
    end
  elsif a1 == "T"
    if a2 == "C"
      ti = 1
    end
  elsif a1 == "G"
    if a2 == "A"
      ti = 1
    end
  end
  return ti
  
end

main()

