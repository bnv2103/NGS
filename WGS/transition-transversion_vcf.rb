#!/usr/bin/env ruby

## calcualte transition/transversion ratio from genome_summary VCF file for 1 sample

def main
 infile = ARGV[0]
#Class 
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
##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  AC2
#21      9411193 .       N       G       90      .       AC1=2;AF1=1;DP=17;DP4=0,0,16,0;FQ=-75;MQ=26     GT:PL:GQ        1/1:123,48,0:93
#21      9411327 .       C       G       225     .       AB=0.663;AC1=1;AF1=0.5;BaseQRankSum=-0.159;DP=89;DP4=29,23,14,11;Dels=0.00;FQ=225;FS=0.000;GC=44.55;HRun=1;HaplotypeScore=13.8593;LowMQ=0.0000,0.0112,89;MQ=46.10;MQ0=0;MQ0Fraction=0.0000;MQRankSum=-0.395;PV4=1,0.34,0.31,1;QD=2.53;ReadPosRankSum=0.023   GT:AB:AD:DP:FA:GQ:MQ0:PL        0/1:0.66:59,30:89:0.337:99:0:255,0,255

lastpos = -1 

bname= File.basename(infile)
homo, het, coding, intergenic, intronic, downstream,upstream,splicing, exonic, utr3, utr5,ncrna = 0,0,0,0,0,0,0,0,0,0,0,0
indelKnown, indelNovel, indel1KG, indelDBSNP, indelNOTDBSNP ,indelESP =0,0,0,0,0,0
snvKnown,snvNovel, knownti, knowntv ,novelti, noveltv, snvDBSNP,snvNOTDBSNP,snv1KG,snvESP  = 0,0,0,0,0,0,0,0,0,0
sid=""

File.new(infile, 'r').each do |line|
    cols=line.chomp.split(/\t/)
    next if ( cols[0] =~ /##/) == 0
    if ( cols[0].downcase =~ /#chr/ ) == 0
	sid = cols[9]
	next
    end

#	bclass, fclass,pos,name,ref,alt,onekname =  cols[0].gsub(/"/,''),cols[2].gsub(/"/,''),cols[22].to_i,cols[8],cols[24].gsub(/"/,''),cols[25].gsub(/"/,''), cols[7]
		
	pos,name,ref,alt,infoall,genotype = cols[1],cols[2],cols[3],cols[4],cols[7],cols[9].split(":")[0]
	info = infoall.split(/;/)
	
	next if pos == lastpos       #ignore same variant i.e. same position in tht chromose

   if ref.size == alt.size and (not ref.include? "-")  and (not alt.include? "-")   # it is SNV
      	ti = judge(ref,alt)
	if not  infoall.include? "dbSNP" and not  infoall.include? "1KG" and not infoall.include? "ESP"
		snvNovel += 1
	else
		snvKnown += 1
	end
      	if infoall.include? "dbSNP"  # knownDBSNP ti/tv
        	if ti == 1
	          knownti +=1
        	else
	          knowntv +=1
        	end
		snvDBSNP +=1
	else  # novel ti/tv
        	if ti == 1
	          novelti +=1
        	else
	          noveltv +=1
        	end
		snvNOTDBSNP +=1
	end
	if infoall.include? "1KG"
		snv1KG +=1 
	end
	if infoall.include? "ESP"
		snvESP +=1 
	end
      ## SYN AND NONSYN
	if infoall.include? "intergenic"
		intergenic+=1	
	end
	if infoall.include? "intronic"
		intronic+=1
	end
	if infoall.include? "downstream"
		downstream+=1
	end
	if infoall.include? "upstream"
		upstream+=1
	end
	if infoall.include? "exonic"
		exonic+=1
 		coding += 1
	end
	if infoall.include? "UTR3"
		utr3+=1
        end
	if infoall.include? "UTR5"
                utr5+=1
	end
	if infoall.include? "splicing"
		splicing +=1
	end
	if infoall.include? "ncRNA"
		ncrna+=1
	end 
        if genotype == '0/1'  ## het
          	het += 1
        elsif genotype == '1/1' ## homo
	        homo += 1
        end
	##SNV end
   else		#INDEL  begin
        if not  infoall.include? "dbSNP" and not  infoall.include? "1KG" and not infoall.include? "ESP"
                indelNovel += 1
	else
		indelKnown += 1
        end
	if infoall.include? "1KG"
		indel1KG += 1 
	end	
        if infoall.include? "dbSNP"  
		indelDBSNP += 1 
	else
		indelNOTDBSNP += 1
	end
        if infoall.include? "ESP"
                indelESP +=1 
	end
	##Indel end
   end
   lastpos = pos
end

  puts "\Sample\t#{sid}"
  print "SNV_all\t#{snvNovel+snvKnown}\n"
  print "SNV_known_any\t#{snvKnown}"
  print "\n"
  print "SNV_novel\t#{snvNovel}"
  print "\n"
  print "SNV_known_dbSNP"
  print "\t#{snvDBSNP}"
  print "\n"
  print "SNV_not_dbSNP"
  print "\t#{snvNOTDBSNP}\n"
  print "SNV_known_1KG\t#{snv1KG}\n"
  print "SNV_known_ESP\t#{snvESP}\n"
  print "SNV_known%_dbSNP\t#{(snvDBSNP/(snvNovel+snvKnown).to_f).round(2)}\n"
  print "SNV_known%_1KG\t#{(snv1KG/(snvNovel+snvKnown).to_f).round(2)}\n"
  print "SNV_known%_ESP\t#{(snvESP/(snvNovel+snvKnown).to_f).round(2)}\n"

  print "INDEL_all\t#{indelNovel+indelKnown}\n"
  print "INDEL_known_any\t#{indelKnown}\n"
  print "INDEL_novel\t#{indelNovel}\n"
  print "INDEL_known_dbSNP\t#{indelDBSNP}\n"
  print "INDEL_not_dbSNP\t#{indelNOTDBSNP}\n"
  print "INDEL_known_1KG\t#{indel1KG}\n"
  print "INDEL_known_ESP\t#{indelESP}\n"
  print "INDEL_known%_dbSNP\t#{(indelDBSNP/(indelNovel+indelKnown).to_f).round(2)}\n"
  print "INDEL_known%_1KG\t#{(indel1KG/(indelNovel+indelKnown).to_f).round(2)}\n"
  print "INDEL_known%_ESP\t#{(indelESP/(indelNovel+indelKnown).to_f).round(2)}\n"
  
  print "SNV Stats\n"
  print "known_DBSNP:ti/tv"
  print "\t#{knownti}/#{knowntv}"
  print "\n"
  print "known_DBSNP:ti/tv-ratio"
  print "\t#{(knownti/knowntv.to_f).round(3)}"
  print "\n"
  print "unknown_DBSNP:ti/tv"
  print "\t#{novelti}/#{noveltv}"
  print "\n"
  print "unknown_DBSNP:ti/tv-ratio"
  print "\t#{(novelti/noveltv.to_f).round(3)}"
  print "\n"
  print "Coding\t#{coding}\n"
  print "NonCoding\t#{snvNovel+snvKnown-coding}\n"
  print "Intronic\t#{intronic}\n"
  print "Intergenic\t#{intergenic}\n"
  print "Exonic\t#{exonic}\n"
  print "Splicing\t#{splicing}\n"
  print "Downstream\t#{downstream}\n"
  print "Upstream\t#{upstream}\n"
  print "5UTR\t#{utr5}\n"
  print "3UTR\t#{utr3}\n"
  print "ncRNA\t#{ncrna}\n"
  print "Homozygous\t#{homo}\n"
  print "Heterozygous\t#{het}\n"

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

