#!/usr/bin/env ruby

## calcualte transition/transversion ratio from genome_summary VCF file for 1 sample

def main
 infile = ARGV[0]

lastpos = -1 

bname= File.basename(infile)
homo, het,  intergenic, intronic, downstream,upstream,splicing, exonic, utr3, utr5 = 0,0,0,0,0,0,0,0,0,0,0
ncrna_exonic, ncrna_intronic, ncrna_splicing, ncrna_utr3,ncrna_utr5 = 0,0,0,0,0
indelKnown, indelNovel, indel1KG, indelDBSNP, indelNOTDBSNP ,indelESP =0,0,0,0,0,0
anyti,anytv,nonetv,noneti,snvKnown,snvNovel, knownti, knowntv ,novelti, noveltv, snvDBSNP,snvNOTDBSNP,snv1KG,snvESP  = 0,0,0,0,0,0,0,0,0,0,0,0,0,0
codingUnknown, codingFrameDel,   codingFrameIns ,codingNonFrameDel, codingNonFrameIns,  codingNonFrameSub, codingMissense, codingNonsense, codingReadthrough, codingSynonymous = 0,0,0,0,0,0,0,0,0,0

sid=""

File.new(infile, 'r').each do |line|
    cols=line.chomp.split(/\t/)
    next if ( cols[0] =~ /##/) == 0
    if ( cols[0].downcase =~ /#chr/ ) == 0
	sid = cols[9]
	next
    end
		
	pos,name,ref,alt,infoall,genotype = cols[1],cols[2],cols[3],cols[4],cols[7],cols[9].split(":")[0]
	info = infoall.split(/;/)
	
	next if pos == lastpos       #ignore same variant i.e. same position in tht chromose

   if ref.size == alt.size and (not ref.include? "-")  and (not alt.include? "-")  and  ref.size == 1 and alt.size == 1  # it is SNV
      	ti = judge(ref,alt)
	if not  infoall.include? "dbSNP" and not  infoall.include? "1KG" and not infoall.include? "ESP"
		snvNovel += 1
                if ti == 1
                  noneti += 1
                else
                  nonetv += 1
                end
	else
		snvKnown += 1
                if ti == 1
                  anyti += 1
                else
                  anytv += 1
                end
	end
      	if infoall.include? "dbSNP"  # knownDBSNP ti/tv
        	if ti == 1
	          knownti += 1
        	else
	          knowntv += 1
        	end
		snvDBSNP += 1
	else  # novel ti/tv
        	if ti == 1
	          novelti += 1
        	else
	          noveltv += 1
        	end
		snvNOTDBSNP += 1
	end
	snv1KG += 1 if infoall.include? "1KG"
      ## SYN AND NONSYN
	intergenic += 1 if infoall.include? "function=intergenic"
	intronic += 1 if infoall.include? "function=intronic"
	downstream += 1 if infoall.include? "function=downstream"
	upstream += 1 if infoall.include? "function=upstream"
	if infoall.include? "function=exonic"
		exonic += 1
		snvESP += 1 if infoall.include? "ESP"
		codingUnknown += 1 if infoall.include? "functionalClass=unknown"
		codingFrameDel += 1 if infoall.include? "functionalClass=frameshift deletion"
		codingFrameIns += 1 if infoall.include? "functionalClass=frameshift insertion"
		codingNonFrameDel += 1 if infoall.include? "functionalClass=nonframeshift deletion"
		codingNonFrameIns += 1 if infoall.include? "functionalClass=nonframeshift insertion"
		codingNonFrameSub += 1 if infoall.include? "functionalClass=nonframeshift substitution"
		codingMissense += 1 if infoall.include? "functionalClass=nonsynonymous SNV"
		codingNonsense += 1 if infoall.include? "functionalClass=stopgain SNV"
		codingReadthrough += 1 if infoall.include? "functionalClass=stoploss SNV"
		codingSynonymous += 1 if infoall.include? "functionalClass=synonymous SNV"
	end
	
	utr3 += 1 if infoall.include? "function=UTR3"
	utr5 += 1 if infoall.include? "function=UTR5"
	if infoall.include? "=splicing" and (not infoall.include? "=ncRNA_splicing")
		        splicing += 1
	end
	ncrna_exonic += 1 if infoall.include? "function=ncRNA_exonic"
        ncrna_intronic += 1 if infoall.include? "function=ncRNA_intronic"
        ncrna_splicing += 1 if infoall.include? "function=ncRNA_splicing"
        ncrna_utr3 += 1 if infoall.include? "function=ncRNA_UTR3"
        ncrna_utr5 += 1 if infoall.include? "function=ncRNA_UTR5"

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
	indel1KG += 1 if infoall.include? "1KG"
        if infoall.include? "dbSNP"  
		indelDBSNP += 1 
	else
		indelNOTDBSNP += 1
	end
	indelESP +=1 if infoall.include? "ESP"
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
  print "SNV_known%_ESP\t#{(snvESP/(exonic).to_f).round(2)}\n"

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
  print "known_DBSNP:ti/tv\t#{knownti}/#{knowntv}\n"
  print "known_DBSNP:ti/tv-ratio\t#{(knownti/knowntv.to_f).round(3)}\n"
  print "unknown_DBSNP:ti/tv\t#{novelti}/#{noveltv}\n"
  print "unknown_DBSNP:ti/tv-ratio\t#{(novelti/noveltv.to_f).round(3)}\n"
  print "known_ANY:ti/tv\t#{anyti}/#{anytv}\n"
  print "known_ANY:ti/tv-ratio\t#{(anyti/anytv.to_f).round(3)}\n"
  print "NOVEL:ti/tv\t#{noneti}/#{nonetv}\n"
  print "NOVEL:ti/tv-ratio\t#{(noneti/nonetv.to_f).round(3)}\n"
  print "\n#AnnotatedCodingSNV\t#{exonic}\n"
  print "Missense\t#{codingMissense}\n"
  print "Nonsense\t#{codingNonsense}\n"  
  print "ReadThrough\t#{codingReadthrough}\n"    
  print "Synonymous\t#{codingSynonymous}\n"      
  print "Frameshift_Deletion\t#{codingFrameDel}\n"       
  print "Frameshift_Insertion\t#{codingFrameIns}\n"       
  print "NonFrameshift_Deletion\t#{codingNonFrameDel}\n"   
  print "NonFrameshift_Insertion\t#{codingNonFrameIns}\n"   
  print "NonFrameshift_Substitution\t#{codingNonFrameSub}\n" 
  print "Function_Unknown\t#{codingUnknown}\n"
  print "\n#UnannotatedSNV\t#{ncrna_exonic + ncrna_intronic + ncrna_splicing+ ncrna_utr3+ncrna_utr5}\n"
  print "Unannotated_Coding\t#{ncrna_exonic}\n"
  print "Unannotated_Splicing\t#{ncrna_splicing}\n"
  print "Unannotated_Intronic\t#{ncrna_intronic}\n"
  print "Unannotated_5UTR\t#{ncrna_utr5}\n"
  print "Unannotated_3UTR\t#{ncrna_utr3}\n"

  print "\n#NonCodingSNV\t#{snvNovel+snvKnown-exonic}\n"
  print "Intronic\t#{intronic}\n"
  print "Intergenic\t#{intergenic}\n"
  print "Splicing\t#{splicing}\n"
  print "Downstream\t#{downstream}\n"
  print "Upstream\t#{upstream}\n"
  print "5UTR\t#{utr5}\n"
  print "3UTR\t#{utr3}\n"
  print "HomozygousSNV\t#{homo}\n"
  print "HeterozygousSNV\t#{het}\n"

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

