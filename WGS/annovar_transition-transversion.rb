#!/usr/bin/env ruby

## calcualte transition/transversion ratio from VCF file
## per sample
## also count the number of synonymous, missense, silent var per sample

def main
  infile = ARGV[0]
  somatic = ARGV[1]
#Class @collumn0
#"downstream"
#"exonic"
#"exonic;splicing"
#Func
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
##Categories in the infile at column2
#ExonicFunc
#frameshift deletion
#frameshift insertion
#nonframeshift deletion
#nonframeshift insertion
#nonsynonymous SNV	= missense
#stopgain SNV		= nonsense
#stoploss SNV		= readthrough
#synonymous SNV	 	= silent
#unknown


  novelti = {}
  noveltv = {}
  knownti = {}
  knowntv = {}
  synon = {}
  missense = {}
  nonsense = {}
  readthrough = {}
  indel = {}
  homo = {}
  het  = {}
  onekdb = {}
  sid = [] # sample array

  lastpos = -1 

#     samples.each do |sname|
if( somatic.nil? ) 
	bname= File.basename(infile)
	infilearr =  bname.chomp.split(/\./)
#        cc = sname.sub(" ","_")
	cc = infilearr[0].sub(" ","_")        
        sid << cc
        novelti[cc] = 0
        noveltv[cc] = 0
        knownti[cc] = 0 
        knowntv[cc] = 0
        synon[cc] = 0
        missense[cc] = 0
        nonsense[cc] = 0
        readthrough[cc] = 0
	indel[cc] = 0
        homo[cc] = 0
        het[cc] = 0
	onekdb[cc] = 0
else 
	cc = "somatic"
        sid << cc
        novelti[cc] = 0
        noveltv[cc] = 0
        knownti[cc] = 0
        knowntv[cc] = 0
        synon[cc] = 0
        missense[cc] = 0
        nonsense[cc] = 0
        readthrough[cc] = 0
        indel[cc] = 0
        homo[cc] = 0
        het[cc] = 0
        onekdb[cc] = 0
	cc = "germline"
        sid << cc
        novelti[cc] = 0
        noveltv[cc] = 0
        knownti[cc] = 0
        knowntv[cc] = 0
        synon[cc] = 0
        missense[cc] = 0
        nonsense[cc] = 0
        readthrough[cc] = 0
        indel[cc] = 0
        homo[cc] = 0
        het[cc] = 0
        onekdb[cc] = 0
	fs = File.new(infile+"_somatic", 'w+')
	fg = File.new(infile+"_germline",'w+')
end

  
#0       Func    exonic
#1       Gene    ZMYND11
#2       ExonicFunc      synonymous SNV
#3       AAChange        NM_001202465:c.A1356G:p.E452E
#4       Conserved       709;Name=lod=1016
#5       SegDup
#6       ESP5400_ALL     0.94841
#7       1000g2012feb_ALL        0.94
#8       dbSNP132        rs1017361
#9       AVSIFT
#10      LJB_PhyloP
#11      LJB_PhyloP_Pred
#12      LJB_SIFT
#13      LJB_SIFT_Pred
#14      LJB_PolyPhen2
#15      LJB_PolyPhen2_Pred
#16      LJB_LRT
#17      LJB_LRT_Pred
#18      LRT_MutationTaster
#19      LRT_MutationTaster_Pred
#20      LJB_GERP++
#21      Chr     10
#22      Start   294953
#23      End     294953
#24      Ref     A
#25      Obs     G
#26      Otherinfo       10
#27              294953
#28              .
#29              A
#30              G
#31              222
#32              .
#33              "DP=30;AF1=1;AC1=2;DP4=0,0,17,10;MQ=39;FQ=-108"
#34              GT:PL:GQ
#35              "1/1:255,81,0:99"
#
File.new(infile, 'r').each do |line|
    cols=line.chomp.split(/\t/)
    next if cols[0].include? "Func"
    bclass, fclass,pos,name,ref,alt,onek =  cols[0].gsub(/"/,''),cols[2].gsub(/"/,''),cols[22].to_i,cols[8],cols[24].gsub(/"/,''),cols[25].gsub(/"/,''), cols[7]
if( somatic.nil? ) 
      info,gt = cols[-3].gsub(/"/,''),  cols[-1..-1]
#    elsif samples.size == 2
    else	
      info,gt = cols[-4].gsub(/"/,''),  cols[-2..-1]
#    else 
#      print "Cannot have >2 samples"
    end

      ## SYN AND NONSYN
      synonFlag, missenseFlag, nonsenseFlag, knowntiFlag, knowntvFlag ,noveltiFlag, noveltvFlag, readthroughFlag,indelFlag,noncodingFlag,onekFlag =0,0, 0,0,0, 0, 0, 0,0, 0,0
      coding = 0

    if bclass == "exonic"
	if  fclass == "nonsynonymous SNV"	# missense
          missenseFlag = 1
          coding = 1
	elsif fclass == "synonymous SNV" 	#silent
          synonFlag = 1
          coding = 1
        elsif fclass == "stoploss SNV"        #readthrough
          readthroughFlag = 1
          coding = 1
        elsif fclass == "stopgain SNV"        #nonsense
          nonsenseFlag = 1
          coding = 1
	elsif fclass.include? "frameshift" 
          coding = 1 
	  indelFlag  =1
	end
    end
     
      next if pos == lastpos       #ignore same variant i.e. same position in tht chromose
      next if coding == 0
        
   if indelFlag == 0
      ti = judge(ref,alt)
      if name =~ /^rs/  # known
        if ti == 1
          knowntiFlag = 1
        else
          knowntvFlag = 1
        end
      else  # novel
        if ti == 1
          noveltiFlag = 1
        else
          noveltvFlag = 1
        end
      end
      if onek != ''
	onekFlag = 1
	end
   end

if( somatic.nil? ) 
#      gt.each do |genotypeinfo|
#      genotype = genotypeinfo.gsub(/"/,'').split(":")[0]
      genotype = gt[0].gsub(/"/,'').split(":")[0]
      subject = sid[0]
        if genotype == '0/1'  ## het
          het[subject] += 1
        elsif genotype == '1/1' ## homo
          homo[subject] += 1
        end
        
        if genotype == '0/1' or genotype == '1/1'
          novelti[subject] += noveltiFlag
          noveltv[subject] += noveltvFlag
          knownti[subject] += knowntiFlag
          knowntv[subject] += knowntvFlag
          synon[subject] += synonFlag
          missense[subject] += missenseFlag
          nonsense[subject] += nonsenseFlag
          readthrough[subject] += readthroughFlag
          indel[subject] += indelFlag
          onekdb[subject] += onekFlag
        end
  else
 	g1 = gt[0].gsub(/"/,'').split(":")[0]
        g2 = gt[1].gsub(/"/,'').split(":")[0]
        if g1 != g2
	  subject = sid[0]
	  fs.puts "#{line}"
	else
	  subject = sid[1]
          fg.puts "#{line}"
	end
          novelti[subject] += noveltiFlag
          noveltv[subject] += noveltvFlag
          knownti[subject] += knowntiFlag
          knowntv[subject] += knowntvFlag
          synon[subject] += synonFlag
          missense[subject] += missenseFlag
          nonsense[subject] += nonsenseFlag
          readthrough[subject] += readthroughFlag
          indel[subject] += indelFlag
	  onekdb[subject] += onekFlag
  end
  lastpos = pos
end

  puts "\Samples\t#{sid.join("\t")}"
  print "all_SNV"
  sid.each {|s| print "\t#{knownti[s] + knowntv[s] + novelti[s] + noveltv[s]}"}
  print "\n"
  
  print "all_known_SNV(dbSNP132)"
  sid.each {|s| print "\t#{knownti[s] + knowntv[s]}"}
  print "\n"

  print "all_novel_SNV"
  sid.each {|s| print "\t#{novelti[s] + noveltv[s]}"}
  print "\n"

  print "all_known%_SNV(dbSNP132)"
  sid.each {|s| print "\t#{((knownti[s] + knowntv[s])/(knownti[s] + knowntv[s] + novelti[s] + noveltv[s]).to_f).round(2)}"}
  print "\n"

  print "all_1KG_SNV"
  sid.each {|s| print "\t#{onekdb[s]}"}
  print "\n"

  print "all_1KG%_SNV"
  sid.each {|s| print "\t#{(onekdb[s]/(knownti[s] + knowntv[s] + novelti[s] + noveltv[s]).to_f).round(2)}"}
  print "\n"

  print "SNV known:ti/tv"
  sid.each {|s| print "\t#{knownti[s]}/#{knowntv[s]}"}
  print "\n"

  print "SNV known:ti/tv-ratio"
  sid.each {|s| print "\t#{(knownti[s]/knowntv[s].to_f).round(3)}"}
  print "\n"


  print "SNV novel:ti/tv"
  sid.each {|s| print "\t#{novelti[s]}/#{noveltv[s]}"}
  print "\n"

  print "SNV novel:ti/tv-ratio"
  sid.each {|s| print "\t#{(novelti[s]/noveltv[s].to_f).round(3)}"}
  print "\n"

  print "silent"
  sid.each {|s| print "\t#{synon[s]}" }
  print "\n"
  
  print "missense"
  sid.each {|s| print "\t#{missense[s]}" }
  print "\n"

  print "nonsense"
  sid.each {|s| print "\t#{nonsense[s]}" }
  print "\n"

  print "readthrough"
  sid.each {|s| print "\t#{readthrough[s]}" }
  print "\n"

  print "indel"
  sid.each {|s| print "\t#{indel[s]}" }
  print "\n"
  print "homozygous"
  sid.each {|s| print "\t#{homo[s]}" }
  print "\n"

  print "heterozygous"
  sid.each {|s| print "\t#{het[s]}" }
  print "\n"

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

