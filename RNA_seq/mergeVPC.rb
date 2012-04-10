#!/usr/bin/ruby

def main()
  # dir = "/ifs/scratch/c2b2/ngs_lab/xs2182/rubyCode/2T/"
  infile = File.new("var_flt_vcf.exonic_variant_function", "r")
  # outfile = File.new("junk.jnk", "w")
  snv = Hash.new
  infile.each {
    |line|
    cols = line.chomp.split(/\t/)
    snv[cols[4]] = {"type"=>cols[1], "exon"=>cols[2], "chrom"=>cols[3], "end"=>cols[5], "ref"=>cols[6], "alt"=>cols[7], "genotype"=>cols[8], "Qual"=>cols[9], "DP"=>cols[10], "MQ"=>cols[11]}
  }
  
  infile.close()
  # keys = snv.keys
  # puts snv.length
  # infile.close
  var = Hash.new
  infile1 = File.new("var_flt.vcf", "r")
  infile1.each {
    |line|
    if line.match("#")
      #puts line
    else
      # puts line
      cols=line.chomp.split(/\s+/)
      #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  111222_lane_4_1N_1
      #snv[cols[4]] = {"type"=>cols[1], "exon"=>cols[2], "chrom"=>cols[3], "end"=>cols[5], "ref"=>cols[6], "alt"=>cols[7]
      
        qual, info, pl = cols[5].to_f, cols[7].split(';'), cols[9].split(':')
        var[cols[1]] = {"chrom"=>cols[0], "DP4"=>info[3], "GT"=>pl[0], "PL"=>pl[1], "GQ"=>pl[2]}
    end
  }

  infile1.close()
  puts var.length

  outfile = File.new("var_summary.txt", "w")
  outfile.write "Type\tExon\tChrom\tStat\tEnd\tRef\tAlt\tGenotype\tQUAL\tDP\tDP4\tGT\tPL\tGQ\n"
  keys = snv.keys
  for key in 0...keys.length
    if snv.has_key?(keys[key])
      temp = keys[key]
      if var.has_key?(temp)
        # temp = keys[key]
        print "key  : ", keys[key], "\n"
        outfile.write  snv[keys[key]]["type"] + "\t"
        outfile.write snv[keys[key]]["exon"]+ "\t"
        outfile.write snv[keys[key]]["chrom"]+ "\t"
        outfile.write keys[key]+ "\t"
        outfile.write snv[keys[key]]["end"]+ "\t"
        outfile.write snv[keys[key]]["ref"]+ "\t"
        outfile.write snv[keys[key]]["alt"]+ "\t"
        outfile.write snv[keys[key]]["genotype"]+ "\t"
        outfile.write snv[keys[key]]["Qual"]+ "\t"
        outfile.write snv[keys[key]]["DP"]+ "\t"
        #outfile.write snv[keys[key]]["MQ"]+ "\t"
        outfile.write var[temp]["DP4"]+ "\t"        
        outfile.write var[temp]["GT"]+ "\t"
        outfile.write var[temp]["PL"]+ "\t"
        outfile.write var[temp]["GQ"]+ "\n"
      end
    end
   # print "\n\n"
  end
  outfile.close

end


main()
