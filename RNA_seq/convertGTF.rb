#!/usr/bin/ruby

def main()
  infile = File.new("genes1.gtf", "r")
  outfile = File.new("genes.gtf", "w")
  count = 1
  infile.each {
    |line|
    cols = line.chomp.split(/\t/)
    
    array_length= cols.length
    for key in 0...array_length-1
      outfile.write cols[key]+ "\t"
    end
    
    info = cols[array_length-1].split(';')
    names = info[0].split('"')
    # gene_id "CG11023"; gene_name "CG11023"; p_id "P17192"; transcript_id "NM_175941"; tss_id "TSS854";
    
    transID = " "
    if cols[2]=="exon"
      transID_term = info[1].split(' ')
      # outfile.write transID_term[1]
      transID = transID_term[1].chomp
      # outfile.write "transcript_id " + "\"" + transID[0] + "\"" + ";"
      outfile.write "gene_id \"" + names[1] + "\";" + " gene_name \"" + names[1] + "\";" + " transcript_id  " + "\"" + transID + "\"" + ";" +  " p_id " + "\"" + transID + "\"" + ";" + " tss_id \"" + count.to_s() + "\";" 
    end
    ptID = " "
    if cols[2]=="CDS"
      transID_term = info[1].split(' ')
      # outfile.write transID_term[1]
      ptID = transID_term[1].chomp
      # outfile.write "transcript_id " + "\"" + transID[0] + "\"" + ";"
      outfile.write "gene_id \"" + names[1] + "\";" + " gene_name \"" + names[1] + "\";" + " transcript_id " + "\"" + ptID + "\"" + ";" + " p_id " + "\"" + ptID + "\"" + ";" + " tss_id \"" + count.to_s() + "\";"
    end
    
    if cols[2]=="start_codon"
      outfile.write "gene_id \"" + names[1] + "\";" + " gene_name \"" + names[1] + "\";"+ " tss_id \"" + count.to_s() + "\";"
      # outfile.write "gene_id \"" + names[1] + "\";" + " gene_name \"" + names[1] + "\";" + " p_id " + "\"" + ptID + "\"" + ";" + " transcript_id  " + "\"" + transID + "\"" + ";" + " tss_id \"" + " " + "\";"
    end 
    count = count + 1
    outfile.write "\n"
  }

  infile.close()
  outfile.close()
  # array_length= cols.length
  # puts array_length

end


main()
