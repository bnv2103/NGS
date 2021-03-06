#!/usr/bin/env ruby
require 'getoptlong'

### !!! need to take care of missing values

## convert VCF to a VT Test input files

def main
  optHash = getopt()
  vcf, pheno = optHash["--vcf"], optHash["--phenotype"]

  switch = 1
  if optHash.key?("--inclusive") ## 
    switch = 0
  end

  
  samples = readPheno(pheno)
  gt = {}
#  snp = {}
  readVCF(gt, samples,  vcf,switch)  # gt: gene -> pos -> sample -> genotype, 
  $stderr.puts "# of genes: #{gt.keys.size}"
  printout(gt, vcf)
  
#  doublehits(gt,samples)
  
end




def printout(gt, prefix)
  snpo = File.new(prefix + ".csnp", "w")
  gto = File.new(prefix + ".cgeno", "w")
  gt.keys.sort.each do |gene|
    num = 0
    gt[gene].keys.sort.each do |pos|
      num += 1
      snv = gt[gene][pos]
      # snpo.puts "#{gene}\t#{num}\t0.5\t#{snv[:chr]}\t#{pos}\t#{snv[:fclass]}\t#{snv[:id]}"  # polyphen score is set to 0.5
      majorAllele = 0
      alleles = {0=> 0, 1=>0, 2=> 0, -9 => 0}
      gt[gene][pos][:controls].values.each do |a|
        if a >= 0
          alleles[a] += 1
        end
      end
      
      if alleles[2] > alleles[0]
        majorAllele = 2
      end
      
      f0=0
      f1=0
      g1 = [0,0,0,0]
      g0 = [0,0,0,0]
      snv[:cases].keys.sort.each do |sample|
        if  snv[:cases][sample] !=  majorAllele # carrier of a minor allele
          allele = (snv[:cases][sample] - majorAllele).to_i.abs
          gto.puts "#{gene}\t#{sample}\t#{num}\t#{allele}"
          f1 += 1
          ## if snv[:cases][sample] 
        end
      end
      snv[:controls].keys.sort.each do |sample|
        if  snv[:controls][sample] !=  majorAllele # carrier of a minor allele
          allele = (snv[:controls][sample] - majorAllele).to_i.abs
          gto.puts "#{gene}\t#{sample}\t#{num}\t#{allele}"
          f0 += 1
        end

      end

      if f0 >0 or f1 > 0 
        snpo.puts "#{gene}\t#{num}\t0.5\t#{snv[:chr]}\t#{pos}\t#{snv[:fclass]}\t#{snv[:id]}\t#{f1}\t#{f0}"  
      end

    end
  end
  snpo.close
  gto.close
end

def readVCF(gt,samples,  vcf, switch)
  sid = []
  flag = 0  # by default it is vcf
  $stderr.puts vcf
  File.new(vcf, 'r').each do |line|
    next if line.match("^##")
    
    cols = line.chomp.split(/\s+/) 

    if cols.size < 2  # a list of VCF
      flag = 1
    else
      flag = 0
    end

    if flag == 1 # 
#      $stderr.puts line
#      $stderr.puts cols.size
      readVCF(gt, samples, cols[0], switch)
    else
      if line.match("^#CHROM")  # header of VCF
        sid=cols[9..-1]
      else
        chr,pos,id,qual,filter,info, gtdetails = cols[0], cols[1].to_i, cols[2], cols[5].to_f, cols[6], cols[7], cols[9..-1]
        fclass = 0  # 0: synonymous; >0: non-syn
        geneName = ""

        info.split(';').each do  |item|
          k,v = item.split('=')[0..1]
          if k =~ /refseq.changesAA/
            if v == "true"  # non-syn
              fclass = 1
            end
          elsif k =~ /refseq.name2/
            geneName = v
            if !gt.key?(geneName) 
              gt[geneName] = {}

            end
          elsif k=~ /refseq.functionalClass/
            if v == "nonsense" or v == "readthrough" # set it as damaging by default
              fclass = 2
            end
          end
        end
           
  
        
        if geneName != "" and ( fclass > 0 or switch == 0 )
          if !gt[geneName].key?(pos)
            snv = {}
            snv[:chr] = chr
            if id.match("^rs")
              snv[:id] = id 
            else
              snv[:id] = "#{chr}_#{pos}"
            end
            #snv[:pos] = pos
            snv[:fclass] = fclass
            snv[:controls] = {}
            snv[:cases] = {}
            # snv[:gene] = geneName
            gt[geneName][pos]= snv
          else
            snv = gt[geneName][pos]
          end

          0.upto(gtdetails.size-1) do |i|
            gtc = gtdetails[i].split(':')
            sname = sid[i]
#            $stderr.puts sname
#            $stderr.puts line
            next unless samples.key?(sname)
            a = encodeGT(gtc[0])
            # next if a < 0 # missing value
#            gt[geneName][pos][sname] = a
            if samples[sname] < 0 # control
              snv[:controls][sname] = a
            else
              snv[:cases][sname] = a
            end
          end
        end
      end
    end
  end
  return gt
end

def encodeGT(vcfcode)
  if vcfcode == '0/0'
    return 0
  elsif vcfcode == "0/1"
    return 1
  elsif vcfcode == "1/1"
    return 2
  else
    return -9
  end
end

def readPheno(pheno)
  samples = {}
  File.new(pheno, 'r').each do |line|
    next if line.match("^#") # header line

    cols = line.chomp.split(/\s+/)
    sid, status = cols[0], cols[-1].to_i
    if status >= 1 # cases
      samples[sid] = 1
    else
      samples[sid] = -1
    end
  end
  return samples 
end


def getopt

  opts = GetoptLong.new(
                        ["--vcf", "-v", GetoptLong::REQUIRED_ARGUMENT],
                        ["--phenotype", "-p", GetoptLong::REQUIRED_ARGUMENT],
                        ["--inclusive", "-i", GetoptLong::NO_ARGUMENT],
                        ["--help", "-h", GetoptLong::NO_ARGUMENT]
                        )
  optHash = {}
  opts.each do |opt, arg|
    optHash[opt] = arg
  end
  if optHash.key?("--help") or !optHash.key?("--vcf") or !optHash.key?("--phenotype")
    $stderr.puts "Usage: ruby __.rb -v VCF_or_list -p phenotype [-i]"
    $stderr.puts " if -n is provided, only non-synonymous variants will be considered"
    exit
  end
  return optHash
  
end

main()
