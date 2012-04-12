BEGIN { OFS="\t" }
# filter snp qual > 20/30 if not a homozygote
# cap the coverage at 1000
$8 <= 1000 {
  if (tolower($3) == tolower($4)) # homozygote
    print $1, $2, $3, $4, $5, $6, $7, $8;
  else if ($6 >= 20) # heterozygote with high SNP quality
    print $1, $2, $3, $4, $5, $6, $7, $8;
}
$8 > 1000 {
  if (tolower($3) == tolower($4)) # homozygote
    print $1, $2, $3, $4, $5, $6, $7, 1000;
  else if ($6 >= 20) # heterozygote with high SNP quality
    print $1, $2, $3, $4, $5, $6, $7, 1000;
}