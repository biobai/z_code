vcf=$1
coor=$2
outdir=$3

BASEDIR=$(dirname "$0")
[ -d $outdir ] || mkdir $outdir

if [ ! -s $outdir/query.vcf ] ; then

zcat $vcf |sed -n "/^[^#]/q;p" >$outdir/query.vcf
cut -f2 $coor |while read line; do 
  tabix $vcf $line >>$outdir/query.vcf
done

cd $outdir
#if [[ ! -s exon.maf_filtered.vcf.annovar.hg19_multianno.txt && ! -s exon.maf_filtered.vcf.annovar.final.tsv ]]; then 
  convert2annovar.pl -format vcf4old query.vcf | cut -f1-7 | awk '{gsub(",\\*", "", $0); print}'> exon.maf_filtered.vcf.avinput 
  if [ -s exon.maf_filtered.vcf.avinput ]; then
    table_annovar.pl exon.maf_filtered.vcf.avinput /data1/baiy7/tools/humandb/ -buildver hg19 -protocol refGene,avsnp147,cosmic70,exac03,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,gnomad_genome,clinvar_20180603 -operation g,f,f,f,f,f,f,f,f,f,f,f --remove --outfile exon.maf_filtered.vcf.annovar --remove
  fi
#fi

if [[ -s exon.maf_filtered.vcf.annovar.hg19_multianno.txt && ! -s exon.maf_filtered.vcf.annovar.splicing.hg19_multianno.txt ]]; then 
  python /home/shengq2/program/ngsperl/lib/Annotation/annovarSplicing.py -i exon.maf_filtered.vcf.annovar.hg19_multianno.txt -d /data1/baiy7/tools/humandb/ -o exon.maf_filtered.vcf.annovar.splicing.hg19_multianno.txt -b hg19 
fi

if [[ -s exon.maf_filtered.vcf.annovar.splicing.hg19_multianno.txt && ! -s exon.maf_filtered.vcf.annovar.final.tsv ]]; then
  sed -n '/^[^#]/q;p' query.vcf|sed '$ d' > exon.maf_filtered.vcf.annovar.final.tsv.header
  cat query.vcf | grep -v "^##" | cut -f7- > exon.maf_filtered.vcf.clean
  grep -v "^##" exon.maf_filtered.vcf.annovar.splicing.hg19_multianno.txt > exon.maf_filtered.vcf.annovar.splicing.hg19_multianno.txt.clean
  paste exon.maf_filtered.vcf.annovar.splicing.hg19_multianno.txt.clean exon.maf_filtered.vcf.clean > exon.maf_filtered.vcf.annovar.final.tsv.data
  cat exon.maf_filtered.vcf.annovar.final.tsv.header exon.maf_filtered.vcf.annovar.final.tsv.data > exon.maf_filtered.vcf.annovar.final.tsv
  rm exon.maf_filtered.vcf.clean exon.maf_filtered.vcf.annovar.splicing.hg19_multianno.txt.clean exon.maf_filtered.vcf.annovar.final.tsv.header exon.maf_filtered.vcf.annovar.final.tsv.data
fi

python /home/shengq2/program/ngsperl/lib/Annotation/filterAnnovar.py  -i exon.maf_filtered.vcf.annovar.final.tsv -t 0.1 -o exon.freq0.1

python $BASEDIR/HeteHomo/VCF_HeteHomo_reduce.py -i exon.freq0.1.filtered.missense.tsv -o exon.freq0.1.filtered.missense.tsv.HeteHomo_reduce.xls

python $BASEDIR/HeteHomo/VCF_HeteHomo_reduce.py -i exon.freq0.1.filtered.tsv -o exon.freq0.1.filtered.tsv.HeteHomo_reduce.xls
