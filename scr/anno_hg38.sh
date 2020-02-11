# /usr/bin/bash
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

  convert2annovar.pl -format vcf4old query.vcf | cut -f1-7 | awk '{gsub(",\\*", "", $0); print}'> Sporadic_WGS.maf_filtered.vcf.avinput 
  table_annovar.pl Sporadic_WGS.maf_filtered.vcf.avinput /scratch/cqs/references/annovar/humandb/ -buildver hg38 -protocol refGene,avsnp147,cosmic70,exac03,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_eur,1000g2015aug_sas,gnomad_genome,clinvar_20180603,topmed05 -operation g,f,f,f,f,f,f,f,f,f,f,f,f --remove --outfile Sporadic_WGS.maf_filtered.vcf.annovar --remove

if [[ -s Sporadic_WGS.maf_filtered.vcf.annovar.hg38_multianno.txt && ! -s Sporadic_WGS.maf_filtered.vcf.annovar.splicing.hg38_multianno.txt ]]; then 
  python /scratch/cqs_share/softwares/ngsperl/lib/Annotation/annovarSplicing.py -i Sporadic_WGS.maf_filtered.vcf.annovar.hg38_multianno.txt -d /scratch/cqs/references/annovar/humandb/ -o Sporadic_WGS.maf_filtered.vcf.annovar.splicing.hg38_multianno.txt -b hg38 
fi

if [[ -s Sporadic_WGS.maf_filtered.vcf.annovar.splicing.hg38_multianno.txt && ! -s Sporadic_WGS.maf_filtered.vcf.annovar.final.tsv ]]; then
  sed -n '/^[^#]/q;p' query.vcf |sed '$ d' > Sporadic_WGS.maf_filtered.vcf.annovar.final.tsv.header
  cat query.vcf | grep -v "^##" | cut -f7- > Sporadic_WGS.maf_filtered.vcf.clean
  grep -v "^##" Sporadic_WGS.maf_filtered.vcf.annovar.splicing.hg38_multianno.txt > Sporadic_WGS.maf_filtered.vcf.annovar.splicing.hg38_multianno.txt.clean
  paste Sporadic_WGS.maf_filtered.vcf.annovar.splicing.hg38_multianno.txt.clean Sporadic_WGS.maf_filtered.vcf.clean > Sporadic_WGS.maf_filtered.vcf.annovar.final.tsv.data
  cat Sporadic_WGS.maf_filtered.vcf.annovar.final.tsv.header Sporadic_WGS.maf_filtered.vcf.annovar.final.tsv.data > Sporadic_WGS.maf_filtered.vcf.annovar.final.tsv
  rm Sporadic_WGS.maf_filtered.vcf.clean Sporadic_WGS.maf_filtered.vcf.annovar.splicing.hg38_multianno.txt.clean Sporadic_WGS.maf_filtered.vcf.annovar.final.tsv.header Sporadic_WGS.maf_filtered.vcf.annovar.final.tsv.data
fi

python /home/shengq2/program/ngsperl/lib/Annotation/filterAnnovar.py  -i Sporadic_WGS.maf_filtered.vcf.annovar.final.tsv -t 0.1 -o Sporadic_WGS.freq0.1

python $BASEDIR/HeteHomo/VCF_HeteHomo_reduce.py -i Sporadic_WGS.freq0.1.filtered.missense.tsv -o Sporadic_WGS.freq0.1.filtered.missense.tsv.HeteHomo_reduce.xls

fi
