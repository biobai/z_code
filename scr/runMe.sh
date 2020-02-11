##
BASEDIR=$(dirname "$0")

query=$1
outdir=$2

bash $BASEDIR/coordinate/query/q2c.sh -i $query -g $BASEDIR/coordinate/geneList/Homo_sapiens.GRCh37.87.geneList -o $outdir/genes_coordinates_hg19.txt
bash $BASEDIR/coordinate/query/q2c.sh -i $query -g $BASEDIR/coordinate/geneList/Homo_sapiens.GRCh38.87.geneList -o $outdir/genes_coordinates_hg38.txt
sed -i 's/\t/\tchr/' $outdir/genes_coordinates_hg38.txt
bash $BASEDIR/anno_hg19.sh /scratch/cqs/baiy7/Tim_proj/WES/analysis/result_exome/bwa_refine_gatk4_hc_gvcf_vqsr/result/exon.indels.snp.recal.pass.norm.nospan.vcf.gz $outdir/genes_coordinates_hg19.txt $outdir/hg19
bash $BASEDIR/anno_hg38.sh /scratch/cqs/baiy7/Tim_proj/Sporadic_WGS/analysis/Sporadic_WGS_pipeline_result/bwa_refine_gatk4_hc_gvcf_vqsr/result/Sporadic_WGS.indels.snp.recal.pass.norm.nospan.vcf.gz $outdir/genes_coordinates_hg38.txt $outdir/hg38

