input_vcf=$1
input_bam=$2
input_barcode=$3

# inputvcf: xxxxxx.T.merge.filter3.snp_indel.vcf

vcfdir=`dirname $input_vcf`
xls_annotated=`basename $input_vcf .vcf`.xls.with.simple.cluster
output_name=$4

awk '{FS="\t";OFS="\t";if($(NF-1)>0.02){OFS="\t";print $63, $64, $65, $0}}' $vcfdir/$xls_annotated |awk '{if($0~"Clone"||$0~"exonic"){print $0}}'   |sort -k1,1n > $output_name.candidate.snp.tsv
cut -f 13-15 $output_name.candidate.snp.tsv |sort|uniq  > $output_name.candidate.snp.bed
bedtools intersect -a $input_vcf -b $output_name.candidate.snp.bed | sed -e 's/^/chr/g' |awk 'BEGIN{start=0}{if($2!=start){print;start=$2}}' > $output_name.candidate.snp.vcf.body
grep "#" $input_vcf | sed -e 's/contig=<ID=/contig=<ID=chr/g' > $output_name.candidate.snp.vcf.head

grep -v "#" $vcfdir/seg.tab.should.remain.germline.vcf  | sed -e 's/^/chr/g' > $output_name.candidate.germline.snp.vcf.body

cat $output_name.candidate.snp.vcf.body $output_name.candidate.germline.snp.vcf.body | sort -k1,1V -k2,2n |cut -f1-8 > $output_name.candidate.firstcols.vcf.body 

cat $output_name.candidate.snp.vcf.head  $output_name.candidate.firstcols.vcf.body |cut -f1-8 |uniq  > $output_name.candidate.snp.vcf
/gpfs/bin/anaconda3-loompy/bin/CrossMap.py vcf /gpfs/genomedb/chains/hg19ToHg38.over.chain $output_name.candidate.snp.vcf /gpfs/genomedb/cellranger/refdata-cellranger-atac-GRCh38-1.2.0/fasta/genome.fa $output_name.candidate.snp.hg38.vcf


# /gpfs/output/Bladder_Cancer_Landscape_Project/CT/dir//2004920S1OWC-2004920S3OWC/seg.tab.should.remain.germline.vcf




cat $output_name.candidate.snp.hg38.vcf |egrep -v "liftOverProgram=|liftOverChainFile=|originalFile=|targetRefGenome=|liftOverDate=" > $output_name.candidate.snp.hg38.vcf1
mv $output_name.candidate.snp.hg38.vcf1 $output_name.candidate.snp.hg38.vcf

/gpfs/bin/vartrix/vartrix_linux   --cell-barcodes $input_barcode --fasta /gpfs/genomedb/cellranger/refdata-cellranger-atac-GRCh38-1.2.0/fasta/genome.fa --bam $input_bam --vcf $output_name.candidate.snp.hg38.vcf --threads `nproc` --scoring-method coverage --out-matrix $output_name.alt.mtx --ref-matrix $output_name.ref.mtx


