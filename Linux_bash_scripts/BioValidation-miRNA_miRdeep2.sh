##############################################################################
# miRNA-seq Biological Validation of Blood Serum from Field-infected Animals #
#    --- Linux bioinformatics workflow for known and novel miRNAs  ---       #
#                     -- Method 2: miRDeep2 software --                      #
##############################################################################
# Authors: Nalpas, N.C.; Correia, C.N. (2014) 
# Zenodo DOI badge: http://dx.doi.org/10.5281/zenodo.16164
# Version 1.2.0
# Last updated on: 14/07/2016

########################################################################
# Merge and uncompress miRNA-seq FASTQ files to be used with miRDeep2  #
########################################################################

# Create and enter temporary workiing directory:
mkdir $HOME/scratch/miRNAseqValidation/fastq_sequence/tmp
cd !$

# Create bash script to uncompress and merge trim.fast.gz from lanes
# 005 and 006 for each library:
for file in `find $HOME/scratch/miRNAseqValidation/fastq_sequence \
-name *_L005_R1_001_trim.fastq.gz`; \
do file2=`echo $file | perl -p -e 's/(_L005_)/_L006_/'`; \
sample=`basename $file | perl -p -e 's/(E\d+_)\w+_L\d+_R\d_\d*_trim.fastq.gz/$1/'`; \
echo "zcat $file $file2 > ${sample}trim.fastq" \
>> uncompress_merge.sh; \
done

# Run script:
chmod 755 uncompress_merge.sh
nohup ./uncompress_merge.sh &

#######################################
# Index reference genome using Bowtie #
#######################################

# Required software is Bowtie 1.1.0, consult manual for details:
# http://bowtie-bio.sourceforge.net/manual.shtml

# Create and enter the Index reference genome directory:
mkdir -p /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/bowtie1.1.0
cd !$

# Index the reference genome UMD3.1 using Bowtie:
nohup bowtie-build \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/source_file/Btau_UMD3.1_multi.fa \
Btau_UMD3.1_multi &

#############################################################
# Preprocessing of miRNA-seq data using miRDeep2: mapper.pl #
#############################################################

# Required software is miRDeep2 v.2.0.0.8, consult manual/tutorial for details:
# https://www.mdc-berlin.de/8551903/en/

# Create and enter the miRDeep2 directory for mapping work:
mkdir -p $HOME/scratch/miRNAseqValidation/mirdeep2/mapper
cd !$

# Create symbolic links to FASTQ files:
for file in \
`ls $HOME/scratch/miRNAseqValidation/fastq_sequence/tmp/*_trim.fastq`; \
do ln -s \
$file $HOME/scratch/miRNAseqValidation/mirdeep2/mapper/`basename $file`; \
done

# Run mapper.pl in one FASTQ file to see if it's working well:
mapper.pl E10_trim.fastq -e -h -m -o 3 -l 17 -r 50 -q -v -p \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/bowtie1.1.0/Btau_UMD3.1_multi \
-s test_collapsed.fa -t test.arf

# Create bash script to map miRNA reads to the reference genome:
for file in `ls *_trim.fastq`; \
do outfile=`basename $file | perl -p -e 's/_trim\.fastq//'`; \
echo "mapper.pl $file -e -h -m -o 3 -l 17 -r 50 -q -v -p \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/bowtie1.1.0/Btau_UMD3.1_multi \
-s ${outfile}_collapsed.fa -t ${outfile}.arf" \
>> mapper.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 12 mapper.sh mapper.sh.
for script in `ls mapper.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

################################################################
# Quantification of known miRNAs using miRDeep2: quantifier.pl #
################################################################

# Create and enter the working directory:
mkdir -p $HOME/scratch/miRNAseqValidation/mirdeep2/quantifier
cd !$

# Run quantifier.pl in one FASTA file to see if it's working well:
quantifier.pl -p \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/bta_hairpin-miRNA.fa \
-m /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/bta_mature-miRNA.fa \
-r $HOME/scratch/miRNAseqValidation/mirdeep2/mapper/E10_collapsed.fa -t bta

# Create a shell script to quantify the mapped mature miRNAs:
# [You will need the mature and precursor (hairpin) miRNA FASTA files for
# Bos taurus sequences only. Please refer to the 
# BioValidation-miRNA-seq_QC_filter.sh script]
for file in \
`ls $HOME/scratch/miRNAseqValidation/mirdeep2/mapper/*_collapsed.fa`; \
do outfile=`basename $file | perl -p -e 's/_collapsed.fa//'`; \
echo "mkdir -p $HOME/scratch/miRNAseqValidation/mirdeep2/quantifier/$outfile; \
cd $HOME/scratch/miRNAseqValidation/mirdeep2/quantifier/$outfile; \
quantifier.pl -p \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/bta_hairpin-miRNA.fa \
-m /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/bta_mature-miRNA.fa \
-r $file -t bta" >> quantifier.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 8 quantifier.sh quantifier.sh.
for script in `ls quantifier.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Create and enter the working directory for high confidence seqs from miRBase:
mkdir $HOME/scratch/miRNAseqValidation/mirdeep2/quantifier/high_confidence
cd !$
# Create a shell script to quantify the mapped high confidence mature and
# precursor (hairpin) miRNAs:
for file in \
`ls $HOME/scratch/miRNAseqValidation/mirdeep2/mapper/*_collapsed.fa`; \
do outfile=`basename $file | perl -p -e 's/_collapsed.fa//'`; \
echo "mkdir -p $HOME/scratch/miRNAseqValidation/mirdeep2/quantifier/high_confidence/$outfile; \
cd $HOME/scratch/miRNAseqValidation/mirdeep2/quantifier/high_confidence/$outfile; \
quantifier.pl -p \
/workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/high_conf_mature.fa \
-m /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/miRBase_fasta/high_conf_hairpin.fa \
-r $file -t bta" >> quantifier_high_conf.sh; \
done

# Split and run all scripts on Stampede:
split -d -l 8 quantifier_high_conf.sh quantifier_high_conf.sh.
for script in `ls quantifier_high_conf.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

########################################################################
# Identification of known and novel miRNAs using miRDeep2: miRDeep2.pl #
######################################################################## 

# Create and enter the miRdeep2 directory for miRdeep discovery work
mkdir -p $HOME/scratch/miRNAseqTimeCourse/mirdeep2/mirdeep
cd !$

# Copy and modify the required fasta files since miRdeep2 software requires no space in headers
cp /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/source_file/Btau_UMD3.1_multi.fa ./Btau_UMD3.1_multi.fa
cp /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/annotation_file/Btau_mature-miRNA.fa ./Btau_mature-miRNA.fa
cp /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/annotation_file/Other_mature-miRNA.fa ./Other_mature-miRNA.fa
cp /workspace/storage/genomes/bostaurus/UMD3.1_NCBI/annotation_file/Btau_pre-miRNA.fa ./Btau_pre-miRNA.fa
perl -p -i -e 's/^(>.*?)\s.*$/$1/' $HOME/scratch/miRNAseqTimeCourse/mirdeep2/mirdeep/Btau_UMD3.1_multi.fa
perl -p -i -e 's/^(>.*?) (.*?) .*$/$1_$2/' $HOME/scratch/miRNAseqTimeCourse/mirdeep2/mirdeep/Btau_mature-miRNA.fa
perl -p -i -e 's/^(>.*?) (.*?) .*$/$1_$2/' $HOME/scratch/miRNAseqTimeCourse/mirdeep2/mirdeep/Other_mature-miRNA.fa
perl -p -i -e 's/^(>.*?) (.*?) .*$/$1_$2/' $HOME/scratch/miRNAseqTimeCourse/mirdeep2/mirdeep/Btau_pre-miRNA.fa

# Create a shell script for identification and quantification of known and novel miRNAs
for file in `ls $HOME/scratch/miRNAseqTimeCourse/mirdeep2/mapper/*_collapsed.fa`; do outfile=`basename $file | perl -p -e 's/_collapsed.fa//'`; echo "mkdir -p $HOME/scratch/miRNAseqTimeCourse/mirdeep2/mirdeep/$outfile; cd $HOME/scratch/miRNAseqTimeCourse/mirdeep2/mirdeep/$outfile; miRDeep2.pl $file $HOME/scratch/miRNAseqTimeCourse/mirdeep2/mirdeep/Btau_UMD3.1_multi.fa $HOME/scratch/miRNAseqTimeCourse/mirdeep2/mapper/${outfile}.arf $HOME/scratch/miRNAseqTimeCourse/mirdeep2/mirdeep/Btau_mature-miRNA.fa $HOME/scratch/miRNAseqTimeCourse/mirdeep2/mirdeep/Other_mature-miRNA.fa $HOME/scratch/miRNAseqTimeCourse/mirdeep2/mirdeep/Btau_pre-miRNA.fa -t Cow" >> miRdeep2.sh; done;

# Split and run all scripts on Stampede
split -d -l 7 miRdeep2.sh miRdeep2.sh.
for script in `ls miRdeep2.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Collect all read counts files for transfer on desktop computer
mkdir -p $HOME/scratch/miRNAseqTimeCourse/Counts/mirdeep2/ 
cd !$
for file in `find $HOME/scratch/miRNAseqTimeCourse/mirdeep2/quantifier -name miRNAs_expressed_all_samples*`; do outfile=`echo $file | perl -p -e 's/^.*quantifier\/(.*)\/.*$/$1/'`; cp $file $HOME/scratch/miRNAseqTimeCourse/Counts/mirdeep2/${outfile}_expressed.csv; done;



# Perform subsequent miRNA analyses in R, follow R pipelines

