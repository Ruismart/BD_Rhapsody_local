######
## To run BD Rhapsody pipeline locally (2022/01, the whole process fitfully took me one month~)
##
##     test data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182001
##
##     BD ref. 
##       https://igor.sbgenomics.com/public/apps/jiewho/bd-public-project/bd-rhapsody-wta-analysis-pipeline  
##       https://scomix.bd.com/hc/en-us/articles/360047408451-BD-Rhapsody-Analysis-Pipeline-Updates
##       https://scomix.bd.com/hc/en-us/articles/360034192672-Rhapsody-WTA-Demo-Datasets
##       https://scomix.bd.com/hc/en-us/articles/360019763251-Bioinformatics-Guides
##           with CWL-runner
##             input
##               1. a CWL workflow file
##               2. a YML input specification file    
##       http://bd-rhapsody-public.s3-website-us-east-1.amazonaws.com/Rhapsody-WTA/?prefix=Rhapsody-WTA/
##           could check those files used on AWS, SampleMultiplexDemo - how do they look like
##
##     https://bitbucket.org/CRSwDev/cwl/src/master/
##         raw CWL and YML files     
##
##     env: WLSC   
##       issue1: docker is forbidden on the cluster server  
##               using singularity instead
##       issue2: NCBI GEO would automatically rewrite the fastq header info with 'SRR- accession id + order + seq length' 
##               not recognized by BD WTA but causing error running the pipeline
##               changing all fastq header to recognizable format
##       issue3: local mm10 reference files for BD WTA
##               Nmur1-fixed for truncated final exon 
## PS
##   thx to several discussions with KaiyueZhou (BD support)
##   thx to software help from LingYang, KangyongHu (WLSC)
##
######



## build local mm10, Nmur1-fixed reference index for BD WTA
# build Nmur1-fixed gtf as https://github.com/Ruismart/Nmur1-EOS/blob/main/single/Nmur1-fixed.index.sh
# code to submit
#   (STAR_2.5.2b to match formats in demo_genome used on AWS)
STAR \
--runMode genomeGenerate \
--runThreadN 10 \
--genomeDir path-to-STAR-index-output/STAR_index_vM25 \
--genomeFastaFiles path-to-mm10-refgenome/GRCm38.p6.genome.fa \
--sjdbGTFfile path-to-mm10-refgenome/gencode.VM25.chr_patch_hapl_scaff.annotation.Nmur1_fixed.gtf \
--sjdbOverhang -74 

# BD WTA requires .tar.gz as input ~
tar czvf \
STAR_index_vM25.tar.gz \
path-to-STAR-index-output/STAR_index_vM25
 
# the index tag.gz must contain one sub-directory and all index files in it, or get error~
tar -tf  STAR_index_vM25.tar.gz


## mod fastq
# re-written header from NCBI GEO, like:
#    @SRR15427821.1 1 length=60/62
# modify:
#    a. replace 'length=**' with '1/2:N:0:CGAGGCTG'  (R1/R2) (index sequence from NCBI GEO corresponding GSE->GSM->SRA)
#    b. replace 'SRR15427821.* ' with 'INST123:001:HABCDEFXX:4:1101:3325:'
#    c. 1and2
##
cd path-to-fastq/fastq_mod

zcat ../raw_link/SRR15427821_1.fastq.gz |\
sed 's/^+.*/+/g' |\
sed 's/length=../1:N:0:CGAGGCTG/g' |\
sed 's/ .* / /g' |\
sed 's/SRR15427821./INST123:001:HABCDEFXX:4:1101:3325:/g' |\
gzip > SRR15427821mod_R1.fastq.gz

zcat ../raw_link/SRR15427821_2.fastq.gz |\
sed 's/^+.*/+/g' |\
sed 's/length=../2:N:0:CGAGGCTG/g' |\
sed 's/ .* / /g' |\
sed 's/SRR15427821./INST123:001:HABCDEFXX:4:1101:3325:/g' |\
gzip > SRR15427821mod_R2.fastq.gz

# parallelly
# sample_info.txt has two columns:
#    (accession_id,index)
#     SRR15427821,CGAGGCTG
#     SRR15427822,GCTACGCT
#             ...
#

for SP in $(cat sample_info.txt)
do
SP1=$(echo ${SP} |cut -d "," -f1)
SP2=$(echo ${SP} |cut -d "," -f2)
SP3=$(echo ${SP1:0-4:4})

for NN in {1,2}
do
zcat ../raw_link/${SP1}_${NN}.fastq.gz |\
sed 's/^+.*/+/g' |\
sed "s/length=../${NN}:N:0:${SP2}/g" |\
sed 's/ .* / /g' |\
sed "s/${SP1}./INST123:001:HABCDEFXX:4:1101:${SP3}:/g" |\
gzip > ${SP1}mm_R${NN}.fastq.gz
done
done


## cwl
conda create -n cwl_py27
source activate cwl_py27

conda install -n cwl_py27 python==2.7
which pip
pip install cwlref-runner
which cwl-runner
conda deactivate


#####
## final code submit to WLSC Slurm
##    5944/5945/5946
#####

#!/bin/bash

#SBATCH -p intel-e5,amd-ep2  
#SBATCH -J BD_5944
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task 8
#SBATCH --mem-per-cpu=12GB
#SBATCH -o %j_slurm.out
#SBATCH -e %j_slurm.err

#
source /home/xuhepingLab/liushaorui/.bashrc

# load cwl, conda installed
source activate cwl_py27

# load singularity, alread installed in cluster-modules
module load singularity

# specify some path for singularity
#   that may be helpful for next runs
#   https://cwl.discourse.group/t/working-offline-with-singularity/246/18
export SINGULARITY_CACHEDIR=/storage/xuhepingLab/0.share/pipelines/BD_Rhapsody/test_script/CIDFILE
export SINGULARITY_TMPDIR=$SINGULARITY_CACHEDIR/tmp
export SINGULARITY_PULLDIR=$SINGULARITY_CACHEDIR/pull
export CWL_SINGULARITY_CACHE=$SINGULARITY_PULLDIR

# if not exist
#mkdir -p ${SINGULARITY_CACHEDIR}
#mkdir -p ${SINGULARITY_TMPDIR}
#mkdir -p ${SINGULARITY_PULLDIR}
#mkdir -p ${CWL_SINGULARITY_CACHE}

# specify some url to make the internet available
#     (special/temporary setting for WLSC clusters)
export http_proxy=http://192.168.103.207:3128
export https_proxy=http://192.168.103.207:3128

#
cwl-runner \
--outdir ./GSM5515944 \
--tmpdir-prefix ./tmp \
--cachedir ./tmp \
--cidfile-dir $SINGULARITY_CACHEDIR \
--singularity \
rhapsody_WTA_1.9.1.cwl template_GSM5515944_1.9.1.yml

## 

####
## this script took over 5 days to run ~
##    I didn't check all the details (e.g. possible inefficient resource application in Slurm-singularity), 
##    but if the env was fine, it might have taken x10 time more than Cellranger to run a similar datasize,
##      (I did run a few small demos before that ~)
##    there's large space for its performance improvement ...
##    however, it's done.
#### 


## after BD_WTA pipeline
# check matrix files: output/GSM5515944/***.st
#     manually compare a few cells-genes with online processed matrix
#         matched
#     output combined bam
#         IGV check
          
cd output/GSM5515944
samtools view -h Combined_SRR15427821mm_final.BAM chr1:86383500-86428228 |\
samtools view -bS |samtools sort - -o Nmur1.bam
samtools index Nmur1.bam

####
## processing is done
####















