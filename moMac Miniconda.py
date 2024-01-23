#check your working directory
pwd

#Installing Miniconda
curl -O https://repo.anaconda.com/miniconda/Miniconda3-py311_23.5.2-0-Linux-x86_64.sh
sha256sum Miniconda3-py311_23.5.2-0-Linux-x86_64.sh
bash Miniconda3-py311_23.5.2-0-Linux-x86_64.sh -u

#Click through TandCs slowly, not to miss the “yes” promt

yes

#restart terminal

source ~/.bashrc
conda create --name my_env python=3
conda activate my_env
conda update --all

#Installing TRIMGALORE
sudo apt-get update -y
conda install -c bioconda trim-galore

#Installing KALLISTO
conda install -c bioconda kallisto
 
@Copy files into instance
gsutil cp -r gs://bucketname/ /home/user/

@Run trimming and fastqc
trim_galore /home/user/bucketname/sample1L1_1.fq.gz /home/user/bucketname/sample1L1_2.fq.gz --fastqc -output_dir /home/user/trimgalore_results/ --phred33 --paired --length 20 --quality 30 --stringency 3 --cores 8

@Download files to make Kallisto index
 
sudo apt-get install wget

wget 
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_rna_from_genomic.fna.gz

wget 
https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/latest_assembly_versions/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz 

@Make kallisto index
kallisto index -i GRCh38.p14.idx /home/user/GCF_000001405.40_GRCh38.p14_rna_from_genomic.fna.gz

#note standard k-mer length 31
 

@Run kallisto
mkdir kallisto_out

kallisto quant -i /home/user/GRCh38.p14.idx --o /home/user/kallisto_out/sample1/ --verbose --threads=8 --gtf /home/user/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz /home/user/trimgalore_results/sample1val_1.fq.gz /home/user/trimgalore_results/sample1val_2.fq.gz

cd kallisto_out

#check your files are in this directory
ls

#to download the files from the instance into your desired bucket
gsutil cp -r sample1 gs://bucketname

#Proceed with R


