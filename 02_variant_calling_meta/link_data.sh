. ./vars.sh

# Reference sequence
cp /n/home08/arturrc/genomes/GCF_000146045.2_R64_genomic.fna data/

# Import fastqgz
for i in {0..69}
do
	label=${meta_fastqgz_labels[$i]}
	fastqgz_dir=${meta_fastqgz_dirs[$i]}
	ln -nsf ${fastqgz_dir}/${label}*_R*_*.fastq.gz data/
done
