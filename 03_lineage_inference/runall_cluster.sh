# 0) Initiate log
touch batch.log

# 1) Import variables
# Defines scratch and outdir paths that are used by all scripts.
# scratch is for temporary files.
# outdir is for the final output of inferred frequencies.
. ./vars.sh

# 2) Set up tmp folder
mkdir -p tmp

# 3) Submit jobs for recursive fits
# UCP 2018
expt=UCP_2018
for tp_id in {1..21}
do
	echo \> ${expt} tp_id=${tp_id}
	sbatch --export=expt=${expt},tp_id=${tp_id} batch/estimate_lineage_freqs_recursively.sbatch
	sleep 0.2
done

# UCP 2019
expt=UCP_2019
for tp_id in {1..15}
do
	echo \> ${expt} tp_id=${tp_id}
	sbatch --export=expt=${expt},tp_id=${tp_id} batch/estimate_lineage_freqs_recursively.sbatch
	sleep 0.2
done

# UIR 2018
expt=UIR_2018
for tp_id in {1..19}
do
	echo \> ${expt} tp_id=${tp_id}
	sbatch --export=expt=${expt},tp_id=${tp_id} batch/estimate_lineage_freqs_recursively.sbatch
	sleep 0.2
done

# UIR 2019
expt=UIR_2019
for tp_id in {1..15}
do
	echo \> ${expt} tp_id=${tp_id}
	sbatch --export=expt=${expt},tp_id=${tp_id} batch/estimate_lineage_freqs_recursively.sbatch
	sleep 0.2
done

# 4) Bind results together, per experiment
module load R/3.3.3-fasrc01

expt=UCP_2018
Rscript bind_results.r ${outdir}/freqs_${expt}.txt ${scratch}/freqs_${expt}_{1..21}.txt
Rscript bind_results.r ${outdir}/fit_stats_${expt}.txt ${scratch}/fit_stats_${expt}_{1..21}.txt

expt=UCP_2019
Rscript bind_results.r ${outdir}/freqs_${expt}.txt ${scratch}/freqs_${expt}_{1..15}.txt
Rscript bind_results.r ${outdir}/fit_stats_${expt}.txt ${scratch}/fit_stats_${expt}_{1..15}.txt

expt=UIR_2018
Rscript bind_results.r ${outdir}/freqs_${expt}.txt ${scratch}/freqs_${expt}_{1..19}.txt
Rscript bind_results.r ${outdir}/fit_stats_${expt}.txt ${scratch}/fit_stats_${expt}_{1..19}.txt

expt=UIR_2019
Rscript bind_results.r ${outdir}/freqs_${expt}.txt ${scratch}/freqs_${expt}_{1..15}.txt
Rscript bind_results.r ${outdir}/fit_stats_${expt}.txt ${scratch}/fit_stats_${expt}_{1..15}.txt

