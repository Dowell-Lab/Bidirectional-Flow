----------2022-06-22-----------

Made SAF files for counting bidirs:

awk '{FS=OFS="\t"}; {print $1":"$2+1"-"$3-1,$1,$2+1,$3-1,"+",$4}' mm10_qc4_gc50perc_tss50perc_MASTER_tfit_dreg.bed \
    > mm10_qc4_gc50perc_tss50perc_MASTER_tfit_dreg.saf
awk '{FS=OFS="\t"}; {print $1":"$2+1"-"$3-1,$1,$2+1,$3-1,"+",$4}' hg38_qc4_gc50perc_tss50perc_MASTER_tfit_dreg.bed \
    > hg38_qc4_gc50perc_tss50perc_MASTER_tfit_dreg.saf

cat saf_header.txt hg38_qc4_gc50perc_tss50perc_MASTER_tfit_dreg.saf > hg38_master_tfit_dreg.saf
cat saf_header.txt mm10_qc4_gc50perc_tss50perc_MASTER_tfit_dreg.saf > mm10_master_tfit_dreg.saf


----------2023-08-22----------

Ru remade the merged bidirectional call files, so remade safs
  -Same commands as above, except the input filename was hg38_tfit_dreg_bidirectionals.bed \
    and the output was hg38_tfit_dreg_bidirectionals.saf (same for mouse except mm10)

Annotations were also redownloaded and different regions extracted, so \
    replaced previous annotation files with those from:
  -Human: /scratch/Shares/dowell/genomes/hg38/ncbi/
  -Mouse: /scratch/Shares/dowell/genomes/mm10/ncbi/

Config files also changed accordingly
