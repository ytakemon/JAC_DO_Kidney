samples= "$(cat samples.txt)"
church_dir="/projects/churchill-lab/data/JAC/DO_crosssectional/kidney/JAC_DO_Kidney_RNASeq/fastq"


for sample in ${samples[*]}
do
  cp ${church_dir}/${sample}* ./fastq_files/
done

cd fastq_files
md5sum * > md5check.txt
