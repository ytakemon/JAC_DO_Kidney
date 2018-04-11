batch_size = 100
N = 6716
shell_script = 'scan1snps_diffint_p.sh'
f = open('scan1snps_diffint_p_ALL.sh', 'w+')


i = 0
while i<N:
    batch = range(i+1, min(i+batch_size, N)+1)
    batch_str = " ".join([str(j) for j in batch])
    i += batch_size
    line = 'qsub -o log_output_scan1snps_diffint_p.txt -e log_error_scan1snps_diffint_p.txt -v col="' + batch_str + '" ' + shell_script + "\n"
    f.write(line)
