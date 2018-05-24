batch_size = 500
N = 22312 # or 6145
shell_script = 'scanBestChr_m.sh'
f = open('scanBestChr_all_m.sh', 'w+')


i = 0
while i<N:
    batch = range(i+1, min(i+batch_size, N)+1)
    batch_str = " ".join([str(j) for j in batch])
    i += batch_size
    line = 'qsub -o log_output_scanBestChr_m.txt -e log_error_scanBestChr_m.txt -v col="' + batch_str + '" ' + shell_script + "\n"
    f.write(line)
