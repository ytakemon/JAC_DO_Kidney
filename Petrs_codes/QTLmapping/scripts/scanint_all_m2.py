batch_size = 1000
N = 22312 # or 6145
shell_script = 'scanint_one_m2.sh'
f = open('scanint_all_m2.sh', 'w+')


i = 0
while i<N:
  batch = range(i+1, min(i+batch_size, N)+1)
  batch_str = " ".join([str(j) for j in batch])
  i += batch_size

  line = 'qsub -o log_output.txt -e log_error.txt -v col="' + batch_str + '" ' + shell_script + "\n"
  f.write(line)
