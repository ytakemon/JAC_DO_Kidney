batch_size = 500
N = 6716 # or 6145
shell_script = 'scan_one_p.sh'
f = open('scan_all_p.sh', 'w+')


i = 0
while i<N:
  batch = range(i+1, min(i+batch_size, N)+1)
  batch_str = " ".join([str(j) for j in batch])
  i += batch_size

  line = 'qsub -o log_output.txt -e log_error.txt -v col="' + batch_str + '" ' + shell_script + "\n"
  f.write(line)
