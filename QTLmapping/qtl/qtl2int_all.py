k = 500
shell_script = 'qtl2int_one.sh'
f = open('qtl2int_all.sh', 'w+')

for i in range(1,k+1):
  line = 'qsub -o log_output.txt -e log_error.txt -v nparts=' + str(k) + ',part=' + str(i) + " " + shell_script + "\n"
  f.write(line)
