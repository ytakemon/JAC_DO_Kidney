k = 250
shell_script = 'qtl2blup_one.sh'
f = open('qtl2blup_all.sh', 'w+')

for i in range(1,k+1):
  line = 'qsub -o log_blup_output.txt -e log_blup_error.txt -v nparts=' + str(k) + ',part=' + str(i) + " " + shell_script + "\n"
  f.write(line)
