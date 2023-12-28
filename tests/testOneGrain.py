#!/usr/bin/env python
import sys
import subprocess
import os

print("Test CA single grain...")

nargs=len(sys.argv)

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
for i in range(4,nargs-2):
  mpicmd = mpicmd + " "+sys.argv[i]
exe = sys.argv[nargs-2]
inp = sys.argv[nargs-1]

#run CA code
command = "{} {} {} {} {}".format(mpicmd,exe,"-i","-c",inp)
print("Run command: {}".format(command))

output = subprocess.check_output(command,shell=True)

#analyse standard output
lines=output.split(b'\n')

av = 0.000737875
ng = 1

for line in lines:
  if line.count(b'Average'):
    print(line)
    words=line.split()
    ngrains = eval(words[0])
    if (ngrains != ng):
      print("Unexpected number of grains: Expected {}, got {}".format(ng,ngrains))
      sys.exit(1)
    avg = eval(words[9])
    if abs(avg-av)>1.e-6:
      print("Unexpected average grain size: Expected {}, got {}".format(av,avg))
      sys.exit(1)

sys.exit(0)
