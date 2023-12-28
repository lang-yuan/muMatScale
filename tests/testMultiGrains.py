#!/usr/bin/env python
import sys
import subprocess
import os

print("Test CA multiple grains...")

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

expected_size =  7.14e-5
for line in lines:
  if line.count(b'Average'):
    print(line)
    words=line.split()
    ngrains = eval(words[0])
    if (ngrains != 336):
      print("Unexpected number of grains: Expected 336, got {}".format(ngrains))
      sys.exit(1)
    avg = eval(words[9])
    if abs(avg-expected_size)>1.e-6:
      print("Unexpected average grain size: {}".format(avg))
      sys.exit(1)

sys.exit(0)
