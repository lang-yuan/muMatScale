#!/usr/bin/env python
import sys
import subprocess
import os

print("Test AM track...")

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

egrains = 1889
eavg = 8.8e-7

flag = 0
for line in lines:
  if line.count(b'Wall'):
    flag = 1
  if line.count(b'Average'):
    flag = 0
  if flag > 0:
    print(line)
  words=line.split()
  if len(words) == 4 and flag>0:
    fs = eval(words[2][:-1])

efs = 87.3
if abs(efs-fs)>0.1:
  print("Unexpected fs of {}, expected {}".format(fs,efs))
  sys.exit(1)

for line in lines:
  if line.count(b'Average'):
    print(line)
    words=line.split()
    ngrains = eval(words[0])
    if (ngrains != egrains):
      print("Unexpected number of grains: Expected {}, got {}".format(egrains,ngrains))
      sys.exit(1)
    avg = eval(words[9])
    if abs(avg-eavg)>1.e-8:
      print("Unexpected average grain size: {}".format(avg))
      sys.exit(1)

sys.exit(0)
