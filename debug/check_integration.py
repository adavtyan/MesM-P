# Author: Aram Davtyan

import sys

class Frame:
  type = []
  dphi_m = []
  dphi_b = []
  phi_m = []
  phi_b = []

  def __init__(self):
    self.type = []
    self.dphi_m = []
    self.dphi_b = []
    self.phi_m = []
    self.phi_b = []

filename = sys.argv[1]

sname = ["TIMESTEP", "NUMBER OF ATOMS", "BOX BOUNDS", "ATOMS"]

file = open(filename)

iframe = -1
step = 0
natom = 0
A = []
frames = []
file_stat = "NONE"
for line in file:
  if line[:5]=="ITEM:":
    for s in sname:
      if line[6:6+len(s)]==s:
        file_stat = s
        break
  else:
    if file_stat=="TIMESTEP":
      step = int(line)
      A = []
      iframe += 1
      frames.append(Frame())
    elif file_stat=="NUMBER OF ATOMS":
      natom = int(line)
    elif file_stat=="BOX BOUNDS":
      l = line.split()
      A.append([float(l[0]), float(l[1])])
    elif file_stat=="ATOMS":
      l = line.split()
      id = int(l[0])
      type = int(l[1])
      x = float(l[2])
      y = float(l[3])
      z = float(l[4])
      dphi_m = float(l[11])
      dphi_b = float(l[12])
      phi_m = float(l[13])
      phi_b = float(l[14])

      frames[iframe].type.append(type)
      frames[iframe].phi_m.append(phi_m)
      frames[iframe].phi_b.append(phi_b)
      frames[iframe].dphi_m.append(dphi_m)
      frames[iframe].dphi_b.append(dphi_b)

n_frame = iframe + 1

file.close()

print "n_frame", n_frame
print "natom", natom

delta_phi_m = []
delta_phi_b = []
for i in range(n_frame-1):
  delta_phi_m.append([])
  delta_phi_b.append([])
  for j in range(natom):
    delta_phi_m[i].append(frames[i+1].phi_m[j] - frames[i].phi_m[j] - 0.5*(frames[i+1].dphi_m[j] + frames[i].dphi_m[j]))
    delta_phi_b[i].append(frames[i+1].phi_b[j] - frames[i].phi_b[j] - 0.5*(frames[i+1].dphi_b[j] + frames[i].dphi_b[j]))

imax_m = []
imax_b = []
max_m = []
max_b = []
for i in range(len(delta_phi_m)):
  imax_m.append(-1)
  imax_b.append(-1)
  max_m.append(0.0)
  max_b.append(0.0)
  for j in range(len(delta_phi_m[i])):
    if abs(delta_phi_m[i][j])>max_m[i]:
      max_m[i] = abs(delta_phi_m[i][j])
      imax_m[i] = j+1
      max_b[i] = abs(delta_phi_b[i][j])
      imax_b[i] = j+1

print imax_m
print max_m
print imax_b
print max_b

sum_m = 0.0
sum_b = 0.0
for j in range(natom):
  sum_m += frames[0].phi_m[j] + frames[0].dphi_m[j]
  sum_b += frames[0].phi_b[j] + frames[0].dphi_b[j]
print "sum_m", sum_m
print "sum_b", sum_b
