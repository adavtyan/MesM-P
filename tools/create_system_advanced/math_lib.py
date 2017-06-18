# Author: Aram Davtyan

import numpy as np

def vecmat(v, m):
  ans = np.zeros(3)

  ans[0] = v[0]*m[0][0] + v[1]*m[1][0] + v[2]*m[2][0];
  ans[1] = v[0]*m[0][1] + v[1]*m[1][1] + v[2]*m[2][1];
  ans[2] = v[0]*m[0][2] + v[1]*m[1][2] + v[2]*m[2][2];

  return ans

def quat_to_mat(quat):
  mat = np.zeros((3,3))

  w2 = quat[0]*quat[0]
  i2 = quat[1]*quat[1]
  j2 = quat[2]*quat[2]
  k2 = quat[3]*quat[3]
  twoij = 2.0*quat[1]*quat[2]
  twoik = 2.0*quat[1]*quat[3]
  twojk = 2.0*quat[2]*quat[3]
  twoiw = 2.0*quat[1]*quat[0]
  twojw = 2.0*quat[2]*quat[0]
  twokw = 2.0*quat[3]*quat[0]

  mat[0][0] = w2+i2-j2-k2;
  mat[0][1] = twoij-twokw;
  mat[0][2] = twojw+twoik;

  mat[1][0] = twoij+twokw;
  mat[1][1] = w2-i2+j2-k2;
  mat[1][2] = twojk-twoiw;

  mat[2][0] = twoik-twojw;
  mat[2][1] = twojk+twoiw;
  mat[2][2] = w2-i2-j2+k2;

  return mat

def quat_to_mat_trans(quat):
  mat = np.zeros((3,3))

  w2 = quat[0]*quat[0]
  i2 = quat[1]*quat[1]
  j2 = quat[2]*quat[2]
  k2 = quat[3]*quat[3]
  twoij = 2.0*quat[1]*quat[2]
  twoik = 2.0*quat[1]*quat[3]
  twojk = 2.0*quat[2]*quat[3]
  twoiw = 2.0*quat[1]*quat[0]
  twojw = 2.0*quat[2]*quat[0]
  twokw = 2.0*quat[3]*quat[0]

  mat[0][0] = w2+i2-j2-k2;
  mat[1][0] = twoij-twokw;
  mat[2][0] = twojw+twoik;

  mat[0][1] = twoij+twokw;
  mat[1][1] = w2-i2+j2-k2;
  mat[2][1] = twojk-twoiw;

  mat[0][2] = twoik-twojw;
  mat[1][2] = twojk+twoiw;
  mat[2][2] = w2-i2-j2+k2;

  return mat
