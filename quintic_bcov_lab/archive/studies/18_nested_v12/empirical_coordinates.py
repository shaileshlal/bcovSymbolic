"""Purely data-driven PT normalization. No DT, GV, modularity, wall crossing.
Phase and quadratic envelope are previously fitted empirical features (v0.9).
"""
from fractions import Fraction as F

PHASE={
0:(1,1,-1,-1),
1:(1,-1,-1,1),
2:(1,1,-1,-1),
3:(1,-1,-1,1),
4:(1,1,-1,-1),
}
C=(10,6,4,4,6)

def W(s,t,k):return 5*k*k+(2*s-5)*k+C[s]-2*t

def K(s,t):return t-1 if s==0 else t+3-s

def group(s):return 0 if s==0 else (1 if s in (1,4) else 2)

def sigma(s,k):return PHASE[s][k%4]

def terminal_t(s,r):return 2*r-1 if s==1 else 2*r+3

def normalized_pt(pt,d,t,B):
 s=d%5;k=d//5;return F(pt,sigma(s,k)*W(s,t,k)) - B

def Z_from_PT(pt,d,t,B):
 s=d%5;k=d//5;return (-1 if t%2 else 1)*(F(pt,sigma(s,k))-W(s,t,k)*B)
