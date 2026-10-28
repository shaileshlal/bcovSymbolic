from experiment import *
from fractions import Fraction as Q

def full_response(g,N):
 K,m,D=dimensions(g);vals,k,_=frame('mum',N+10);H=[[Q(0)]*D for _ in range(N)]
 for r in range(D):
  gw=q_coefficients(evaluate_poly(basis_poly(g,r),vals)*k**(g-1),N)
  dn=[Q(0)]*(N+1)
  for d in range(1,N+1):
   dn[d]=gw[d]-sum((Q(a)**(2*g-3)*dn[d//a] for a in range(2,d+1) if d%a==0),Q(0))
   H[d-1][r]=dn[d]
 return H

def main():
 for g,N in [(4,10),(5,10),(6,10),(7,10),(8,10),(9,9),(10,9),(11,9),(12,9)]:
  t=time.monotonic();K,m,D=dimensions(g);H=full_response(g,N)
  B=[solve(H[:D],[Q(i==j) for i in range(D)]) for j in range(D)]
  A=[[sum((H[d][r]*B[j][r] for r in range(D)),Q(0)) for j in range(D)] for d in range(N)]
  assert all(A[i][j]==(i==j) for i in range(D) for j in range(D))
  print('g',g,'D',D,'N',N,'elapsed',round(time.monotonic()-t,3))
  for j in range(D):
   dd=[(d+1, A[d][j].denominator) for d in range(D,N)]
   print(' unit n',j+1,'higher denominator:', dd, 'lcm',reduce(math.lcm,[v for _,v in dd],1))
  print(flush=True)
if __name__=='__main__':main()
