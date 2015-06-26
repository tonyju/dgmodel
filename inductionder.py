'''induction generator'''
from sympy import *
import numpy as np
#dynamic Rr induction generator
uare = Symbol("x1",real=True)
uaim = Symbol("x2",real=True)
ubre = Symbol("x3",real=True)
ubim = Symbol("x4",real=True)
ucre = Symbol("x5",real=True)
ucim = Symbol("x6",real=True)
iare = Symbol("y1",real=True)
iaim = Symbol("y2",real=True)
ibre = Symbol("y3",real=True)
ibim = Symbol("y4",real=True)
icre = Symbol("y5",real=True)
icim = Symbol("y6",real=True)
Rr=Symbol("z1",real=True)
psp=Symbol("u1",real=True)
rs=Symbol("u2",real=True)
xs=Symbol("u3",real=True)
xm=Symbol("u4",real=True)
#islimit=Symbol("u10",real=True)
are=Symbol("u5",real=True)
aim=Symbol("u6",real=True)
real1_3=Symbol("u7",real=True)
slip=Symbol("u8",real=True)
Xr=Symbol("u9",real=True)


T=(real1_3)*Matrix([[1,1,1],[1,are+I*aim,are-I*aim],[1,are-I*aim,are+I*aim]])
iseq=T*Matrix([iare+I*iaim,ibre+I*ibim,icre+I*icim])
useq=T*Matrix([uare+I*uaim,ubre+I*ubim,ucre+I*ucim])
id0=re(iseq[0])
iq0=im(iseq[0])
ud0=re(useq[0])
uq0=im(useq[0])
id1=re(iseq[1])
iq1=im(iseq[1])
ud1=re(useq[1])
uq1=im(useq[1])
id2=re(iseq[2])
iq2=im(iseq[2])
ud2=re(useq[2])
uq2=im(useq[2])
r0=0
r1=re(rs+I*xs+I*xm*(Rr+I*slip*Xr)/(I*xm+Rr+I*slip*Xr))
x0=0
x1=im(rs+I*xs+I*xm*(Rr+I*slip*Xr)/(I*xm+Rr+I*slip*Xr))
r2=re(rs+I*xs+I*xm*(Rr+I*(2-slip)*Xr)/(I*xm+Rr+I*(2-slip)*Xr))
x2=im(rs+I*xs+I*xm*(Rr+I*(2-slip)*Xr)/(I*xm+Rr+I*(2-slip)*Xr))

f1=re((uare+I*uaim)*(iare-I*iaim)+(ubre+I*ubim)*(ibre-I*ibim)+(ucre+I*ucim)*(icre-I*icim))+(id0**2+iq0**2)*r0+(id1**2+iq1**2)*r1+(id2**2+iq2**2)*r2-psp
f2=re((id0+I*iq0))
f3=im((id0+I*iq0))
f4=re((id2+I*iq2)*(r2+I*x2)+ud2+I*uq2)
f5=im((id2+I*iq2)*(r2+I*x2)+ud2+I*uq2)
f6=re((id1+I*iq1)*(r1+I*x1)+(ud1+I*uq1))
f7=im((id1+I*iq1)*(r1+I*x1)+(ud1+I*uq1))
f=Matrix([f1,f2,f3,f4,f5,f6,f7])
print f1
print f2
print f3
print f4
print f5
print f6
print f7
print latex(f1)
print latex(f2)
print latex(f3)
print latex(f4)
print latex(f5)
print latex(f6)
print latex(f7)


#print latex(f)
var=Matrix([uare,uaim,ubre,ubim,ucre,ucim,iare,iaim,ibre,ibim,icre,icim,Rr])
J=f.jacobian(var)
for i in range(7):
    for j in range(13):
        #pass
        print (i,j),latex(J[i,j])
