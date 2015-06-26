'''synchronous generator'''
from sympy import *
import numpy as np
#pv model
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
vfre = Symbol("z1",real=True)
vfim = Symbol("z2",real=True)
psp=Symbol("u1",real=True)

#qsp=Symbol("u2",real=True)
r0=Symbol("u2",real=True)
x0=Symbol("u3",real=True)
rneg=Symbol("u4",real=True)
xneg=Symbol("u5",real=True)
r2=rneg
x2=xneg
rpos=Symbol("u6",real=True)
xpos=Symbol("u7",real=True)
r1=rpos
x1=xpos
#islimit=Symbol("u10",real=True)
are=Symbol("u8",real=True)
aim=Symbol("u9",real=True)
real1_3=Symbol("u10",real=True)
qsp=Symbol("u11",real=True)


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
f1=re((uare+I*uaim)*(iare-I*iaim)+(ubre+I*ubim)*(ibre-I*ibim)+(ucre+I*ucim)*(icre-I*icim))+(id0**2+iq0**2)*r0+(id1**2+iq1**2)*r1+(id2**2+iq2**2)*r2-psp
f2=im((uare+I*uaim)*(iare-I*iaim)+(ubre+I*ubim)*(ibre-I*ibim)+(ucre+I*ucim)*(icre-I*icim))-qsp
#vlim or ilim equation
#when limit is update,f2=sqrt(iare**2+ia**2)-islimit
#when limit vflimit is update, 
f3=re((id0+I*iq0)*(r0+I*x0)+ud0+I*uq0)
f4=im((id0+I*iq0)*(r0+I*x0)+ud0+I*uq0)
f5=re((id2+I*iq2)*(r2+I*x2)+ud2+I*uq2)
f6=im((id2+I*iq2)*(r2+I*x2)+ud2+I*uq2)
f7=re((id1+I*iq1)*(r1+I*x1)+(ud1+I*uq1-vfre-I*vfim))
f8=im((id1+I*iq1)*(r1+I*x1)+(ud1+I*uq1-vfre-I*vfim))
f=Matrix([f1,f2,f3,f4,f5,f6,f7,f8])
print f1
print f2
print f3
print f4
print f5
print f6
print f7
print f8
print latex(f1)
print latex(f2)
print latex(f3)
print latex(f4)
print latex(f5)
print latex(f6)
print latex(f7)
print latex(f8)

#print latex(f)
var=Matrix([uare,uaim,ubre,ubim,ucre,ucim,iare,iaim,ibre,ibim,icre,icim,vfre,vfim])
J=f.jacobian(var)
for i in range(8):
    for j in range(14):
        #pass
        print (i,j),latex(J[i,j])
