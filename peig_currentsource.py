'''power electronic interfaced generator'''
from sympy import *
import numpy as np
from scipy import arctan


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
idc=Symbol("z1",real=True)
psp=Symbol("u1",real=True)
XC=Symbol("u2",real=True)
are=Symbol("u3",real=True)
aim=Symbol("u4",real=True)
real1_3=Symbol("u5",real=True)
cosalphasp=Symbol("u6",real=True)
sinalphasp=Symbol("u7",real=True)
sqrt2_pi=Symbol("u8",real=True)
real1_2=Symbol("u9",real=True)
realsqrt3_2=Symbol("u10",real=True)

uab_complex=(uare+I*uaim)-(ubre+I*ubim)
ubc_complex=(ubre+I*ubim)-(ucre+I*ucim)
uac_complex=(uare+I*uaim)-(ucre+I*ucim)
uba_complex=-(uare+I*uaim)+(ubre+I*ubim)
ucb_complex=-(ubre+I*ubim)+(ucre+I*ucim)
uca_complex=-(uare+I*uaim)+(ucre+I*ucim)

uabmag=sqrt(re(uab_complex)**2+im(uab_complex)**2)
ubcmag=sqrt(re(ubc_complex)**2+im(ubc_complex)**2)
uacmag=sqrt(re(uac_complex)**2+im(uac_complex)**2)
ubamag=uabmag
ucbmag=ubcmag
ucamag=uacmag
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
#cos(c1+)=cos(A+pi/3)
#cos(alpha)=cos(alpha-2 pi /3)
#cos(pi/3)=cos(5pi/3)
cosuA=uare/sqrt(uare**2+uaim**2)
sinuA=uaim/sqrt(uare**2+uaim**2)
cosuB=ubre/sqrt(ubre**2+ubim**2)
sinuB=ubim/sqrt(ubre**2+ubim**2)
cosuC=ucre/sqrt(ucre**2+ucim**2)
sinuC=ucim/sqrt(ucre**2+ucim**2)
cosc1pos=cosuA*real1_2-realsqrt3_2*sinuA
sinc1pos=sinuA*real1_2+realsqrt3_2*cosuA
#c2+=C-pi/3
cosc2pos=cosuC*real1_2+realsqrt3_2*sinuC
sinc2pos=sinuC*real1_2-realsqrt3_2*cosuC
#c3+=B-pi/3
cosc3pos=cosuB*real1_2+realsqrt3_2*sinuB
sinc3pos=sinuB*real1_2-realsqrt3_2*cosuB

cosc3=re(uca_complex)/sqrt(re(uca_complex)**2+im(uca_complex)**2)#ca
sinc3=im(uca_complex)/sqrt(re(uca_complex)**2+im(uca_complex)**2)#ca

cosc2=re(ucb_complex)/sqrt(re(ucb_complex)**2+im(ucb_complex)**2)#cb
sinc2=im(ucb_complex)/sqrt(re(ucb_complex)**2+im(ucb_complex)**2)#cb

cosc1=re(uab_complex)/sqrt(re(uab_complex)**2+im(uab_complex)**2)#ab
sinc1=im(uab_complex)/sqrt(re(uab_complex)**2+im(uab_complex)**2)#ca

cosia=iare/sqrt(iare**2+iaim**2)
sinia=iaim/sqrt(iare**2+iaim**2)

vdc=sqrt2_pi*(uacmag*((-cosc1pos*cosc3-sinc1pos*sinc3)*cosalphasp+(sinc1pos*cosc3-cosc1pos*sinc3)*sinalphasp)
             -uacmag*((-cosc2pos*cosc3-sinc2pos*sinc3)*cosalphasp+(sinc2pos*cosc3-cosc2pos*sinc3)*sinalphasp) 
            +ubamag*((cosc2pos*cosc1+sinc2pos*sinc1)*cosalphasp-(sinc2pos*cosc1-cosc2pos*sinc1)*sinalphasp)
            -ubamag*((cosc3pos*cosc1+sinc3pos*sinc1)*cosalphasp-(sinc3pos*cosc1-cosc3pos*sinc1)*sinalphasp)
            +ubcmag*((cosc3pos*cosc2+sinc3pos*sinc2)*cosalphasp-(sinc3pos*cosc2-cosc3pos*sinc2)*sinalphasp)
            -ubcmag*((-cosc1pos*cosc2-sinc1pos*sinc2)*cosalphasp+(sinc1pos*cosc2-cosc1pos*sinc2)*sinalphasp) )
f1=vdc*idc-psp
#vlim or ilim equation
#when limit is update,f2=sqrt(iare**2+ia**2)-islimit
#when limit vflimit is update, 
f2=re((id0+I*iq0))
f3=im((id0+I*iq0))
f4=re((id2+I*iq2))
f5=im((id2+I*iq2))
f6=cosia-(cosc1pos*cosc1+sinc1pos*sinc1)*cosalphasp+(sinc1pos*cosc1-cosc1pos*sinc1)*sinalphasp
f7=sinia-(sinc1pos*cosc1-cosc1pos*sinc1)*cosalphasp-(cosc1pos*cosc1+sinc1pos*sinc1)*sinalphasp
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
var=Matrix([uare,uaim,ubre,ubim,ucre,ucim,iare,iaim,ibre,ibim,icre,icim])
J=f.jacobian(var)
for i in range(6):
    for j in range(12):
        #pass
        print (i,j),latex(J[i,j])
