#!/usr/bin/env python
'''
STU analisis for SDFDM model
'''
import numpy as np
import sys
def PI0(m1,m2,LAMBDA=1e+16):
    if m1==m2:
        pizero=0.
    elif m1==-m2:
        pizero=1.0/(16.0*np.pi**2)*( 4.0*m1**2*(-1+np.log(LAMBDA**4/m1**4)) )
    else:
        pizero = 1./(16*np.pi**2)*( (m1-m2)**2*np.log(LAMBDA**4/(m1**2*m2**2))-2*m1*m2\
            +( 2*m1*m2*(m1**2+m2**2)-m1**4-m2**4 )/( m1**2-m2**2 )*np.log(m1**2/m2**2) )
    return pizero

def A(m1,m2,v=246.2,alpha_em=1./128.,LAMBDA=1e+16,):
    aa=1./(alpha_em*v**2)*PI0(m1,m2,LAMBDA)   
    return 2*aa  #Note the factor 2.0 because we are dealing with Majorana fermions

def stu_new(ms,md,ll,tbeta,vev=246.2,sw=np.sqrt(0.23),gg=0.653,LAMBDA=1E16):
    alpha=(gg*sw)**2/(4.0*np.pi)
    yy=ll*np.cos(np.arctan(tbeta))  
    yyp=ll*np.sin(np.arctan(tbeta)) 
    gamma=(yyp+yy)/2.0
    beta=(yyp-yy)/2.0
    
    #Mass matrix from arXiv:1411.1335, eq(12). 
    Mmass=[[-md,              0,               -beta*vev ],
           [0.0,              md,              -gamma*vev], 
           [-beta*vev, -gamma*vev, ms]]
    #Mass matrix from PRD by ERAMO, eq(6) but with md->-md, ms->-2ms in 
    #order to be compatible with arXiv:1411.1335, eq(12).   
    #These changes were also implemented in the T parameter formulae.

    eigsys=np.linalg.eigh(Mmass)
    m=eigsys[0]
    V=eigsys[1]
    
    TT1=0
    TT2=0
    for i in range(3):
        TT1=TT1+(V[0,i])**2*A(-md,m[i],vev,alpha,LAMBDA)+(V[1,i])**2*A(-md,-m[i],vev,alpha,LAMBDA)
        for j in range(3):
            TT2=TT2+(V[0,i]*V[1,j]+V[1,i]*V[0,j])**2*A(m[i],-m[j],vev,alpha,LAMBDA)

    TT=TT1-0.5*TT2

    return TT

vstu_new=np.vectorize(stu_new,excluded={'vev':246.2,'sw':np.sqrt(0.23),'gg':0.653,'LAMBDA':1E16})