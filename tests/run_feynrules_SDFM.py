#1. CHOOSE A BENCHMARK POINT  
#if 1==1:
#    MDF = 110.;MN = 101.;tanb = 10.0;lam = 0.15;v=246.2196; 
def run_feynrules_SDFDM(MDF = 110.,MN = 101.,lu = 0.1,ld = 0.1,v=246.2196,path='../micromegas/SDFDM/main'):
    import pandas as pd
    import numpy as np
    import commands
    o={}
    #2. File to run micrOMEGAS installed in galcen
    
    #lu=lam*np.sin(np.arctan(tanb))
    #ld=lam*np.cos(np.arctan(tanb))
    
    M=np.matrix([[ MN,                -ld*v/np.sqrt(2.),  lu*v/np.sqrt(2.)],
             [ -ld*v/np.sqrt(2.),  0.,                MDF ],
             [ lu*v/np.sqrt(2.),  MDF,               0. ]])

    (Mchi,N)=np.linalg.eig(M)
    
    pd.Series({'MDF':MDF,'MN':MN,'ld':ld,'lu':lu,\
               'N11':N[0,0],'N12':N[0,1],'N13':N[0,2],\
               'N21':N[1,0],'N22':N[1,1],'N23':N[1,2],\
               'N31':N[2,0],'N32':N[2,1],'N33':N[2,2]}).to_csv('mo.dat',sep='\t')

    #3. Run micromegas
    mo=commands.getoutput('%s mo.dat' %path)

    #4. Extrac some of the output (Dependence of the micrOMEGAs vertion)
    o['Full']=mo
    o['Mchi']=Mchi
    o['N']=N
    o['Omega']=eval(mo.split('Omega=')[1].split('\n')[0])
    o['proton_SI']=eval(mo.split('proton  SI')[1].split('[')[0])
    o['proton_SD']=eval(mo.split('proton  SI')[1].split('SD')[1].split('[')[0])
    #o['neutron_SI']=eval(mo.split('proton  SI')[1].split('[')[0])
    #o['neutron_SD']=eval(mo.split('proton  SI')[1].split('SD')[1].split('[')[0])
    o['sigmav']=eval(mo.split('annihilation cross section')[1].split('cm^3/s\n')[0])
    return o