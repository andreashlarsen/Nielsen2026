import numpy as np
import matplotlib.pyplot as plt
from smooth import *

PLOT_ALL = 0
PLOT_AV = 1
IDP = 0

#skip_first,skip_last=20,15
#skip_first,skip_last=5,15
skip_first,skip_last=0,0
if IDP:
    pep = ['idp_hec23','idp_b2','idp_atr','idp_kor2','idp_d22','idp_d1']
    rep = [3,3,3,4,3,3]
#    rep = [0,0,0,4,0,0]
    label = ['Hecate','B2AR','AT1R','KOR','D2R','D1R']
else:
    pep = ['hecate','b2','atr','kor2','d22','d1','div4a','ashuffle','gshuffle']
    rep = [3,3,3,3,3,3,0,0,0]
    label = ['Hecate','B2AR','AT1R','kOR','D2R','D1R','DIV4a','Ashuffle','Gshuffle']
color = ['green','blue','cyan','orange','magenta','red','cyan','orange','grey','grey','grey','grey','grey']

plt.rcParams.update({'font.size': 12})

N_bulk = 25
for i in range(len(pep)):

    sum_E = 0
    list_E = []
    matrix_E = []
    
    for r in range(rep[i]):
        if IDP:
            #filename = '%s/CG_rep%d/umb_%s_bE_rep%d/bsResult.xvg' % (pep[i],r+1,pep[i],r+1)
            filename = '%s/CG_rep%d/umb_%s_bE_rep%d/profile.xvg' % (pep[i],r+1,pep[i],r+1)
        else:
            #filename = '%s/rep%d/umb_%s_bE_rep%d/bsResult.xvg' % (pep[i],r+1,pep[i],r+1)
            filename = '%s/rep%d/umb_%s_bE_rep%d/profile.xvg' % (pep[i],r+1,pep[i],r+1)
        #d,E,dE = np.genfromtxt(filename,skip_header=18+skip_first,skip_footer=skip_last,unpack=True)
        d,E = np.genfromtxt(filename,skip_header=17+skip_first,skip_footer=skip_last,unpack=True)
        E -= np.mean(E[-N_bulk:])
        E = smooth(E,5)
        #dE = smooth(dE,5)
        
        sum_E += E
        list_E.append(E)

        #err_bulk = np.sqrt(np.sum(dE[-N_bulk:]**2))/np.sqrt(N_bulk)
        idx = np.argmin(E)
        d_min = d[idx]
        #E_min,err_min = E[idx],dE[idx]
        E_min = E[idx]
        #dG = np.sqrt(err_bulk**2 + err_min**2)
        G = E_min
        #print('%-12s, rep %d: DG =  %1.1f +/- %1.1f' % (label[i],r+1,G,dG))
        print('%-12s, rep %d: DG =  %1.1f' % (label[i],r+1,G))
        
        if PLOT_ALL:
            #plt.plot(d,E,color=color[i],label=r'%s #%d, $\Delta$G = %1.1f $\pm$ %1.1f, d_min = %1.2f' % (label[i],r,G,dG,d_min))
            plt.plot(d,E,color=color[i],label=r'%s #%d, $\Delta$G = %1.1f, d_min = %1.2f' % (label[i],r,G,d_min))
            #plt.fill_between(d,E-dE,E+dE,color=color[i],alpha=0.5)

    if rep[i] > 1:
        mean_E = sum_E/rep[i]
        matrix_E = np.vstack(list_E)
        std_E = matrix_E.std(0)
        err_E = std_E/np.sqrt(rep[i])
        idx = np.argmin(mean_E)
        d_min = d[idx]
        E_min,err_min = mean_E[idx],err_E[idx]
        err_bulk = np.sqrt(np.sum(err_E[-N_bulk:]**2))/np.sqrt(N_bulk)
        dG = np.sqrt(err_bulk**2 + err_min**2)
        G = E_min
        print('%-12s, rep %d: DG =  %1.1f +/- %1.1f, min dist = %1.2f' % (label[i],r+1,G,dG,d_min))
        depth = 4.5/2-d_min     
        if PLOT_AV:
            #plt.plot(d,mean_E,color=color[i],label='%s: %1.1f $\pm$ %1.1f kJ/mol, depth = %1.2f nm' % (label[i],G,dG,depth))
            plt.plot(d,mean_E,color=color[i],label='%s: %1.1f $\pm$ %1.1f kJ/mol' % (label[i],G,dG))
            plt.fill_between(d,mean_E-err_E,mean_E+err_E,color=color[i],alpha=0.5)   

plt.rcParams.update({'font.size': 11})
plt.ylim(-160,5)    
plt.plot(d,d*0,color='grey',linestyle='--')
plt.xlabel('COM distance from helix to membrane (nm)')
plt.ylabel('Potential of mean force [kJ/mol]')
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig('pmf_binding.pdf')
plt.show()
