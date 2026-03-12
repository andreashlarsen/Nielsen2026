import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

PLOT_ALL = 0
PLOT_AV  = 1
IDP = 0
HECATE = 0
PLOT_SUMMARY_HECATE=0

peps = ['hecate','b2','d22','d1','ashuffle','gshuffle','atr','div4a','kor2','hecate11','hecate14','hecate17','hecate20','hecate23','hecate26','hecate29','hecate32','hecate35','hecate38','asyn','alps','hecate_idp','kor2_idp','b2_idp','atr_idp','d22_idp','d1_idp']
if HECATE:
    reps = [10,0,0,0,0,0,0,0,5,5,5,5,5,5,5,5,5,5]
    ylim = [-2,2.5]
    names = ['hecate15C (Paris-GUV construct)',r'H8 from $\beta2$AR',r'H8 from D2','Ashuffle','Gshuffle','H8 from ATR','Helix from DivIVA',r'H8 from $\kappa$OR','hecate11','hecate14','hecate17','hecate20','hecate23','hecate26','hecate29','hecate32','hecate35','hecate38']
elif IDP:
    reps = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,3,3,3,3,3]
    ylim = [-3,3]
    names = peps
else:
    #reps = [10,3,3,5,3,5,3,3,0,0,0,0,0,0,0,0,0,0,3,3]
    reps = [10,3,3,3,0,0,3,0,3,0,0,0,0,0,0,0,0,0,0,0,0]
    ylim = [-2.5,5]
    names = peps
    #names = ['Hecate',r'H8 from $\beta2$AR',r'H8 from D2','Ashuffle','Gshuffle','H8 from ATR','Helix from DivIVA',r'H8 from $\kappa$OR','hecate11','hecate14','hecate17','hecate20','hecate23','hecate26','hecate29','aSyn','ALPS']
    names = ['Hecate','B2AR','D2R','D1R','ashuffle','gshuffle','AT1R','div4a','KOR','hecate11','hecate14','hecate17','hecate20','Hecate','hecate26','hecate29','hecate32','hecate35','hecate38','asyn','alps','hecate_idp','kor_idp','b2_idp','atr_idp','d2_idp','d1_idp']
colors = ['green','blue','magenta','red','grey','orange','cyan','brown','orange','orange','brown','cyan','green','chocolate','deepskyblue','black','purple','skyblue','firebrick','pink','blue','cyan','red','magenta','brown','black','grey']

f = open('pmf_curv.txt','w')
f.write('%15s %10s %10s %5s\n' % ('# name','E [kT]','dE [kT]','reps'))

fontsize = 12
plt.rcParams.update({'font.size': fontsize})

for (name,pep,rep,c) in zip(names,peps,reps,colors):
    
    sum_E = 0
    list_E = []
    matrix_E = []
    for i in range(rep):
        K,pmf,pmf_av = np.genfromtxt('%s/rep%d/pmf_100_1000.dat' % (pep,i+1),usecols=[0,1,2],unpack=True)
        idx_0 = np.where(K>0)[0][0]
        pmf_0 = pmf_av[idx_0]
        pmf -= pmf_0
        pmf_av -= pmf_0
        if PLOT_ALL:
            plt.plot(K,pmf_av,color=c,alpha=0.7,zorder=100)
    
        # interpolate
        x = np.linspace(-0.189,0.189,1000)
        f_int = interp1d(K,pmf_av)
        
        sum_E += f_int(x)
        list_E.append(f_int(x))

    if rep > 0:
        mean_E = sum_E/rep
        matrix_E = np.vstack(list_E)
        std_E = matrix_E.std(0)
        err_E = std_E/np.sqrt(rep)
        if PLOT_AV:
            #plt.plot(x,mean_E,color=c,label='%s (n=%d)' % (name,rep))
            plt.plot(x,mean_E,color=c,label='%s' % name)
            plt.fill_between(x,mean_E-err_E,mean_E+err_E,color=c,alpha=0.5)
        idx = np.where(abs(x)<0.1)
        pmf = np.amax(mean_E[idx])-np.amin(mean_E[idx])
        idx_first = idx[0][0]
        idx_last = idx[0][-1]
        err_pmf = np.sqrt(err_E[idx_first]**2 + err_E[idx_last]**2)
        f.write('%15s %10.2f %10.2f %5d\n' % (pep,pmf,err_pmf,rep))

f.close()

xlim = [-0.2,0.2]
plt.plot([0,0],ylim,color='black',alpha=0.1,linestyle='--')
plt.plot(xlim,[0,0],color='black',alpha=0.1,linestyle='--')
plt.fill_between([-0.2,-0.1],[ylim[1],ylim[1]],[ylim[0],ylim[0]],facecolor='grey',edgecolor='grey',alpha=0.2)
plt.fill_between([0.1,0.2],[ylim[1],ylim[1]],[ylim[0],ylim[0]],facecolor='grey',edgecolor='grey',alpha=0.2)
plt.xlim(xlim)
plt.ylim(ylim)
plt.xlabel(r'Curvature [nm$^{-1}$]')
plt.ylabel('Potential of mean force [kT]')

plt.legend(frameon=False,loc='upper right')
plt.savefig('curvature.pdf')
plt.tight_layout()
plt.show()

plt.rcParams['figure.figsize'] = [6,2.5] 
if PLOT_SUMMARY_HECATE and HECATE:
    l = np.linspace(11,38,10)
    E,dE = np.genfromtxt('pmf_curv.txt',skip_header=1,usecols=[1,2],unpack=True)
    plt.errorbar(l,E,yerr=dE,linestyle='none',marker='.',color='black')
    #plt.ylabel('E(K=-0.1) - E(K=0.1)')
    plt.ylabel('Curvature sensing [kJ/mol]')
    plt.xlabel('Hecate length')
    plt.xticks(l)
    plt.ylim(0,2.2)
    plt.tight_layout()
    plt.savefig('hecate_length.pdf',format='pdf')
    plt.show() 
