from pathlib import Path
import json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
H=Path(__file__).resolve().parent;M=json.loads((H/'parameters.json').read_text(encoding='utf-8'));D={'names':M['names'],'families':{f:[q for q in M['conditions'] if q['family']==f] for f in ['AGST','AGSlash','GMM','Laplace']}};S=H/'generated';S.mkdir(exist_ok=True)
families=['AGST','AGSlash','GMM','Laplace'];titles=["(a) Gaussian + Student's t",'(b) Gaussian + Slash','(c) Gaussian mixture','(d) Gaussian + Laplace outliers']
names=[n.replace('alpha',r'$\alpha$') for n in D['names']]
colors=['b','k','g','m','c','b','g',(0,.35,.65),(.85,.45,0),(.2,.2,.2),(.45,.2,.65),(.6,.3,.15),'r','m']
styles=['--','--','--','--','--','-.','-.','-','-','-.','--','-.','-','-'];markers=['^','s','d','+','x',None,None,'v','p','>','*','o','o','s']
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':9,'axes.linewidth':.7,'pdf.fonttype':42,'ps.fonttype':42,'savefig.dpi':300})
fig,axs=plt.subplots(2,2,figsize=(7.15,6.35));fig.subplots_adjust(left=.09,right=.988,bottom=.07,top=.825,hspace=.31,wspace=.24)
arrays={}
for family,title,ax in zip(families,titles,axs.flat):
 y=np.array([q['rmse'] for q in D['families'][family]]);arrays[family]=y
 for j,name in enumerate(names):ax.plot(range(1,9),y[:,j],color=colors[j],linestyle=styles[j],marker=markers[j],lw=1.45 if j>=12 else .9,ms=3.8,markerfacecolor='white',label=name,zorder=5 if j>=12 else 2)
 ax.set_xticks(range(1,9),[f'C{i}' for i in range(1,9)]);ax.set_xlim(.8,8.2);ax.set_ylim(bottom=0);ax.set_xlabel('Noise condition');ax.set_ylabel('RMSE');ax.set_title(title,fontsize=9);ax.grid(alpha=.2);ax.set_axisbelow(True)
 if family in ['AGST','Laplace']:
  zoom=ax.inset_axes([.14,.57,.43,.35]);hi=18 if family=='AGST' else 10
  for j in range(14):zoom.plot(range(1,9),y[:,j],color=colors[j],linestyle=styles[j],marker=markers[j],lw=.95 if j>=12 else .65,ms=2,markerfacecolor='white',zorder=5 if j>=12 else 2)
  zoom.set_xlim(.8,8.2);zoom.set_ylim(0,hi);zoom.set_xticks([1,4,8],['C1','C4','C8']);zoom.set_yticks([0,hi/2,hi]);zoom.tick_params(labelsize=6,pad=1,length=2);zoom.set_title(f'RMSE 0-{hi}',fontsize=6,pad=3);zoom.grid(alpha=.15)
fig.legend(*axs[0,0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.52,.995),ncol=5,fontsize=7.5,columnspacing=.9,handlelength=2.2,frameon=True)
fig.savefig(S/'Synthetic_Condition_RMSE.pdf');fig.savefig(S/'Synthetic_Condition_RMSE.png');np.savez(S/'figure_values.npz',**arrays);plt.close(fig)
