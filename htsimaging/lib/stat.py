from rohan.global_imports import *
from htsimaging.lib.plot import *

from scipy.spatial import distance
def distance_effective(particle,frame1,frame2,t_cor):
    a=t_cor.loc[((t_cor['particle']==particle) & (t_cor['frame']==frame1)),['x','y']]
    b=t_cor.loc[((t_cor['particle']==particle) & (t_cor['frame']==frame2)),['x','y']]
    return distance.euclidean(a.values, b.values)

def get_distance_travelled(frames,t_cor,out_fh,test=False,force=False):
    ddistancesp=f"{out_fh}_distances.tsv"
    if not exists(ddistancesp) or force:    
        t_cor=read_table(f"{out_fh}.t2.tsv")
        for f1,f2 in zip(list(range(0,t_cor['frame'].max())),
                    list(range(1,t_cor['frame'].max()+1))):
            for p in t_cor['particle'].unique():
                a=t_cor.loc[((t_cor['particle']==p) & (t_cor['frame']==f1)),['x','y']]
                b=t_cor.loc[((t_cor['particle']==p) & (t_cor['frame']==f2)),['x','y']]
                if len(a)!=0 and len(b)!=0:
                    t_cor.loc[((t_cor['particle']==p) & (t_cor['frame']==f2)),'distance']=distance.euclidean(a.values, b.values)
        if test:        
            print(t_cor.head())
        if 'distance' in t_cor:
            t_cor_distance=t_cor.groupby('particle').agg({'distance':sum})

            t_cor_rangeframes=t_cor.groupby('particle').agg({'frame':[min,max]})
            t_cor_rangeframes.columns=coltuples2str(t_cor_rangeframes.columns)
            t_cor_rangeframes['distance effective']=t_cor_rangeframes.apply(lambda x : distance_effective(x.name,x['frame min'],x['frame max'],t_cor) ,axis=1)

            t_cor=t_cor.merge(t_cor_distance,
                        left_on='particle',right_index=True,suffixes=[' delta',' total'])
            t_cor=t_cor.merge(t_cor_rangeframes,
                        left_on='particle',right_index=True)
            t_cor['distance delta']=t_cor['distance delta'].fillna(0)
            t_cor['distance total per frame']=t_cor.apply(lambda x : t_cor.set_index(['frame','particle']).loc[[(i,x['particle']) for i in range(int(x['frame'])+1)],'distance delta'].sum() ,axis=1)
            t_cor['distance effective per frame']=t_cor.apply(lambda x : distance_effective(particle=x['particle'],frame1=x['frame min'],frame2=x['frame'],t_cor=t_cor) ,axis=1)
            t_cor['intensity']=t_cor['mass']            
            to_table(t_cor,ddistancesp)
            ## plot dist distances
            plotp=f"{out_fh}_hist_distances.png"
            plt.figure()
            ax=plt.subplot()
            t_cor[['distance delta','distance total','distance effective']].dropna().hist(ax=ax)
            plt.tight_layout()
            savefig(plotp)   
            ## plot image labeled particles
            fig=plt.figure(figsize=[10,10])
            ax=plt.subplot()
            ax=plot_trajectories(t_cor, image=frames[0],label=True,colorby='frame',cmap='hsv',
                        ax=ax,
                        #params_text={'size':5},
                        )
            plotp=f"{out_fh}_trajectories.png"
            plt.tight_layout()
            savefig(plotp)    
        else:
            t_cor=pd.DataFrame(columns=t_cor.columns)
            to_table(t_cor,ddistancesp)
        return t_cor
