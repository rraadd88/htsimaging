from rohan.global_imports import *
from htsimaging.lib.plot import *

from scipy.spatial import distance
def distance_effective(particle,frame1,frame2,t_cor):
    a=t_cor.loc[((t_cor['particle']==particle) & (t_cor['frame']==frame1)),['x','y']]
    b=t_cor.loc[((t_cor['particle']==particle) & (t_cor['frame']==frame2)),['x','y']]
    return distance.euclidean(a.values, b.values)

def get_distance_travelled(t_cor):
    for f1,f2 in zip(list(range(0,t_cor['frame'].max())),
                list(range(1,t_cor['frame'].max()+1))):
        for p in t_cor['particle'].unique():
            a=t_cor.loc[((t_cor['particle']==p) & (t_cor['frame']==f1)),['x','y']]
            b=t_cor.loc[((t_cor['particle']==p) & (t_cor['frame']==f2)),['x','y']]
            if len(a)!=0 and len(b)!=0:
                t_cor.loc[((t_cor['particle']==p) & (t_cor['frame']==f2)),'distance']=distance.euclidean(a.values, b.values)
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
    return t_cor
