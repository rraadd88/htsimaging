from roux.global_imports import *
#from htsimaging.lib.plot import *

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

def get_slope(df,ds):
    return sc.stats.linregress(df.loc[ds.index,'frame'],df.loc[ds.index,'distance effective from centroid per frame']).slope
def get_inflection_point(df,threshold_slope=0.25):
    label=df.name
    df['slope distance effective from centroid versus frame']=df.rolling(6)['y'].apply(lambda x: get_slope(df,x),raw=False)
    if any(df['slope distance effective from centroid versus frame']>threshold_slope) and not all(df['slope distance effective from centroid versus frame']>threshold_slope):
        inflection_point=df.set_index('frame')['slope distance effective from centroid versus frame'].idxmax()
    else:
        inflection_point=None
    df['inflection point']=inflection_point
    return df


def get_distance_from_centroid(df1,center=[75,75]):
    # get distance from center
#     from scipy.spatial import distance
    if 'distance effective from centroid' in df1:
        df1=df1.drop(df1.filter(like='distance effective from centroid',axis=1).columns.tolist(),axis=1)
    df1['distance effective from centroid per frame']=df1.apply(lambda x: distance.euclidean(center,[x['x'],x['y']]),axis=1)
    df=df1.groupby('particle',as_index=False).agg({'distance effective from centroid per frame':[np.min,np.max]})
    df.columns=coltuples2str(df.columns)
    df=df.rename(columns={c:c if not 'per frame a' in c else c.replace('per frame a','') for c in df})
#     print(df.columns)
    df['distance effective from centroid']=df['distance effective from centroid max']-df['distance effective from centroid min']
    return df1.merge(df,on='particle',how='left')
