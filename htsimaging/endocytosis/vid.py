#!/usr/bin/env python
"""To make the video of the timelapse images."""

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

def make_gif(
    cellcfg=None,
    frames: list=None,
    t_cor: pd.DataFrame=None,
    img_bright=None,
    outd: str=None,
    particle2color: dict=None,
    test: bool=False,
    force: bool=False,
    ):
    """
    Make a .gif file out of frames. 
    """
    if not cellcfg is None:
        if frames is None:
            frames=[np.load(p) for p in sorted(cellcfg['cellframes'])]
        if t_cor is None:
            t_cor=read_table(f"{cellcfg['outp']}/d2filter_stubs.tsv")
        if img_bright is None:
            img_bright=np.load(cellcfg['cellbrightp'])
        if outd is None:    
            outd=f"{cellcfg['outp']}/vid"
        if not particle2color is None:
            if exists(f"{cellcfg['outp']}/d4distance.tsv"):
                df1=read_table(f"{cellcfg['outp']}/d4distance.tsv").drop_duplicates(subset=['particle'])
                from roux.stat.norm import rescale
                from roux.viz.colors import saturate_color
                df1['distance effective from centroid']=rescale(df1['distance effective from centroid'])
                df1.loc[(df1['distance effective from centroid']<=0.5),'distance effective from centroid']=0.5
                df1['color']=df1['distance effective from centroid'].apply(lambda x: saturate_color('limegreen',x))
                particle2color=df1.set_index('particle')['color'].to_dict()
            else:
                logging.warning('distance file (d4distance.tsv) not found')
    makedirs(outd,exist_ok=True)
    gifp=f"{dirname(outd)}/{basenamenoext(outd)}.gif"

    if exists(gifp) and not force:
        return
    if not 'frame min' in t_cor:
        df=t_cor.groupby('particle').agg({'frame':[min,max]})
        df.columns=coltuples2str(df.columns)
        t_cor=t_cor.merge(df,left_on='particle',right_index=True)
    t_cor=t_cor.sort_values(['frame','particle'])
    for framei,frame in enumerate(frames):
        plot_moving_particles(t_cor,img_bright=img_bright,frame=frame,
                            framei=framei,particle2color=particle2color,
                             outd=outd,
                             test=False)
    if not test:    
        plt.close('all')
        com=f"convert -delay 10 -loop 0 {outd}/*.jpeg {gifp}"
        from roux.lib.sys import runbashcmd
        runbashcmd(com)
    return gifp
