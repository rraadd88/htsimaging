#!/usr/bin/env python
"""Visualizations."""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def plot_properties_cell(
    cellcfg,
    df2,
    cols_colorby,
    colx='x',
    coly='y',
    ):
    """Plot properties of the cell.

    Args:
        cellcfg (_type_): config of a cell.
        df2 (_type_): input dataframe.
        cols_colorby (_type_): columns to color by.
        colx (str, optional): column with the x values. Defaults to 'x'.
        coly (str, optional): column with the y values. Defaults to 'y'.
    """
    from roux.lib.stat.norm import rescale
    from roux.lib.plot.contour import plot_contourf
    ncols=4
    nrows=int(np.ceil(len(cols_colorby)/4))
    fig,axes=plt.subplots(nrows,ncols,
#                           sharex=True,sharey=True,
                          figsize=[ncols*4,nrows*4]
                         )
#     metric_type='max'
    for axi,(colorby,ax) in enumerate(zip(cols_colorby,np.ravel(axes))):
        fig,ax=plot_contourf(df2[colx],df2[coly],
                      rescale(df2[colorby]),
                     ax=ax,
                     fig=fig,
                      cbar=True if ((axi+1) % 4)==0 else False,
                     params_contourf={'cmap':'binary','vmin':0,'vmax':1,'corner_mask':False},)
        ax=df2.plot.scatter(x=colx,y=coly,color='green',
                         s=10,
                         alpha=0.5,
                         ax=ax)   
        ax.contour(np.load(cellcfg['cellbrightp']), [0.5], linewidths=1, linestyles='dashed',colors='cyan')    
        ax.set_title(colorby)
        ax.invert_yaxis()
        ax.axis('equal')
#         ax.set_ylim(ax.get_ylim()[::-1])
#         ax.set_axis_off()
    plt.tight_layout()

def image_locate_particles(
    df1: pd.DataFrame,
    frame,
    img_region,
    annotate_particles: str=False,
    fig=None,
    ax: plt.Axes=None,
    ) -> plt.Axes:
    """Plot image with particles.

    Args:
        df1 (pd.DataFrame): input dataframe.
        frame (_type_): image frame.
        img_region (_type_): regions in the image.
        annotate_particles (str, optional): annotate the paticles or not. Defaults to False.
        fig (_type_, optional): figure object. Defaults to None.
        ax (plt.Axes, optional): subplot object. Defaults to None.

    Returns:
        plt.Axes: _description_
    """
    import trackpy as tp
    fig=plt.figure(figsize=[20,20]) if fig is None else fig
    ax=plt.subplot(111) if ax is None else ax
    from htsimaging.viz.image import image_background
    ax=image_background(img_region=img_region,img=frame,ax=ax)
    ax=tp.annotate(df1, frame,ax=ax)
    if annotate_particles:
        _=df1.apply(lambda x:ax.text(x['x'],x['y'],int(x['particle']),color='lime'),axis=1)
#     ax.grid(False)
    return ax

def image_trajectories(
    dtraj: pd.DataFrame,
    img_gfp=None,
    img_bright=None,
    label: bool=True,
    fig=None,
    ax: plt.Axes=None,
    ) -> plt.Axes:
    """Plot trajectories.

    Args:
        dtraj (pd.DataFrame): input dataframe with the trajectories.
        img_gfp (_type_, optional): channel image e.g. GFP. Defaults to None.
        img_bright (_type_, optional): segmentation image e.g. bright field. Defaults to None.
        label (bool, optional): label. Defaults to True.
        fig (_type_, optional): figure object. Defaults to None.
        ax (plt.Axes, optional): subplot object. Defaults to None.

    Returns:
        plt.Axes: subplot
    """
    import trackpy as tp    
    fig=plt.figure(figsize=[20,20]) if fig is None else fig
    ax=plt.subplot(111) if ax is None else ax
    from htsimaging.viz.image import image_background
    ax=image_background(
        img_region=img_bright,img=img_gfp,
        cmap='gfp',
        ax=ax)
    ax = tp.plot_traj(dtraj,label=False,ax=ax,lw=2,plot_style={'color':'k'})
    dtrajagg=dtraj.groupby('particle').agg({c:np.median for c in ['x','y']}).reset_index()
#     dtrajagg.columns=coltuples2str(dtrajagg.columns)
#     print(dtrajagg.iloc[0,:])
    if label:
        dtrajagg.reset_index().apply(lambda x: ax.text(x['x'],x['y'],int(x['particle']),
                                    color='magenta'),
                  axis=1)
    ax.grid(False)
    return ax

def plot_moving_particles(
    t_cor: pd.DataFrame,
    img_bright=None,
    frame=None,
    framei: int=0,
    particle2color=None,
    test: bool=False,
    outd: str=None,
    ):
    """Plot moving particles.

    Args:
        t_cor (pd.DataFrame): input table
        img_bright (_type_, optional): segmentation raw image (e.g. bright field). Defaults to None.
        frame (_type_, optional): image frame. Defaults to None.
        framei (int, optional): image frame index. Defaults to 0.
        particle2color (_type_, optional): particle-wise colors. Defaults to None.
        test (bool, optional): test-mode. Defaults to False.
        outd (str, optional): path to the output directory. Defaults to None.
    """
    plotp=f'{outd}/{framei:03d}.jpeg'
    plt.figure(figsize=[5,5])
    ax=plt.subplot(111)
    if not ((img_bright is None) and (frame is None)):
        from htsimaging.viz.image import image_background
        ax=image_background(img_region=img_bright,img=frame,ax=ax)
    ax.text(ax.get_xlim()[1],ax.get_ylim()[1],
            f"frame#{framei:03d}",ha='right',
            va='bottom',
            size='20',
            color='k')
    for p in t_cor.loc[(t_cor['frame']==framei),'particle'].unique():
        framei_min=t_cor['frame min'].unique()[0]
        df=t_cor.loc[((t_cor['particle']==p) \
                      & (t_cor['frame'].isin(range(framei_min,framei+1)))\
                      & (t_cor['x'].between(0,frame.shape[0]))\
                      & (t_cor['y'].between(0,frame.shape[1]))),:].sort_values(by=['frame','particle'])
        if len(df)!=0:
            if framei!=0:    
                df.plot.line(x='x',y='y',lw=4,
                    c='limegreen' if particle2color is None else particle2color[p],
                    legend=False,ax=ax)
            else:
                df.plot.scatter(x='x',y='y',s=4,
                    c='limegreen' if particle2color is None else particle2color[p],
                    legend=False,ax=ax)
            if test:
                d=df.head(1).T.to_dict()
                d=d[list(d.keys())[0]]
                ax.text(x=d['x'],y=d['y'],s=int(d['particle']),
                    color='magenta',
                    )
    ax.set_xlim(0,frame.shape[1])
    if not test:
        plt.axis('off')
        makedirs(dirname(plotp),exist_ok=True)
        plt.savefig(plotp)
        plt.clf();plt.close();
