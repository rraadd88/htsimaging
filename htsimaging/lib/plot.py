from roux.global_imports import *

def _plot(ax, coords, pos_columns, **plot_style):
    """ This function wraps Axes.plot to make its call signature the same for
    2D and 3D plotting. The y axis is inverted for 2D plots, but not for 3D
    plots.

    Parameters
    ----------
    ax : Axes object
        The axes object on which the plot will be called
    coords : DataFrame
        DataFrame of coordinates that will be plotted
    pos_columns : list of strings
        List of column names in x, y(, z) order.
    plot_style : keyword arguments
        Keyword arguments passed through to the `Axes.plot(...)` method

    Returns
    -------
    Axes object
    """
    if len(pos_columns) == 3:
        return ax.plot(coords[pos_columns[0]], coords[pos_columns[1]],
                       zs=coords[pos_columns[2]], **plot_style)
    elif len(pos_columns) == 2:
        return ax.plot(coords[pos_columns[0]], coords[pos_columns[1]],
                       **plot_style)
def plot_trajectories(traj,image, colorby='particle', mpp=None, label=False,
                        cmap=None, ax=None, t_column=None,
                        pos_columns=None, plot_style={},
                        params_text={'ha':'center','va':'center'}, **kwargs):
    """Plot traces of trajectories for each particle.
    Optionally image it on a frame from the video.

    Parameters
    ----------
    traj : DataFrame
        The DataFrame should include time and spatial coordinate columns.
    colorby : {'particle', 'frame'}, optional
    mpp : float, optional
        Microns per pixel. If omitted, the labels will have units of pixels.
    label : boolean, optional
        Set to True to write particle ID numbers next to trajectories.
    image : ndarray, optional
        Background image, default None
    cmap : colormap, optional
        This is only used in colorby='frame' mode. Default = mpl.cm.winter
    ax : matplotlib axes object, optional
        Defaults to current axes
    t_column : string, optional
        DataFrame column name for time coordinate. Default is 'frame'.
    pos_columns : list of strings, optional
        Dataframe column names for spatial coordinates. Default is ['x', 'y'].
    plot_style : dictionary
        Keyword arguments passed through to the `Axes.plot(...)` command

    Returns
    -------
    Axes object
    
    See Also
    --------
    plot_traj3d : the 3D equivalent of `plot_traj`
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from roux.lib.plot.colors import get_cmap_subset
#     if cmap is None:
#         cmap = get_cmap_subset('binary',vmin=0.15,vmax=0.05)
    if t_column is None:
        t_column = 'frame'
    if pos_columns is None:
        pos_columns = ['x', 'y']
    if len(traj) == 0:
        raise ValueError("DataFrame of trajectories is empty.")
    _plot_style = dict(linewidth=1)
    # Background image
    ax=plt.subplot() if ax is None else ax
    ax.imshow(image,cmap=get_cmap_subset('binary',vmin=0.2,vmax=0))
    # Trajectories
    
    
#     if colorby == 'particle':
#         # Unstack particles into columns.
#         unstacked = traj.set_index(['particle', t_column])[pos_columns].unstack()
#         for i, trajectory in unstacked.iterrows():
#             _plot(ax, mpp*trajectory, pos_columns, **_plot_style)
#     elif colorby == 'frame':
#         # Read http://www.scipy.org/Cookbook/Matplotlib/MulticoloredLine
#         x = traj.set_index([t_column, 'particle'])['x'].unstack()
#         y = traj.set_index([t_column, 'particle'])['y'].unstack()
#         color_numbers = traj[t_column].values/float(traj[t_column].max())
# #         logger.info("Drawing multicolor lines takes awhile. "
# #                     "Come back in a minute.")
#         for particle in x:
#             points = np.array(
#                 [x[particle].values, y[particle].values]).T.reshape(-1, 1, 2)
#             segments = np.concatenate([points[:-1], points[1:]], axis=1)
#             lc = LineCollection(segments, cmap=cmap)
#             lc.set_array(color_numbers)
#             ax.add_collection(lc)
# #             ax.set_xlim(x.apply(np.min).min(), x.apply(np.max).max())
# #             ax.set_ylim(y.apply(np.min).min(), y.apply(np.max).max())
#     else:    
        
#     if label:
#         unstacked = traj.set_index([t_column, 'particle'])[pos_columns].unstack()
#         first_frame = int(traj[t_column].min())
#         coords = unstacked.fillna(method='backfill').stack().loc[first_frame]
#         for particle_id, coord in coords.iterrows():
#             ax.text(*coord.tolist(), s="%d" % particle_id,
#                     **params_text)
    ax.set_xlim(0,image.shape[0])
    ax.set_ylim(0,image.shape[1])
    ax.set_ylim(ax.get_ylim()[::-1])
    return ax    

def plot_properties_cell(cellcfg,df2,cols_colorby,colx='x',coly='y'):
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
    
def dist_signal(img,threshold=None,label_threshold=None,
                params_hist={},
                params_axvline={'color':'r','linestyle':'dashed'},
                ax=None):
    ax=plt.subplot() if ax is None else ax 
    a=np.ravel(img)
    a = a[~np.isnan(a)]
    a=np.extract(a<np.quantile(a,0.9),a)
    _=ax.hist(a,linewidth=0,**params_hist)
    # ax.set_xlim(np.min(a),np.quantile(a,0.9))
    if not threshold is None:
        ax.axvline(threshold,label=label_threshold,**params_axvline)
        ax.legend()
    ax.set_xlabel('signal')
    ax.set_ylabel('density')
    return ax

def image_background(img_region=None,img=None,ax=None,cmap='binary_r',alpha=1):
    if cmap=='gfp':
        from rohan.lib.plot.colors import make_cmap
        cmap=make_cmap(["#000000",'#83f52c'],N=50)
    ax=plt.subplot(111) if ax is None else ax
    if not img is None:
        ax.imshow(img,cmap=cmap,alpha=alpha)
    if not img_region is None:
        ax.contour(img_region, [0.5], linewidths=1, linestyles='dashed',colors='cyan')
    ax.grid(False)
    return ax

def image_locate_particles(df1,frame,img_region,fig=None,ax=None,annotate_particles=False):
    import trackpy as tp
    fig=plt.figure(figsize=[20,20]) if fig is None else fig
    ax=plt.subplot(111) if ax is None else ax
    ax=image_background(img_region=img_region,img=frame,ax=ax)
    ax=tp.annotate(df1, frame,ax=ax)
    if annotate_particles:
        _=df1.apply(lambda x:ax.text(x['x'],x['y'],int(x['particle']),color='lime'),axis=1)
#     ax.grid(False)
    return ax

def image_trajectories(dtraj,img_gfp=None,img_bright=None,label=True,
                       fig=None,ax=None):
    import trackpy as tp    
    fig=plt.figure(figsize=[20,20]) if fig is None else fig
    ax=plt.subplot(111) if ax is None else ax
    ax=image_background(img_region=img_bright,img=img_gfp,ax=ax)
    ax = tp.plot_traj(dtraj,label=False,ax=ax,lw=2,plot_style={'color':'lime'})
    dtrajagg=dtraj.groupby('particle').agg({c:np.median for c in ['x','y']}).reset_index()
#     dtrajagg.columns=coltuples2str(dtrajagg.columns)
#     print(dtrajagg.iloc[0,:])
    if label:
        dtrajagg.reset_index().apply(lambda x: ax.text(x['x'],x['y'],int(x['particle']),
                                    color='magenta'),
                  axis=1)
    ax.grid(False)
    return ax

def plot_moving_particles(t_cor,img_bright=None,frame=None,framei=0,particle2color=None,test=False,outd=None):
    plotp=f'{outd}/{framei:03d}.jpeg'
    plt.figure(figsize=[5,5])
    ax=plt.subplot(111)
    if not ((img_bright is None) and (frame is None)):
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
        
def make_gif(cellcfg=None,frames=None,t_cor=None,img_bright=None,
             outd=None,particle2color=None,
             test=False,force=False):
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
                from roux.lib.stat.norm import rescale
                from roux.lib.plot.colors import saturate_color
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
        _framei=framei
    if not test:    
        plt.close('all')
        com=f"convert -delay 10 -loop 0 {outd}/*.jpeg {gifp}"
        from roux.lib.io_sys import runbashcmd
        runbashcmd(com)
    return gifp

def plot_trajectories_stats(df,coly,colx='frame',rescalex=True,
                    label=None,
                    axvlinex=None,
                    params_plot={'color':'k','alpha':0.5},
                   fig=None,ax=None):
    fig=plt.figure() if fig is None else fig
    axin=not ax is None
    ax=plt.subplot() if ax is None else ax
    label=label if not hasattr(df,'label') else df.name       
    df=df.sort_values(colx)
    if rescalex:
        xmin=df[colx].min()
        df[colx]=range(len(df)) 
    ax=df.plot(x=colx,y=coly,ax=ax,label=label,**params_plot)
    ax.set_xlim(df[colx].min(),df[colx].max()+2)
    ax.set_ylim(df[coly].min(),df[coly].max())    
    if not axvlinex is None:
        ax.axvline(axvlinex-(xmin if rescalex else 0),
               label='inflection point',color='gray')
    ax.set_xlabel(colx);ax.set_ylabel(coly)
    return ax
                         