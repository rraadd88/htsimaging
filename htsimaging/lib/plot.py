from rohan.global_imports import *

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
    from rohan.dandage.plot.colors import get_cmap_subset
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

def plot_properties_cell(cellcfg,df2,cols_colorby):
    from rohan.dandage.stat.norm import rescale
    from rohan.dandage.plot.contour import plot_contourf
    ncols=4
    nrows=int(np.ceil(len(cols_colorby)/4))
    fig,axes=plt.subplots(nrows,ncols,
#                           sharex=True,sharey=True,
                          figsize=[nrows*5,ncols*3])
    metric_type='max'
    for axi,(colorby,ax) in enumerate(zip(cols_colorby,np.ravel(axes))):
        fig,ax=plot_contourf(df2['x median'],df2['y median'],
                      rescale(df2[f'{colorby} {metric_type}']),
                     ax=ax,
                     fig=fig,
                      cbar=True if ((axi+1) % 4)==0 else False,
                     params_contourf={'cmap':'binary','vmin':0,'vmax':1},)
        ax=df2.plot.scatter(x='x median',y='y median',color='lime',
                         s=10,
                         alpha=0.2,
                         ax=ax)   
        ax.contour(np.load(cellcfg['cellbrightp']), [0.5], linewidths=1, linestyles='dashed',colors='cyan')    
        ax.set_title(colorby)
        ax.invert_yaxis()
        ax.axis('equal')
#         ax.set_ylim(ax.get_ylim()[::-1])
#         ax.set_axis_off()

def dist_signal(img,threshold=None,label_threshold=None,params_hist={},ax=None):
    ax=plt.subplot() if ax is None else ax 
    a=np.ravel(img)
    _=ax.hist(np.extract(a<np.quantile(a,0.9),a),**params_hist)
    # ax.set_xlim(np.min(a),np.quantile(a,0.9))
    if not threshold is None:
        ax.axvline(threshold,color='r',label=label_threshold,linestyle='dashed')
        ax.legend()
    ax.set_xlabel('signal')
    ax.set_ylabel('density')
    return ax

def make_gif(frames,t_cor,outd=None,test=False,force=False):
    from rohan.dandage.plot.colors import get_cmap_subset
    if not outd is None:
        makedirs(outd,exist_ok=True)
        gifp=f"{dirname(outd)}/vid.gif"
    else:
        test=True
        gifp=''
    if not exists(gifp) or force:
        for framei,frame in enumerate(frames):
            plotp=f'{outd}/{framei:02d}.png'
            plt.figure(figsize=[10,10])
            ax=plt.subplot(111)
            ax.imshow(frame,cmap=get_cmap_subset('binary_r',vmin=0.2,vmax=1),alpha=0.8,
                     )
            ax.text(frame.shape[0],frame.shape[1],f"frame={framei:05d}",ha='right',va='bottom',size='20',color='y')
            if not framei==0:
                # traj of particle
                for p in t_cor.loc[(t_cor['frame']==framei),'particle'].unique():
                    # traj of particle
                    framei_min=t_cor['frame min'].unique()[0]
                    df=t_cor.loc[((t_cor['particle']==p) \
                                  & (t_cor['frame'].isin(range(framei_min,framei+1)))\
                                  & (t_cor['x'].between(0,frame.shape[0]))\
                                  & (t_cor['y'].between(0,frame.shape[1]))),:].sort_values(by=['frame','particle'])
                    if test:
                        print(df['frame'])
                    if len(df)!=0:
                        if not 'move' in df:
                            df.plot.line(x='x',y='y',lw=7,
                                c='limegreen',
                                legend=False,ax=ax)                        
                        else:
                            df.plot.line(x='x',y='y',lw=7,
                                c='limegreen' if df.loc[:,['particle','move']].drop_duplicates().set_index('particle').loc[p,'move']==1 else 'magenta',
                                legend=False,ax=ax)
                ax.set_xlim(0,frame.shape[1])
                ax.set_ylim(0,frame.shape[1])
                plt.axis('off')
                if test:
                    if framei==3:
                        return ax
                if not test:    
                    makedirs(dirname(plotp),exist_ok=True)
                    plt.savefig(plotp)
                plt.clf();plt.close();
            _framei=framei
        if not test:    
            plt.close('all')
            com=f"convert -delay 10 -loop 0 {outd}/*.png {gifp}"
            runbashcmd(com)
    return gifp
