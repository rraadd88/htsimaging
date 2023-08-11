#!/usr/bin/env python

from os.path import splitext, join, exists, isdir,basename,abspath,dirname
from os import makedirs
import trackpy as tp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pims
# import pims_nd2
# import nd2reader
from glob import glob
import logging
logging.basicConfig(format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..): %(message)s',level=logging.DEBUG) # 

from os.path import dirname, basename
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def plot_msd(
    imsd: pd.DataFrame,
    emsd: pd.DataFrame,
    scale: str="log",
    plot_fh: str=None,
    time_points_max: int = None,
    ax: plt.Axes=None,
    ) -> plt.Axes:
    """Plot MSD.

    Args:
       imsd (pd.DataFrame): mean squared displacement of each particle.
       emsd (pd.DataFrame): ensemble mean squared displacement of particles.
       scale (str, optional): axis scale. Defaults to "log".
       plot_fh (str, optional): output path. Defaults to None.
       params_msd (dict, optional): parameters of MSD. Defaults to { "mpp":0.0645, "fps":0.2, "max_lagtime":100 }.
       ax (plt.Axes, optional): subplot. Defaults to None.

    Returns:
       plt.Axes: subplot
       
    Examples:
        Calculate `time_points_max`:
        ```
            params_msd: dict={
                "mpp":0.0645,
                "fps":0.2,
                "max_lagtime":100
                }
            time_points_max=params_msd["fps"]*params_msd["max_lagtime"]
        ```
    """
    if ax is None:
        ax=plt.gca()
    if not time_points_max is None:
        imsd=imsd.head(time_points_max)
        emsd=emsd.head(time_points_max)

    imsd.plot(
        # x=imsd.index,
        legend=None,alpha=0.75,ax=ax)
    ax.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',
           xlabel='lag time $t$')
    if scale=="log":
        ax.set_xscale('log')
        ax.set_yscale('log')

    emsd.plot(
        # x=emsd.index,
        style='o',legend=None,ax=ax)
    if scale=="log":
        ax.set_xscale('log')
        ax.set_yscale('log')
    ax.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',
           xlabel='lag time $t$')
    plt.tight_layout()
    if not plot_fh is None: 
        ax.figure.savefig(plot_fh);
    return ax

def plot_emsd(
    expt_data: pd.DataFrame,
    color: str='k',
    scale: str="log",
    plot_fh: str=None,
    ax: plt.Axes=None,
    ) -> plt.Axes:
    """
    Plot ensemble mean squared displacement of particles.

    Args:
       expt_data (pd.DataFrame): input data.
       color (str, optional): color. Defaults to 'k'.
       scale (str, optional): scale of the axes. Defaults to "log".
       plot_fh (str, optional): output path. Defaults to None.
       ax (plt.Axes, optional): subplot. Defaults to None.

    Returns:
       plt.Axes: subplot
    """
    if ax is None:
        plt.figure(figsize=(6, 3))
        ax=plt.subplot(121)
    expt_data.plot(x=expt_data.index,style='o',c=color,ax=ax,markeredgecolor = 'none')
    if scale=="log":
        ax.set_xscale('log')
        ax.set_yscale('log')
    ax.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',
                xlabel='lag time $t$')
    ax.legend(loc = 'center left', bbox_to_anchor = (1.0, 0.5)) #bbox_to_anchor=(2.2, 1.0)
    plt.tight_layout()
    if not plot_fh is None: 
        ax.figure.savefig(plot_fh);
#         plt.clf();plt.close()
    return ax 

## kinetics 
def plot_kin(
    traj: pd.DataFrame,
    dparams: pd.DataFrame=None,
    ctime: str='time (s)',
    fit_eqn: str=None,
    ylabel: str='',
    label: str='',
    color: str='b',
    ax1: plt.Axes=None,
    plot_fh: str=None,
    ) -> plt.Axes:
    """Plot kinetics.

    Args:
        traj (pd.DataFrame): table containing the trajectories.
        dparams (pd.DataFrame, optional): table containing the parameters. Defaults to None.
        ctime (str, optional): time point. Defaults to 'time (s)'.
        fit_eqn (str, optional): fit equation name. Defaults to None.

    Returns:
        plt.Axes
    """
    traj_mean=traj.set_index(ctime)
    traj_mean=traj.T.mean()
    traj_std=traj.T.std()
    if fit_eqn=="logistic4":
        params=["Minimum asymptote","Hill's slope","Inflection point","Maximum asymptote"]
    elif fit_eqn=="logistic5":
        params=["Minimum asymptote","Hill's slope","Inflection point","Maximum asymptote","Asymmetry factor"]
    elif fit_eqn=='power':
        params=["power law exponent","amplitude"]
        
    if ax1 is None:
        plt.figure(figsize=(8,3))
        ax1=plt.subplot(131)
    ax1.set_xlabel(ctime)
    ax1.set_ylabel(ylabel)
    ax1.plot(traj_mean.index, traj_mean,'o',zorder=2, 
             markerfacecolor='none',markeredgecolor=color,alpha=0.5,
#             label='sd',
            )
    ax1.legend(label,loc='center', 
               bbox_to_anchor=(1.63, 0.5))
    ax1.fill_between(traj_std.index, \
                     traj_mean-traj_std, \
                     traj_mean+traj_std, \
                     color=color, alpha=0.2,zorder=0)
    
    if not dparams is None:
        params_mean=dparams.mean()
        params_std=dparams.std()
        p0=params_mean
    else:
        p0=np.zeros(len(params))
    if fit_eqn=="logistic4":
        plsq = leastsq(logistic4_residuals, p0, args=(traj_mean, traj_mean.index))
        ax1.plot(traj_mean.index.tolist(),logistic4_peval(traj_mean.index,plsq[0]),label=label,color=color,zorder=2)        
    elif fit_eqn=="logistic5":
        plsq = leastsq(logistic5_residuals, p0, args=(traj_mean, traj_mean.index))
        ax1.plot(traj_mean.index.tolist(),logistic5_peval(traj_mean.index,plsq[0]),label=label,color=color,zorder=2)        
    elif fit_eqn=="power":
        plsq = leastsq(power_residuals, p0, args=(traj_mean, traj_mean.index))
        ax1.plot(traj_mean.index.tolist(),power_peval(traj_mean.index,plsq[0]),label=label,color=color,zorder=2,lw=2)        
                
    ax1.set_ylim([0,ax1.get_ylim()[1]])
    plt.tight_layout()
    if plot_fh!=None:
        plt.savefig(plot_fh, format='pdf')
    return ax1

def plot_kin_all(
    expt_dh: str,
    imsd_fhs: list,
    ):
    """
    Plot multiple kinetics plots.
    """
    fig=plt.figure(figsize=(8,3))
    colors=cm.rainbow(np.linspace(0,1,len(imsd_fhs)))
    ylimmx=[]
    for i,imsd_fh in enumerate(imsd_fhs):
        imsd_flt_fh='%s.imsd_flt' % imsd_fh
        imsd_flt=pd.read_csv(imsd_flt_fh).set_index('lag time [s]')
        params_flt_fh='%s.params_flt' % imsd_fh
        params_flt=pd.read_csv(params_flt_fh)
        ax1=plt.subplot(131)
        ax1=plot_kin(imsd_flt,
                 params_flt.loc[:,['power law exponent','power law constant']],
                 ctime='lag $t$',
                 fit_eqn='power',
                 smp_ctrl=None,
                 ylabel='MSD ($\mu m^{2}$)',
                 color=colors[i],
                 ymargin=0.2,
                 ax1=ax1,
                 label=basename(imsd_flt_fh).split('_replicate', 1)[0].replace('_',' '),
                 plot_fh=None)
        ylimmx.append(np.max(imsd_flt.max()))
    ax1.set_ylim([0,np.max(ylimmx)])
    ax1.legend(bbox_to_anchor=[1,1,0,0],loc=2)
    fig.savefig('%s/plot_flt_%s.pdf' % (expt_dh,dirname(expt_dh)))

    
def plot_trajectories_stats(
    df: pd.DataFrame,
    coly: str,
    colx: str='frame',
    rescalex: bool=True,
    label: str=None,
    axvlinex=None,
    params_plot: dict={'color':'k','alpha':0.5},
    fig=None,
    ax: plt.Axes=None,
    ) -> plt.Axes:
    """
    Plot statistics of the particle trajectories.
    """
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

## on a cell image
def plot_trajectories(
    traj,
    image,
    colorby: str='particle',
    mpp: float=None,
    label: str=False,
    cmap: str=None,
    t_column: str=None,
    pos_columns: list=None,
    plot_style: dict={},
    params_text: dict={'ha':'center','va':'center'},
    ax: plt.Axes=None,
    **kwargs,
    ) -> plt.Axes:
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
    # Background image
    ax=plt.subplot() if ax is None else ax
    ax.imshow(image,cmap=get_cmap_subset('binary',vmin=0.2,vmax=0))
    # Trajectories
    ax.set_xlim(0,image.shape[0])
    ax.set_ylim(0,image.shape[1])
    ax.set_ylim(ax.get_ylim()[::-1])
    return ax