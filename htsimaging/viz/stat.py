#!/usr/bin/env python
"""Visualization of the statistics."""
import numpy as np
import matplotlib.pyplot as plt

def dist_signal(
    img,
    threshold: float=None,
    label_threshold: float=None,
    params_axvline: dict={'color':'r','linestyle':'dashed'},
    ax: plt.Axes=None,
    **kws,
    ) -> plt.Axes:
    """Plot the distribution of intensity.

    Args:
        img (_type_): inpput image
        threshold (float, optional): threshold applied. Defaults to None.
        label_threshold (float, optional): label of the threshold. Defaults to None.
        params_axvline (_type_, optional): parameters provided to the vertical line plot. Defaults to {'color':'r','linestyle':'dashed'}.
        ax (plt.Axes, optional): subplot object. Defaults to None.

    Keyword Args:
        parameters provided to the `hist` function. 
    
    Returns:
        plt.Axes
    """
    ax=plt.subplot() if ax is None else ax 
    a=np.ravel(img)
    a = a[~np.isnan(a)]
    a=np.extract(a<np.quantile(a,0.9),a)
    _=ax.hist(a,linewidth=0,**kws,)
    # ax.set_xlim(np.min(a),np.quantile(a,0.9))
    if not threshold is None:
        ax.axvline(threshold,label=label_threshold,**params_axvline)
        ax.legend()
    ax.set_xlabel('signal')
    ax.set_ylabel('density')
    return ax

def plot_summary_stats(
    input_paths: list,
    ax: plt.Axes=None,
    ) -> plt.Axes:
    """Plot summary stats for a set of images e.g. time-lapse images.

    Args:
        input_paths (list): list of paths of the images.
        ax (plt.Axes, optional): subplot object. Defaults to None.

    Returns:
        plt.Axes
    """
    import pandas as pd
    from htsimaging.lib.io import read_image
    
    if ax is None:
        fig,ax=plt.subplots(figsize=[3,3])
    frames = [read_image(p) for p in input_paths]
    dframes=pd.DataFrame({i:np.ravel(np.array(frame)) for i,frame in enumerate(frames)})
    df=dframes.describe().T
    ax=df.loc[:,['mean','min','50%']].plot(ax=ax)
    return ax
