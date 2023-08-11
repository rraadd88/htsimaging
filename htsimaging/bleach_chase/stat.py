#!/usr/bin/env python
"""Statistical analysis of the bleach-chase"""
import pandas as pd

def exp1(x,lag,amplitude,rate):
    """
    One-term exponential.
    
    Parameters:
        x (list): input vector
        lag (float): lag   
        amplitude (float): amplitude
        rate (float): rate
    """
    return lag+(amplitude*(1-np.exp(-rate*x)))
def get_scores(
    df: pd.DataFrame,
    ) -> pd.DataFrame:
    """
    Calculate the rates other parameters.
    
    Parameters:
        df (pd.DataFrame): input table.
    
    Returns:
        pd.DataFrame: table with the parameters.
    """
    popt, pcov = curve_fit(
        exp1, 
        xdata=df['time'], 
        ydata=df['abundance rescaled'],
        p0=[0,1,0.01],
        bounds=([0,0.9,0],[0.1,1.1,1]),
    )

    perr = np.sqrt(np.diag(pcov))
    df_=pd.DataFrame({'score':popt,'score SD':perr,})
    df_.index=['lag','amplitude','rate']
    df_.index.name='variable'
    return df_
