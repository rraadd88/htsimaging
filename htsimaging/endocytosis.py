from skimage import io 
from glob import glob,iglob
from rohan.global_imports import *
import pims
import trackpy as tp
import logging
import yaml
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

def run_trials(prjd):
    cfg={'prjd':prjd}
    cfgp=f"{cfg['prjd']}/cfg.yml"
    cfg['cfgp']=cfgp
    if not exists(cfgp):
        cfg['trials']={basename(d):{'datad':d} for d in glob(f"{cfg['prjd']}/*") if not (d.endswith('.yml') or d.endswith('.yaml'))}
        for k in cfg['trials']:
            cfg['trials'][k]['gfp']=[p for p in glob(f"{cfg['trials'][k]['datad']}/*tif") if '_t' in p]
            cfg['trials'][k]['bright']=[p for p in glob(f"{cfg['trials'][k]['datad']}/*tif") if not '_t' in p]
            cfg['trials'][k]['plotd']=f"{cfg['trials'][k]['datad']}/plot"
                                                        
        # QC
        ## 1 CHECK BLEACHING
        for k in cfg['trials']:
            frames = pims.ImageSequence(np.sort(cfg['trials'][k]['gfp']), as_grey=True)
            dframes=pd.DataFrame({i:np.ravel(np.array(frame)) for i,frame in enumerate(frames)})
            plotp=f"{cfg['trials'][k]['plotd']}/plot_check_bleaching.png"
            makedirs(dirname(plotp),exist_ok=True)

            df=dframes.describe().T
            ax=plot_mean_std(df,cols=['mean','min','50%'])
            ax.set_xlabel('time points');ax.set_ylabel('intensity')
            plt.savefig(plotp)
        yaml.dump(cfg,open(cfgp,'w'))   
                                                        
    if not 'yeastspotter_srcd' in cfg:
        from htsimaging.lib.segment import run_yeastspotter
        cfg['yeastspotter_srcd']='yeast_segmentation'
        cfg=run_yeastspotter(cfg,test=True)
        yaml.dump(cfg,open(cfgp,'w'))
    #     break
    if not '' in cfg:

        yaml.dump(cfg,open(cfgp,'w'))

# assembling:
parser = argh.ArghParser()
parser.add_commands([run_trials])

if __name__ == '__main__':
    parser.dispatch()