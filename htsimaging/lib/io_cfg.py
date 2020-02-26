from rohan.global_imports import *
from os.path import isdir
import pims

def make_project_cfg(prjd,bright_fn_marker,test,force,cores):
    prjd=abspath(prjd)
    cfgp=f"{prjd}/cfg.yml"
    if not exists(cfgp) or force:    
        print('making cfg')
        cfg={'prjd':prjd}
        cfg['cfgp']=cfgp
        cfg['cores']=cores
        cfg['bright_fn_marker']=bright_fn_marker
        if 'T1C1' in cfg['bright_fn_marker']:
            cfg['magnification']=150
        elif cfg['bright_fn_marker']=='_t':
            cfg['magnification']=100
        else:
            logging.error("unknown bright_fn_marker {cfg['bright_fn_marker']}")
            sys.exit()
        cfg['trials']={basename(d):{'datad':d} for d in glob(f"{cfg['prjd']}/*") if (isdir(d) and basename(d).replace('/','')!='segmentation_cell' and not basename(d).startswith('_'))}
        trials_bad=[]
        for k in cfg['trials']:
            cfg['trials'][k]['bright']=[abspath(p) for p in glob(f"{cfg['trials'][k]['datad']}/*tif") if cfg['bright_fn_marker'] in p and not p.endswith('.segmented.tif')]
            cfg['trials'][k]['gfp']=[abspath(p) for p in glob(f"{cfg['trials'][k]['datad']}/*tif") if not (cfg['bright_fn_marker'] in p or 'segmented' in p) ]
            cfg['trials'][k]['plotd']=f"{cfg['trials'][k]['datad']}/plot"
            makedirs(cfg['trials'][k]['plotd'],exist_ok=True)
            if len(cfg['trials'][k]['bright'])==0 or len(cfg['trials'][k]['gfp'])==0:
                trials_bad.append(k)
        for trial in trials_bad:
            del cfg['trials'][trial]
        yaml.dump(cfg,open(cfgp,'w'))   
        # QC
        ## 1 CHECK BLEACHING
        from rohan.dandage.plot.line import plot_mean_std                                                        
        for k in cfg['trials']:
            frames = pims.ImageSequence(np.sort(cfg['trials'][k]['gfp']), as_grey=True)
            dframes=pd.DataFrame({i:np.ravel(np.array(frame)) for i,frame in enumerate(frames)})
            plotp=f"{cfg['trials'][k]['plotd']}/plot_check_bleaching.png"
            makedirs(dirname(plotp),exist_ok=True)

            df=dframes.describe().T
            ax=plot_mean_std(df,cols=['mean','min','50%'])
            ax.set_xlabel('time points');ax.set_ylabel('intensity')
            savefig(plotp)
    #         yaml.dump(cfg,open(cfgp,'w'))
    else:
        cfg=yaml.load(open(cfgp,'r'))
    return cfg

def make_cell_cfg(cfg,frames,cells,trial,celli,cellbox,test,force):
    outp=f"{cfg['trials'][trial]['datad']}/cells/cell{celli+1:08d}/"
    cfgp=f"{outp}/cfg.yml"
    if not exists(cfgp) or force:
        cellcfg={}
        cellcfg['outp']=outp
        cellcfg['plotp']=f"{outp}/plot"
        cellcfg['cfgp']=cfgp
        cellcfg['test']=test
        cellcfg['force']=force
        cellcfg['cellbrightp']=f"{cellcfg['outp']}/cellbright.npy"       
        cellbright=cells[cellbox[2]:cellbox[3],cellbox[0]:cellbox[1]]
                                                              
        if not exists(dirname(cellcfg['cellbrightp'])): 
            makedirs(dirname(cellcfg['cellbrightp']),exist_ok=True)
        np.save(cellcfg['cellbrightp'], cellbright) 
        makedirs(dirname(cellcfg['outp']),exist_ok=True)
        makedirs(dirname(cellcfg['plotp']),exist_ok=True)
                                                              
        # only one cell per box
        cellcfg['cellbrightmaskp']=f"{cellcfg['outp']}/cellbrightmask.npy"
        from htsimaging.lib.utils import filter_regions
        cellbrightmask=filter_regions(cellbright.astype(int),prop_type='centroid_x',mn=70,mx=80)==0
        np.save(cellcfg['cellbrightmaskp'], cellbrightmask)
        cellframeps=[]
        cellframesmaskedps=[]
#         cellframes=[]
        cellframesmaskeds=[]
        for framei,frame in enumerate(frames):
            cellframe=frame[cellbox[2]:cellbox[3],cellbox[0]:cellbox[1]]
            cellframep=f"{cellcfg['outp']}/cellframe/frame{framei:08d}.npy"
            if not exists(dirname(cellframep)): 
                makedirs(dirname(cellframep),exist_ok=True)
            np.save(cellframep, cellframe)
            cellframeps.append(cellframep);#cellframes.append(cellframe)

            cellframemasked=cellframe.copy()
            cellframemasked[cellbrightmask]=0
            cellframemaskedp=f"{cellcfg['outp']}/cellframesmasked/frame{framei:08d}.npy"
            if not exists(dirname(cellframemaskedp)): 
                makedirs(dirname(cellframemaskedp),exist_ok=True)
            np.save(cellframemaskedp, cellframemasked)
            cellframesmaskedps.append(cellframemaskedp);cellframesmaskeds.append(cellframemasked)

        cellcfg['cellframeps']=cellframeps
        cellcfg['cellframesmaskedps']=cellframesmaskedps
        #gfp min max
        cellcfg['cellgfpminp']=f"{cellcfg['outp']}/cellgfpmin.npy"
        cellcfg['cellgfpmaxp']=f"{cellcfg['outp']}/cellgfpmax.npy"
        cellgfpmin=np.amax(cellframesmaskeds,axis=0)
        cellgfpmax=np.amin(cellframesmaskeds,axis=0)
        np.save(cellcfg['cellgfpminp'], cellgfpmin)
        np.save(cellcfg['cellgfpmaxp'], cellgfpmax)

        from htsimaging.lib.utils import get_signal_summary_by_roi
        cellcfg['signal_cytoplasm']=get_signal_summary_by_roi(cellframesmaskeds,
                                 xy_center=None,
                                width=20,
                                fun_summary_frame='min',
                                fun_summary_frames='median',)

        cellframesmaskedsubstractedps=[]
        for cellframesmasked in cellframesmaskeds:
            cellframesmaskedsubstractedp=f"{cellcfg['outp']}/cellframesmaskedsubstracted/frame{framei:08d}.npy"
            cellframesmaskedsubstracted=np.where(cellframesmasked<cellcfg['signal_cytoplasm']*1.5,
                                                              cellcfg['signal_cytoplasm'],
                                                              cellframesmasked)
            np.save(cellframesmaskedsubstractedp, cellframesmaskedsubstracted)
            cellframesmaskedsubstractedps.append(cellframesmaskedsubstractedp)
        cellcfg['cellframesmaskedsubstractedps']=cellframesmaskedsubstractedps
        yaml.dump(cellcfg,open(cellcfg['cfgp'],'w'))
    else:
        cellcfg=read_dict(cfgp)
    return cellcfg
                                                              
