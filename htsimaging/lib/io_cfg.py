from roux.global_imports import *
from os.path import isdir
import pims
import sys

def make_project_cfg(prjd,bright_fn_marker,cores=1,
                    test=False,force=False):
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
        from roux.lib.plot.line import plot_mean_std                                                        
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

def make_cell_cfg(cfg,frames,cells,trial,celli,cellbox,
                 params_get_signal_summary_by_roi={'xy_center':None,'width':20,
                                  'fun_summary_frame':'min',
                                  'fun_summary_frames':'median'
                                 },
                  test=False,force=False,
                        ):
                                                              
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
        cellframe_type2frames={'cellframes':{'frames':[],'ps':[]},
                               'cellframes_masked':{'frames':[],'ps':[]},
                               'cellframes_masked_substracted':{'frames':[],'ps':[]},
                              }
        for dn in sorted(list(cellframe_type2frames.keys())):
            if not exists(f"{cellcfg['outp']}/{dn}"): 
                makedirs(f"{cellcfg['outp']}/{dn}",exist_ok=True)
        from htsimaging.lib.utils import get_signal_summary_by_roi
        for framei,frame in enumerate(frames):
            cellframe=frame[cellbox[2]:cellbox[3],cellbox[0]:cellbox[1]]
            
            cellframe_masked=cellframe.copy()
            cellframe_masked[cellbrightmask]=0
                         
            cellframe_masked_substracted=cellframe_masked.copy()
            signal_cytoplasm=get_signal_summary_by_roi([cellframe_masked_substracted],**params_get_signal_summary_by_roi)
            cellframe_masked_substracted=np.where(cellframe_masked_substracted<signal_cytoplasm*1.5,
                                                              signal_cytoplasm,
                                                              cellframe_masked_substracted) 
            cellframe_masked_substracted[cellbrightmask]=0

            cellframe_type2frames['cellframes']['frames'].append(cellframe)
            cellframe_type2frames['cellframes_masked']['frames'].append(cellframe_masked)
            cellframe_type2frames['cellframes_masked_substracted']['frames'].append(cellframe_masked_substracted)
            for dn in sorted(list(cellframe_type2frames.keys())):
                cellframe_type2frames[dn]['ps'].append(f"{cellcfg['outp']}/{dn}/frame{framei:08d}.npy")
                np.save(cellframe_type2frames[dn]['ps'][framei],
                        cellframe_type2frames[dn]['frames'][framei])

        for dn in sorted(list(cellframe_type2frames.keys())):
            cellcfg[dn]=cellframe_type2frames[dn]['ps']
        #gfp min max
        for funn in ['min','max']:
            cellcfg[f'cellgfp{funn}p']=f"{cellcfg['outp']}/cellgfp{funn}.npy"
            frame=getattr(np,f"a{funn}")(cellframe_type2frames['cellframes_masked_substracted']['frames'],axis=0)
            np.save(cellcfg[f'cellgfp{funn}p'], frame)
        del frame

        cellcfg['signal_cytoplasm']=get_signal_summary_by_roi(cellframe_type2frames['cellframes_masked']['frames'],
                                                              **params_get_signal_summary_by_roi)

        yaml.dump(cellcfg,open(cellcfg['cfgp'],'w'))
    else:
        cellcfg=read_dict(cfgp)
    return cellcfg
                                                              
