from rohan.global_imports import *

def make_project_cfg():
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
                                                              
def make_cell_cfg(cfg,cells,cellbox):
    outp=f"{cfg['trials'][trial]['datad']}/cells/cell{celli+1:08d}/"
    cfgp=f"{cellcfg['outp']}/cfg.yml"
    if not exists(cellcfg['cfgp']) or force:
        cellcfg={}
        cellcfg['outp']=outp
        cellcfg['cfgp']=cfgp
        cellcfg['test']=test
        cellcfg['force']=force
        cellcfg['cellbrightp']=f"{cellcfg['outp']}/cellbright.npy"
        cellbright=cells[cellbox[2]:cellbox[3],cellbox[0]:cellbox[1]]

        if not exists(dirname(cellcfg['cellbrightp'])): 
            makedirs(dirname(cellcfg['cellbrightp']),exist_ok=True)
        np.save(cellcfg['cellbrightp'], cellbright) 
        # only one cell per box
        cellcfg['cellbrightmaskp']=f"{cellcfg['outp']}/cellbrightmask.npy"
        cellbrightmask=filter_regions(cellbright.astype(int),prop_type='centroid_x',mn=70,mx=80)==0
        np.save(cellcfg['cellbrightmaskp'], cellbrightmask)

        cellframeps=[]
        cellframesmaskedps=[]
        for framei,frame in enumerate(frames):
            cellframe=frame[cellbox[2]:cellbox[3],cellbox[0]:cellbox[1]]
            cellframep=f"{cellcfg['outp']}/cellframe/frame{framei:08d}.npy"
            if not exists(dirname(cellframep)): 
                makedirs(dirname(cellframep),exist_ok=True)
            np.save(cellframep, cellframe)
            cellframeps.append(cellframep)

            cellframemasked=cellframe.copy()
            cellframemasked[cellbrightmask]=0
            cellframemaskedp=f"{cellcfg['outp']}/cellframesmasked/frame{framei:08d}.npy"
            if not exists(dirname(cellframemaskedp)): 
                makedirs(dirname(cellframemaskedp),exist_ok=True)
            np.save(cellframemaskedp, cellframemasked)
            cellframesmaskedps.append(cellframemaskedp)

        cellcfg['cellframeps']=cellframeps
        cellcfg['cellframesmaskedps']=cellframesmaskedps
        yaml.dump(cellcfg,open(cellcfg['cfgp'],'w'))
    else:
        cellcfg=read_dict(cfgp)
    return cellcfg
                                                              