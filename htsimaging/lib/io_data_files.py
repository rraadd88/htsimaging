




def expt_dh2expt_info(expt_dh):
    expt_info=pd.read_csv(expt_dh+"info")
    cols_del=["dh","fn_lead"]
    for col in expt_info.columns.tolist():
        if "Unnamed" in col:
            cols_del.append(col)
    for col in cols_del:
        if col in expt_info.columns.tolist():
            del expt_info[col]    
    expt_info2=expt_info.drop(["smpn"],axis=1).T
    expt_info2.columns=expt_info.loc[:,"smpn"]

    for col in expt_info2.columns.tolist():
        for i in range(len(expt_info2)):
            if not pd.isnull(expt_info2.loc["replicate %d" % (i+1),col]):
                expt_info2.loc[("replicate %d" % (i+1)),col]=['replicate %d' % (i+1),expt_info2.loc[("replicate %d" % (i+1)),col]]
            else:
                expt_info2.loc[("replicate %d" % (i+1)),col]=[np.nan,np.nan]
    expt_info2=expt_info2.reset_index()
    del expt_info2["index"]
    return expt_info2

def createinfo(expt_dh):
    info=pd.read_csv(expt_dh+"info")
    for i in range(len(info)):
        reps=glob("%s/%s*" % (info.loc[i,"dh"],info.loc[i,"fn_lead"]))
        for repi,rep in enumerate(reps):
            info.loc[i,"replicate %d" % (repi+1)]=rep
    #     break
    info.to_csv(expt_dh+"info")

