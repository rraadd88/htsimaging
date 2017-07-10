def test_configure():
    from dms2dfe import configure    
    configure.main('prj')

    import htsimaging
    from imaging.lib.spt import expt_dh2expt_info,nd2msd
