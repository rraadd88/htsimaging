def test_configure():
    from dms2dfe import configure    
    configure.main('prj')

    from htsimaging import trackinginfo2emsd
    trackinginfo2emsd.main(expt_dh)