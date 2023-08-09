PRO FMwLASCO20130424_13
@gcs_config
;To load images from CORs A&B plus LASCO and call rtsccguicloud with forward model from Thernisien et al. (2009)
SAVE_FILE=GCSD_PATH+'/GCS_20130424(1)/3.sav'
secchipath=DATA_PATH+'/stereo/secchi/L0'
lascopath=DATA_PATH+'/soho/lasco/level_1/c2'
;opens images of interest and backgrounds
;SECCHI A
ima=sccreadfits(secchipath+'/a/seq/cor1/20130424(1)/preped/20130424_064000_0B4c1A.fts', hdreventa)
imaprev=sccreadfits(secchipath+'/a/seq/cor1/20130424(1)/preped/20130424_063000_0B4c1A.fts',hdreventpa)
;SECCHI B
imb=sccreadfits(secchipath+'/b/seq/cor1/20130424(1)/preped/20130424_064000_0B4c1B.fts', hdreventb)
imbprev=sccreadfits(secchipath+'/b/seq/cor1/20130424(1)/preped/20130424_063000_0B4c1B.fts', hdreventpb)
;gets masks and scales
ma=get_smask(hdreventa)
mb=get_smask(hdreventb)
;a=bytscl(rebin(alog10(ma*(ima-imaprev) > 1e-12 < 1e-10),512,512))
;b=bytscl(rebin(alog10(mb*(imb-imbprev) > 1e-12 < 1e-10),512,512))
a=sigrange(rebin((ma*(ima)-ma*(imaprev)),512,512))   ;better for COR1
b=sigrange(rebin((mb*(imb)-mb*(imbprev)),512,512)) ;better for COR1
;a=sigrange(rebin(ma*sigrange(ima)-ma*sigrange(imaprev),512,512)) ;better for COR2
;b=sigrange(rebin(mb*sigrange(imb)-mb*sigrange(imbprev),512,512)) ;better for COR2
;opens LASCO of interest and background
lasco1=readfits(lascopath+'/20130424(1)/25456655.fts', lasco1hdr);+'/C3/L0/20130328/32334111.fts', lasco1hdr)
lasco0=readfits(lascopath+'/20130424(1)/25456654.fts', lasco0hdr);+'/C3/L0/20130328/32334110.fts', lasco0hdr)
lascohdr = LASCO_FITSHDR2STRUCT(lasco1hdr)
print, 'DATE_OBS:  ', lascohdr.date_obs
lasco1=rebin(lasco1, 512, 512)
lasco0=rebin(lasco0, 512, 512)
;lasco=rebin(((lasco1-lasco0)>(-2.e-11)<4.0e-11),512,512)   ;para C3:  >(-3.e-12)<3.e-12), 512,512)
lasco=sigrange(lasco1-lasco0)
;lasco= lasco>(-0.83e-9) <7.5e-10
;EUVIs
imeuvib=sccreadfits(secchipath+'/b/img/euvi/20130424(1)/preped/20130424_054530_14euB.fts', heuvib)
imeuvia=sccreadfits(secchipath+'/a/img/euvi/20130424(1)/preped/20130424_054530_14euA.fts', heuvia)
;This one replaces one of the EUVIs for SDO AIA
;imeuvib=sccreadfits(DATA_PATH+'/sdo/aia/L1/193/20130424(1)/preped/AIA20130424_054506_0193.fits', heuvib)
;secchi_prep,eveuvia,heuvia,imeuvia,/PRECOMMCORRECT_ON
;secchi_prep,eveuvib,heuvib,imeuvib,/PRECOMMCORRECT_ON
imea=alog10(rebin(imeuvia,512,512) > 1)
imeb=alog10(rebin(imeuvib,512,512) > 1)
;calls application
;calls application, loads and saves the fit parameters in SAVE_FLIE
if file_test(SAVE_FILE) eq 1 then restore, SAVE_FILE
rtsccguicloud, a, b, hdreventa, hdreventb, imlasco=lasco, hdrlasco=lascohdr, imeuvia=imea, hdreuvia=heuvia, imeuvib=imeb, hdreuvib=heuvib, sgui=sgui, sparaminit=sgui
save, filename=SAVE_FILE, sgui
END

