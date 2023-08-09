;Loads input data to run rtsccguicloud with the forward model from Thernisien et al. (2009)
pro FMwLASCO20130531_10

@gcs_config
;Differential image scaling [SECCHIA,SECCHIA, SECCHIB,SECCHIB, AIA, AIA]
RNG=0.5*[-1,1,-1,1,-1,1]

;Global Paths
SAVE_FILE=DATA_PATH+'/Polar_Observations/Polar_Documents/GCSs/GCS_20130531/0.sav'
secchipath=DATA_PATH+'/stereo/secchi/L0'
lascopath=DATA_PATH+'/soho/lasco/level_1/c2'
sdopath=DATA_PATH+'/sdo'

;Reads and prepeares data

;date for both stereos
SDATE0='20130531'
SDATE ='20130531'
STIME0='051530'
STIME ='053030'

;SECCHIA/EUVI (to be used as CORA)
peuvia0=secchipath+'/a/img/euvi/'+SDATE0+'/'+SDATE0+'_'+STIME0+'_n4euA.fts'
peuvia =secchipath+'/a/img/euvi/'+SDATE+'/'+SDATE+'_'+STIME+'_n4euA.fts'
euvia0=sccreadfits(peuvia0,euviahdr0)
euvia =sccreadfits(peuvia,euviahdr)
secchi_prep,peuvia0,euviahdr0,euvia0,/PRECOMMCORRECT_ON ; preprocessing 
secchi_prep,peuvia,euviahdr,euvia,/PRECOMMCORRECT_ON ; preprocessing 
euvia=sigrange(rebin(euvia-euvia0,512,512))

;SECCHIB/EUVI (to be used as CORB)
peuvib0=secchipath+'/b/img/euvi/'+SDATE0+'/'+SDATE0+'_'+STIME0+'_n4euB.fts'
peuvib =secchipath+'/b/img/euvi/'+SDATE+'/'+SDATE+'_'+STIME+'_n4euB.fts'

euvib0=sccreadfits(peuvib0,euvibhdr0)
euvib =sccreadfits(peuvib,euvibhdr)
secchi_prep,peuvib0,euvibhdr0,euvib0,/PRECOMMCORRECT_ON; preprocessing 
secchi_prep,peuvib,euvibhdr,euvib,/PRECOMMCORRECT_ON ; preprocessing 
euvib=sigrange(rebin(euvib-euvib0,512,512))

;SDO/AIA (to be used as LASCO)
DATE0=['2013','05','31']
TIME0=['05','15','06']
DATE=['2013','05','31']
TIME=['05','30','06']

; AIA preprocessing if necesary
oaia0=sdopath+'/aia/L1/193/'+strjoin(DATE0)+'/preped/AIA'+strjoin(DATE0)+'_'+strjoin(TIME0)+'_0193.fits'
oaia =sdopath+'/aia/L1/193/'+strjoin(DATE)+'/preped/AIA'+strjoin(DATE)+'_'+strjoin(TIME)+'_0193.fits'
if ~file_test(oaia0) then begin 
  paia0=sdopath+'/aia/L1/193/'+strjoin(DATE0)+'/aia.lev1.193A_'+DATE0[0]+'-'+DATE0[1]+'-'+$ 
   DATE0[2]+'T'+TIME0[0]+'_'+TIME0[1]+'_'+TIME0[2]+'.84Z.image_lev1.fits'
  aia_prep, paia0, -1, outdir=file_dirname(oaia0),/do_write_fits ; preprocessing
endif
if ~file_test(oaia) then begin
  paia =sdopath+'/aia/L1/193/'+strjoin(DATE)+'/aia.lev1.193A_'+DATE[0]+'-'+DATE[1]+'-'+$ 
   DATE[2]+'T'+TIME[0]+'_'+TIME[1]+'_'+TIME[2]+'.84Z.image_lev1.fits'
  aia_prep, paia, -1, outdir=file_dirname(oaia),/do_write_fits ; preprocessing
endif

;reads AIA preprocessed data
aia0=sccreadfits(oaia0, aiahdr)
aia =sccreadfits(oaia, aiahdr)
aia =sigrange(rebin(aia-aia0, 512, 512))

;smoothing
;euvib=smooth(euvib,3)

;diff images scaling
euvia=bytscl(euvia,min=mean(euvia)+RNG[0]*stddev(euvia), max=mean(euvia)+RNG[1]*stddev(euvia))
euvib=bytscl(euvib,min=mean(euvib)+RNG[2]*stddev(euvib), max=mean(euvib)+RNG[3]*stddev(euvib))
aia=bytscl(aia,min=mean(aia)+RNG[4]*stddev(aia), max=mean(aia)+RNG[5]*stddev(aia))


;calls application, loads and saves the fit parameters in SAVE_FLIE
Result = file_test(SAVE_FILE)
if Result eq 1 then restore, SAVE_FILE
rtsccguicloud, euvia, euvib, euviahdr, euvibhdr, imlasco=aia, hdrlasco=aiahdr, sgui=sgui, sparaminit=sgui
save, filename=SAVE_FILE, sgui

end

