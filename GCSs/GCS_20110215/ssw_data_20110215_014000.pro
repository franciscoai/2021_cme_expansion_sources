;+
;*************************************************************************
; Francisco A. Iglesias - UTN-FRM/GEAA - franciscoaiglesias@frm.utn.edu.ar
;
; Solar data download and preprocessing script.
; Uses SSW tootls (mainly from vso and jsoc) to download and preprocess data for various missions/instruments
;
; OUTPUT:
; 	-The output file structure is the same as in the oficial repository
;   -Existing ouput files ar NOT download again unless specified

; History:
;	20180409 - First version
;	20180510 - Changed the input format for AIA, added checking for file existence
; 20180511 - Added AIA preprocessing
; 20180514 - Added movie making and playing for AIA
; 20180524 - Added HMI data download
; 20180528 - Added SHARPS data download
; 20180808 - Added EUVIA data downloading, preprocessing and movie making
; 20181002 - CHANGED vso_get for simple spawns because the SSL is not working anymor in IDL 7.1
; 20181002 - Added do_euvia/b cases
; 20181002 - Added do_cora/b cases
;
; TODO:
;     - Change all cases to use the file ID returned by vso_search to replicate the repo structure.
;     - Add movie makeing to all cases
;*************************************************************************
;-
pro ssw_data_20110215_014000
time=SYSTIME(/SECONDS) ; just to check the script running time

; General Constants
 ODIR0='/media/hebe/ADATA838/data/francisco' ; Main output data directory

;Tasks
;For each Instrument, set to (not all opt are avilable for all inst.):
;[download, preprocess, make movie, play movie]
;  * 1 to do it
;  * 0 to do nothing
;  * 2 to overwrite ouput (when relevant)
;  * 3 to save individual frames instead of makeing a movie (only for make movie)
;SOHO
do_lasco=1*[1,1]  
;SDO
do_aia=1*[1,1,0,0]   
do_hmi=[0]
do_sharp=[0]   
;STEREO
do_euvia=0*[1,1,0,0]  
do_euvib=0*[1,1,0,0] 
do_cora=1*[1,1] 
do_corb=1*[1,1]
  
;**********Main

;***********************SOHO/LASCO using VSO
if total(do_lasco) gt 0 then begin
;download
 DATES=['2011/02/15','2011/02/15']
 TIMES=['02:15:00','03:45:00']
 TIMER=1*60.  ; Time resolution [s]
 ;WAVE='193'  ; Wavelenght [nm]
 INSTR='lasco' ; Instrument
 SOURCE='soho'; Observatory
 DETECTOR='c2'
 PHYSOBS='intensity'
 LEVEL='L1'  ; Data reduction level
 ;Only for preprocessing
 PDIR='level_1' ; Preprocessed data output dir is ODIR+'/'+'/'+PDIR
 FRMDIR='png_frames' ; diur to save the png frames 
 ;Only to make movie
 MNAME='movie.mvi'
 LOG=0 ; Applies ALOG10() function to image before byte scaling  
 FOV=[1536.+[0, 1024], 768.+[0, 1024]] ; size of the square FOV to use
 AVG=1 ; Num of frms to average 
 CT=[0,0,2000,0.6] ; Color table settings [table ID,range low(DN), range high (DN), gamma correction]

 ;play move
 RDIFF=0 ; set to play a runnig difference movie


;*****
 ; checks repository
 message, 'INFO: Checking repository.....' ,/informational 
 bus=vso_search(DATES[0]+' '+TIMES[0],DATES[1]+' '+TIMES[1], $ 
 source=SOURCE,instr=INSTR, detector=DETECTOR,physobs=PHYSOBS, sample=TIMER) ;Gets metadata from VSO
 fnamet=strjoin(strsplit(bus[0].time.start,':',/extract),'_') ; only the times
 for i=1, n_elements(bus.time.start)-1 do fnamet=[fnamet,strjoin(strsplit(bus[i].time.start,':',/extract),'_')]

 ; downloads files only if they do not exist in ODIR
 if do_lasco[0] gt 0 then begin
  message, 'INFO: Downloading files.....' ,/informational 
  for i=0, n_elements(bus.time.start)-1 do begin
    sopa = vso_get(bus[i], /RICE, out_dir=ODIR,/NODOWNLOAD ) ; Downloads the data
    fname = fnamet[i]+'_'+file_basename(sopa.fileid)
    odir1=ODIR0+file_dirname(sopa.fileid)
    if ~file_test(odir1+'/'+fname) or do_lasco[0] eq 2 then begin
      if file_test(odir1+'/'+fname) eq 0 then file_mkdir, odir1   ; creates the output directory
      spawn,'wget '+string(34B)+sopa.url[0]+string(34B)+' -O '+ odir1+'/'+fname
    endif else begin
      message,/informational, ' Existing file, it will not be downloaded...' + fname
    endelse
  endfor
 endif

 ; preprocesses files if they do not exist in ODIR+'/'+PDIR
 if do_lasco[1] gt 0 then begin
   message, 'INFO: Preprocessing files.....' ,/informational    
    for i=0, n_elements(bus.time.start)-1 do begin 
      sopa = vso_get(bus[i], /RICE, out_dir=ODIR,/NODOWNLOAD ) ; Downloads the data
      fname = file_basename(sopa.fileid)
      odir1=ODIR0+file_dirname(sopa.fileid)
      odir2=strjoin(strsplit(odir1,'level_05',/extract,/regex),PDIR)
      if file_test(odir2) eq 0 then file_mkdir, odir1+'/'+PDIR   ; creates the output directory     
      if ~file_test(odir1+'/'+PDIR+'/'+fname+'*') or do_lasco[1] eq 2 then begin
        file=file_Search(odir1+'/*'+fname+'*', count=count)
        if count eq 1 then $
          reduce_level_1, file, savedir=odir2; preprocessing the data
        endif else begin
          message,/informational, ' Existing file, it will not be preprocessed...' + fname
      endelse
    endfor
  endif
 endif

;***********************SDO/AIA using VSO
if total(do_aia) gt 0 then begin
;download
 DATES=['2011/02/15','2011/02/15']
 TIMES=['01:35:00','02:00:00']
 TIMER=1*60.  ; Time resolution [s]
 WAVE='193'  ; Wavelenght [nm]
 INSTR='aia' ; Instrument
 SOURCE='sdo'; Observatory
 LEVEL='L1'  ; Data reduction level
 ;Only for preprocessing
 PDIR='preped' ; Preprocessed data output dir is ODIR+'/'+'/'+PDIR
 FRMDIR='png_frames' ; diur to save the png frames 
 ;Only to make movie
 MNAME='movie.mvi'
 LOG=0 ; Applies ALOG10() function to image before byte scaling  
 FOV=[1536.+[0, 1024], 768.+[0, 1024]] ; size of the square FOV to use
 AVG=1 ; Num of frms to average 
 CT=[0,0,2000,0.6] ; Color table settings [table ID,range low(DN), range high (DN), gamma correction]

 ;play move
 RDIFF=0 ; set to play a runnig difference movie

 ;
 ; odris 
 ODIR = ODIR0 + '/'+SOURCE+'/'+INSTR+'/'+LEVEL+'/'+WAVE+'/'+strjoin(strsplit(DATES[0],'/',/extract))
 if file_test(ODIR) eq 0 then file_mkdir, ODIR   ; creates the output directory
 if file_test(ODIR+'/'+PDIR) eq 0 then file_mkdir, ODIR+'/'+'/'+PDIR   ; creates the output directory

 ; checks repository
 message, 'INFO: Checking repository.....' ,/informational 
  bus=vso_search(DATES[0]+' '+TIMES[0],DATES[1]+' '+TIMES[1], $ 
 	source=SOURCE,instr=INSTR, wave=WAVE, sample=TIMER) ;Gets metadata from VSO
  fname=strjoin(strsplit(bus[0].time.start,':',/extract),'_') ; only the times
  for i=1, n_elements(bus.time.start)-1 do fname=[fname,strjoin(strsplit(bus[i].time.start,':',/extract),'_')]

 ; downloads files only if they do not exist in ODIR
 if do_aia[0] gt 0 then begin
  message, 'INFO: Downloading files.....' ,/informational 
   for i=0, n_elements(fname)-1 do $ 
   	if ~file_test(ODIR+'/*'+fname[i]+'*') or do_aia[0] eq 2 then begin
        sopa = vso_get(bus[i], /RICE, out_dir=ODIR,/NODOWNLOAD ) ; Downloads the data
        spawn,'wget '+string(34B)+sopa.url[0]+string(34B)+' -O '+ODIR+'/'+fname[i]+'.fits' 
      endif else begin
        message,/informational, ' Existing file, it will not be downloaded...' + fname[i]
      endelse
 endif

 ; preprocesses files if they do not exist in ODIR+'/'+PDIR
 if do_aia[1] gt 0 then begin
   message, 'INFO: Preprocessing files.....' ,/informational    
    pfname=strjoin(strsplit(strjoin(strsplit(strjoin(strsplit(fname[0],'_',/extract)),'-',/extract)),'T',/extract),'_')
    for i=1, n_elements(bus.time.start)-1 do $ 
      pfname=[pfname,strjoin(strsplit(strjoin(strsplit(strjoin(strsplit(fname[i],'_',/extract)),'-',/extract)),'T',/extract),'_')] 
    for i=0, n_elements(pfname)-1 do begin 
    	if ~file_test(ODIR+'/'+PDIR+'/*'+pfname[i]+'*') or do_aia[1] eq 2 then begin
        bla=file_Search(ODIR+'/*'+fname[i]+'*', count=count)
        if count eq 1 then $
         aia_prep, bla, -1, outdir=ODIR+'/'+PDIR, /do_write_fits ; preprocessing the data
      endif else begin
        message,/informational, ' Existing file, it will not be preprocessed...' + fname[i]
      endelse
    endfor
  endif
  
  ;makes movie
  if do_aia[2] gt 0 then begin
    f=file_search(ODIR+'/'+PDIR+'/*.fits')
    loadct, CT[0]
    gamma_ct, CT[3],/current
   if file_test(ODIR+'/'+PDIR+'/'+FRMDIR) eq 0 then file_mkdir, ODIR+'/'+PDIR+'/'+FRMDIR   ; creates the output directory
   if do_aia[2] eq 3 then begin ; saves individual average_frames
      scc_mkmovie, f,CT[1],CT[2], cam='aia',save= ODIR+'/'+PDIR+'/'+FRMDIR+'/yyyymmdd_hhmmss_'+WAVE+'.png', /times,/dorotate, log_scl=LOG ,average_frames=AVG,coords=FOV;,/getsubfield
    endif else begin
      if ~file_search(ODIR+'/'+PDIR+'/'+MNAME) or do_aia[2] eq 2 then $ 
        scc_mkmovie, f,CT[1],CT[2], save= ODIR+'/'+PDIR+'/'+MNAME, /times,/dorotate, log_scl=LOG ,average_frames=AVG,coords=FOV;,/getsubfield
    endelse
  endif

  ; plays movie
  if do_aia[3] gt 0 then scc_playmovie, ODIR+'/'+PDIR+'/'+MNAME,running_diff=RDIFF
endif

;**********************SDO/HMI using VSO
if total(do_hmi) gt 0 then begin
 DATES=['2010/11/12','2010/11/12']
 TIMES=['12:30:00','13:30:00']
 TIMER=720.  ; Time resolution [s]. 45, 45, 45 and 90 s for Velocity, LOS B, Continuum and Full B respectivelly
 INSTR='hmi' ; Instrument
 SOURCE='sdo'; Observatory
 PHYSOBS='vector_magnetic_field'; 'LOS_magnetic_field', VSO obserbable, for more options see: https://vso.nascom.nasa.gov/cgi-bin/show_details?instrument=HMI

 ;
 ; odris 
 ODIR = ODIR0 + '/'+SOURCE+'/'+INSTR+'/'+strjoin(strsplit(DATES[0],'/',/extract))
 if file_test(ODIR) eq 0 then file_mkdir, ODIR   ; creates the output directory

 ; checks repository
 message, 'INFO: Checking repository.....' ,/informational 
  bus=vso_search(DATES[0]+' '+TIMES[0],DATES[1]+' '+TIMES[1], $ 
  source=SOURCE,instr=INSTR, sample=TIMER, physobs=PHYSOBS) ;Gets metadata from VSO
  fname=bus.time.start
  
 ; downloads files if they do not exist in ODIR
 if do_hmi[0] gt 0 then begin
  message, 'INFO: files found in repository.....' ,/informational 
   message,/informational, strjoin(fname," ; ")
   ef=file_search(ODIR+'/*.fits',count=enf)
   efname = ''
   if enf gt 0 then begin
     read_sdo,ef[0], eh, /nodata
     eh=strsplit(eh.DATE_D$OBS,'.',/extract)
     efname = eh[0]
     for j=1, enf-1 do begin 
      read_sdo,ef[j], eh, /nodata
      eh = strsplit(eh.DATE_D$OBS,'.',/extract)
      efname=[efname,eh[0]]
     endfor
  endif
   for i=0, n_elements(fname)-1 do begin
    if total(strmatch(efname, fname[i])) eq 0 or do_hmi[0] eq 2 then begin
        sopa = vso_get(bus[i], /RICE, out_dir=ODIR,/NODOWNLOAD ) ; Downloads the data
        spawn,'wget '+string(34B)+sopa.url[0]+string(34B)+' -O '+ODIR+'/'+fname[i]+'.fits' 
        stop ; chekc the next line!!!!
    endif else begin
        message,/informational, ' Existing file, it will not be downloaded...' + bus[i].info + ' ; ' + fname[i]
    endelse
  endfor
 endif

endif

;**********************SHARPS using JSOC
if total(do_sharp) gt 0 then begin
 DATES=['2010/11/11','2010/11/11']
 TIMES=['04:50:00','05:20:00']
 INSTR='hmi' ; Instrument
 SOURCE='sdo'; Observatory
 DATA_SERIES='hmi.sharp_cea_720s'; for more options see: http://jsoc.stanford.edu/doc/data/hmi/sharp/sharp.htm
 ;optional (comment out to avoid using) ; For more info see https://link.springer.com/content/pdf/10.1007%2Fs11207-014-0529-3.pdf
 SEGMENT=['Bp','Bt','Br','Bp_err','Bt_err','Br_err','conf_disambig','magnetogram'] ; data products
 ;KEWYWORDS=['totusjh','absnjzh','meanjzh'] ; header kewywords
 HARPNUM=245 ; NOAA AR number
 ;bla='hmi.sharp_cea_720s[245][2010.11.11_01:20_TAI-2010.11.11_14:59_TAI]{magnetogram,Bp,Bt,Br,Bp_err,Bt_err,Br_err,conf_disambig}'
 ; odris 
 ODIR = ODIR0 + '/'+SOURCE+'/'+INSTR+'/sharps/'+strjoin(strsplit(DATES[0],'/',/extract))
 if file_test(ODIR) eq 0 then file_mkdir, ODIR   ; creates the output directory
 
 ndates=igl_dateconv([DATES[0]+' '+TIMES[0], DATES[1]+' '+TIMES[1]],[2,3]) ; changes dir format

 ; Downloads data from JSOC one segment at a time
 for i=0, n_elements(SEGMENT)-1 do begin
   ; gets metadata
   ;see http://www.heliodocs.com/php/xdoc_print.php?file=$SSW/vobs/ontology/idl/jsoc/ssw_jsoc_time2data.pro for extra parameters
   ssw_jsoc_time2data, ndates[0], ndates[1], ind, ds=DATA_SERIES, /copy_only, parent_out=ODIR $
                       ,keywords=KEWYWORDS, segment=SEGMENT[i], harpnum=HARPNUM,/slient
   nofname = ODIR+'/'+DATA_SERIES+'.'+strtrim(ind.harpnum,2)+'.'+igl_dateconv(ind.t_obs,[4,5])+'.'+SEGMENT[i]+'.fits'

   ;Downloads the segment if any file is missing
   ef=file_search(ODIR+'/*.fits')
   flag=0
   for j=0, n_elements(nofname)-1 do flag += total(strmatch(ef, nofname[j])) gt 0
   if (flag eq n_elements(nofname)) and (do_sharp[0] ne 2) then begin
      message,/informational, ' WARNING, all files for segment '+SEGMENT[i]+' already exist, skipping downloads'     
   endif else begin
     ssw_jsoc_time2data, ndates[0], ndates[1], ind, img, ds=DATA_SERIES, /copy_only, parent_out=ODIR $
                         ,keywords=KEWYWORDS, segment=SEGMENT[i], harpnum=HARPNUM, fnames_uncomp=ofname, /comp_delete,/silent
     ; reneames the uncompressed oputput files to match the official nameing and avoid overwriting.                          
     nofname = file_dirname(ofname)+'/'+DATA_SERIES+'.'+strtrim(ind.harpnum,2)+'.'+igl_dateconv(ind.t_obs,[4,5])+'.'+SEGMENT[i]+'.fits'
     file_move, ofname,nofname  
   endelse
 endfor
endif

;******************************STEREOA/EUVI using vso_search and html to download
if total(do_euvia) gt 0 then begin
;download
 DATES=['2011/02/15','2011/02/15']
 TIMES=['01:45:00','02:00:00']
 TIMER=1.*60.   ; Time resolution [s]
 WAVE='195'  ; Wavelenght [nm]
 INSTR='euvi' ; Instrument
 SOURCE='stereo_a'; Observatory
 LEVEL='L0'  ; Data reduction level
 ;preprocessing
 PDIR='preped' ; Preprocessed data output dir is ODIR+'/'+PDIR
  ;make movie
 MNAME='movie.mvi'
 LOG=0 ; Applies ALOG10() function to image before byte scaling  
 FOV=[256.+[0, 512], 350.+[0, 512]] ; size of the square FOV to use, comment ot use all
 ;AVG=1 ; Num of frms to average 
 ;play move
 RDIFF=1 ; set to play a runnig difference movie


 ODIR = ODIR0 + '/stereo/secchi/'+LEVEL+'/a/img/'+INSTR+'/'+WAVE
 ; checks repository
 message, 'INFO: Checking repository.....' ,/informational 
 bus=vso_search(DATES[0]+' '+TIMES[0],DATES[1]+' '+TIMES[1], $ 
 source=SOURCE,instr=INSTR, wave=WAVE, sample=TIMER) ;Gets metadata from VSO
 fname=strjoin(strsplit(bus[0].time.start,':',/extract),'_') ; only the times
 for i=1, n_elements(bus.time.start)-1 do fname=[fname,strjoin(strsplit(bus[i].time.start,':',/extract),'_')]
 fname0=igl_dateconv(fname,[7,6])
 fname=fname0+ '_n4euA.fts'

 ; downloads files only if they do not exist in ODIR
 if do_euvia[0] gt 0 then begin
  message, 'INFO: Downloading files.....' ,/informational 
  for i=0, n_elements(fname)-1 do begin
    date=strsplit(fname[i],'_',/extract)
    date=date[0]
    odir1=ODIR+'/'+date
    if ~file_test(odir1+'/'+fname[i]) or do_euvia[0] eq 2 then begin
      if file_test(odir1) eq 0 then file_mkdir, odir1   ; creates the output directory
       sopa = vso_get(bus[i], /RICE, out_dir=ODIR,/NODOWNLOAD ) ; Downloads the data
       spawn,'wget '+string(34B)+sopa.url[0]+string(34B)+' -O '+odir1+'/'+fname[i]+'.fits'    
    endif else begin
      message,/informational, ' Existing file, it will not be downloaded...' + fname[i]
    endelse
  endfor
 endif

 ; preprocesses files if they do not exist in ODIR+'/'+PDIR
 if do_euvia[1] gt 0 then begin
   message, 'INFO: Preprocessing files.....' ,/informational    
    fname=fname0+ '_14euA.fts'
    for i=0, n_elements(fname)-1 do begin 
      date=strsplit(fname[i],'_',/extract)
      date=date[0]
      odir1=ODIR+'/'+date 
      if file_test(odir1+'/'+PDIR) eq 0 then file_mkdir, odir1+'/'+PDIR   ; creates the output directory     
      if ~file_test(odir1+'/'+PDIR+'/'+fname[i]) or do_euvia[1] eq 2 then begin
        file=file_Search(odir1+'/*'+fname0[i]+'*', count=count)
        if count eq 1 then $
          ;im=sccreadfits(file,h)
          ; preprocessing , see https://hesperia.gsfc.nasa.gov/ssw/stereo/secchi/doc/secchi_prep.html
          secchi_prep,file, /rotate_on, /write_fts, savepath=file_dirname(file)+'/'+PDIR, outssize=512
        endif else begin
          message,/informational, ' Existing file, it will not be preprocessed...' + fname[i]
      endelse
    endfor
  endif

  ;makes movie
  mp= ODIR+'/'+strjoin(strsplit(DATES[0],'/',/extract))+'/'+PDIR+'/'+MNAME
  if do_euvia[2] eq 1 and ~file_search(mp) or do_euvia[2] eq 2 then begin
    f=''
    for i=0, n_elements(fname0)-1 do begin 
      date=strsplit(fname0[i],'_',/extract)
      f=[f,ODIR+'/'+date[0]+'/'+PDIR+'/'+fname0[i]+ '_14euA.fts']
    endfor
    f=f[1:*]
    d=mrdfits(f[0])
    sz=size(d)
    if ~keyword_set(FOV) then FOV=[0,sz[1]-1,0,sz[2]-1]
    ;m=median(d[FOV[0]:FOV[1],FOV[2]:FOV[3]])
    ;sd=stddev(d[FOV[0]:FOV[1],FOV[2]:FOV[3]])
    scc_mkmovie, f,0,0,'euvi',/AUTO,/DCOLOR ,save=mp,/dorotate, log_scl=LOG ,average_frames=AVG,coords=FOV;,/getsubfield
  endif

  ; plays movie
  if do_euvia[3] gt 0 then scc_playmovie, mp,running_diff=RDIFF, /times

endif


;********************************STEREOB/EUVI using VSO
if total(do_euvib) gt 0 then begin
;download
 DATES=['2011/02/15','2011/02/15']
 TIMES=['01:45:00','02:20:00']
 TIMER=1.*60.   ; Time resolution [s]
 WAVE='195'  ; Wavelenght [nm]
 INSTR='euvi' ; Instrument
 SOURCE='stereo_b'; Observatory
 LEVEL='L0'  ; Data reduction level
 ;preprocessing
 PDIR='preped' ; Preprocessed data output dir is ODIR+'/'+PDIR
  ;make movie
 MNAME='movie.mvi'
 LOG=0 ; Applies ALOG10() function to image before byte scaling  
 FOV=[256.+[0, 512], 350.+[0, 512]] ; size of the square FOV to use, comment ot use all
 ;AVG=1 ; Num of frms to average 
 ;play move
 RDIFF=1 ; set to play a runnig difference movie


 ODIR = ODIR0 + '/stereo/secchi/'+LEVEL+'/b/img/'+INSTR+'/'+WAVE
 ; checks repository
 message, 'INFO: Checking repository.....' ,/informational 
 bus=vso_search(DATES[0]+' '+TIMES[0],DATES[1]+' '+TIMES[1], $ 
 source=SOURCE,instr=INSTR, wave=WAVE, sample=TIMER) ;Gets metadata from VSO
 fname=strjoin(strsplit(bus[0].time.start,':',/extract),'_') ; only the times
 for i=1, n_elements(bus.time.start)-1 do fname=[fname,strjoin(strsplit(bus[i].time.start,':',/extract),'_')]
 fname0=igl_dateconv(fname,[7,6])
 fname=fname0+ '_n4euB.fts'

 ; downloads files only if they do not exist in ODIR
 if do_euvib[0] gt 0 then begin
  message, 'INFO: Downloading files.....' ,/informational 
  for i=0, n_elements(fname)-1 do begin
    date=strsplit(fname[i],'_',/extract)
    date=date[0]
    odir1=ODIR+'/'+date
    if ~file_test(odir1+'/'+fname[i]) or do_euvib[0] eq 2 then begin
      if file_test(odir1) eq 0 then file_mkdir, odir1   ; creates the output directory
       sopa = vso_get(bus[i], /RICE, out_dir=ODIR,/NODOWNLOAD ) ; Downloads the data
       spawn,'wget '+string(34B)+sopa.url[0]+string(34B)+' -O '+odir1+'/'+fname[i]+'.fits'       
    endif else begin
      message,/informational, ' Existing file, it will not be downloaded...' + fname[i]
    endelse
  endfor
 endif

 ; preprocesses files if they do not exist in ODIR+'/'+PDIR
 if do_euvib[1] gt 0 then begin
   message, 'INFO: Preprocessing files.....' ,/informational    
    fname=fname0+ '_14euB.fts'
    for i=0, n_elements(fname)-1 do begin 
      date=strsplit(fname[i],'_',/extract)
      date=date[0]
      odir1=ODIR+'/'+date 
      if file_test(odir1+'/'+PDIR) eq 0 then file_mkdir, odir1+'/'+PDIR   ; creates the output directory     
      if ~file_test(odir1+'/'+PDIR+'/'+fname[i]) or do_euvib[1] eq 2 then begin
        file=file_Search(odir1+'/*'+fname0[i]+'*', count=count)
        if count eq 1 then $
          ;im=sccreadfits(file,h)
          ; preprocessing , see https://hesperia.gsfc.nasa.gov/ssw/stereo/secchi/doc/secchi_prep.html
          secchi_prep,file, /rotate_on, /write_fts, savepath=file_dirname(file)+'/'+PDIR, outssize=512
        endif else begin
          message,/informational, ' Existing file, it will not be preprocessed...' + fname[i]
      endelse
    endfor
  endif

  ;makes movie
  mp= ODIR+'/'+strjoin(strsplit(DATES[0],'/',/extract))+'/'+PDIR+'/'+MNAME
  if do_euvib[2] eq 1 and ~file_search(mp) or do_euvib[2] eq 2 then begin
    f=''
    for i=0, n_elements(fname0)-1 do begin 
      date=strsplit(fname0[i],'_',/extract)
      f=[f,ODIR+'/'+date[0]+'/'+PDIR+'/'+fname0[i]+ '_14euA.fts']
    endfor
    f=f[1:*]
    d=mrdfits(f[0])
    sz=size(d)
    if ~keyword_set(FOV) then FOV=[0,sz[1]-1,0,sz[2]-1]
    ;m=median(d[FOV[0]:FOV[1],FOV[2]:FOV[3]])
    ;sd=stddev(d[FOV[0]:FOV[1],FOV[2]:FOV[3]])
    scc_mkmovie, f,0,0,'euvi',/AUTO,/DCOLOR ,save=mp,/dorotate, log_scl=LOG ,average_frames=AVG,coords=FOV;,/getsubfield
  endif

  ; plays movie
  if do_euvib[3] gt 0 then scc_playmovie, mp,running_diff=RDIFF, /times

endif

;********************************STEREOB/CORB using VSO
if total(do_corb) gt 0 then begin
;download
 DATES=['2011/02/15','2011/02/15']
 TIMES=['02:15:00','03:45:00']
 TIMER=1.*60.   ; Time resolution [s]
 INSTR='cor2' ; Instrument
 PHYSOBS='intensity'
 SOURCE='stereo_b'; Observatory
 ;preprocessing
 PDIR='level1' ; Preprocessed data output dir is ODIR+'/'+PDIR

;*****
 ; checks repository
 message, 'INFO: Checking repository.....' ,/informational 
 bus=vso_search(DATES[0]+' '+TIMES[0],DATES[1]+' '+TIMES[1], $ 
 source=SOURCE,instr=INSTR, sample=TIMER, physobs=PHYSOBS) ;Gets metadata from VSO

 ; downloads files only if they do not exist in ODIR
 if do_corb[0] gt 0 then begin
  message, 'INFO: Downloading files.....' ,/informational 
  for i=0, n_elements(bus.time.start)-1 do begin
    sopa = vso_get(bus[i], /RICE, out_dir=ODIR,/NODOWNLOAD ) ; Downloads the data
    fname = file_basename(sopa.fileid)
    odir1=ODIR0+'/stereo/'+file_dirname(sopa.fileid)
    if ~file_test(odir1+'/'+fname) or do_corb[0] eq 2 then begin
      if file_test(odir1) eq 0 then file_mkdir, odir1   ; creates the output directory
      spawn,'wget '+string(34B)+sopa.url[0]+string(34B)+' -O '+ odir1+'/'+fname
    endif else begin
      message,/informational, ' Existing file, it will not be downloaded...' + fname
    endelse
  endfor
 endif

 ; preprocesses files if they do not exist in ODIR+'/'+PDIR
 if do_corb[1] gt 0 then begin
   message, 'INFO: Preprocessing files.....' ,/informational    
    for i=0, n_elements(bus.time.start)-1 do begin 
      sopa = vso_get(bus[i], /RICE, out_dir=ODIR,/NODOWNLOAD ) ; Downloads the data
      fname = file_basename(sopa.fileid)
      odir1=ODIR0+'/stereo/'+file_dirname(sopa.fileid)
      if file_test(odir1+'/'+PDIR) eq 0 then file_mkdir, odir1+'/'+PDIR   ; creates the output directory     
      if ~file_test(odir1+'/'+PDIR+'/'+fname+'*') or do_corb[1] eq 2 then begin
        file=file_Search(odir1+'/*'+fname+'*', count=count)
        if count eq 1 then $
          secchi_prep,file, /rotate_on, /write_fts, savepath=file_dirname(file)+'/'+PDIR, outssize=512
        endif else begin
          message,/informational, ' Existing file, it will not be preprocessed...' + fname
      endelse
    endfor
  endif
endif

;********************************STEREOB/CORA using VSO
if total(do_cora) gt 0 then begin
;download
 DATES=['2011/02/15','2011/02/15']
 TIMES=['02:15:00','03:45:00']
 TIMER=1.*60.   ; Time resolution [s]
 INSTR='cor2' ; Instrument
 PHYSOBS='intensity'
 SOURCE='stereo_a'; Observatory
 ;preprocessing
 PDIR='level1' ; Preprocessed data output dir is ODIR+'/'+PDIR


;*****
 ; checks repository
 message, 'INFO: Checking repository.....' ,/informational 
 bus=vso_search(DATES[0]+' '+TIMES[0],DATES[1]+' '+TIMES[1], $ 
 source=SOURCE,instr=INSTR, sample=TIMER, physobs=PHYSOBS) ;Gets metadata from VSO

 ; downloads files only if they do not exist in ODIR
 if do_cora[0] gt 0 then begin
  message, 'INFO: Downloading files.....' ,/informational 
  for i=0, n_elements(bus.time.start)-1 do begin
    sopa = vso_get(bus[i], /RICE, out_dir=ODIR,/NODOWNLOAD ) ; Downloads the data
    fname = file_basename(sopa.fileid)
    odir1=ODIR0+'/stereo/'+file_dirname(sopa.fileid)
    if ~file_test(odir1+'/'+fname) or do_cora[0] eq 2 then begin
      if file_test(odir1) eq 0 then file_mkdir, odir1   ; creates the output directory
      spawn,'wget '+string(34B)+sopa.url[0]+string(34B)+' -O '+ odir1+'/'+fname
    endif else begin
      message,/informational, ' Existing file, it will not be downloaded...' + fname
    endelse
  endfor
 endif

 ; preprocesses files if they do not exist in ODIR+'/'+PDIR
 if do_cora[1] gt 0 then begin
   message, 'INFO: Preprocessing files.....' ,/informational    
    for i=0, n_elements(bus.time.start)-1 do begin 
      sopa = vso_get(bus[i], /RICE, out_dir=ODIR,/NODOWNLOAD ) ; Downloads the data
      fname = file_basename(sopa.fileid)
      odir1=ODIR0+'/stereo/'+file_dirname(sopa.fileid)
      if file_test(odir1+'/'+PDIR) eq 0 then file_mkdir, odir1+'/'+PDIR   ; creates the output directory     
      if ~file_test(odir1+'/'+PDIR+'/'+fname+'*') or do_cora[1] eq 2 then begin
        file=file_Search(odir1+'/*'+fname+'*', count=count)
        if count eq 1 then $
          secchi_prep,file, /rotate_on, /write_fts, savepath=file_dirname(file)+'/'+PDIR, outssize=512
        endif else begin
          message,/informational, ' Existing file, it will not be preprocessed...' + fname
      endelse
    endfor
  endif
endif

message, "ALL TASKS COMPLETED :-D" ,/informational
end