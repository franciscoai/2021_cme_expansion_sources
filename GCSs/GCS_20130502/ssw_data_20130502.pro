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
;*************************************************************************
;-
pro ssw_data_20130502
@gcs_config

time=SYSTIME(/SECONDS) ; just to check the script running time

; General Constants
 ODIR=DATA_PATH ; Main output data directory

;Tasks
;For each Instrument, set to 1 to:
;[download, preprocess, make movie, play the movie] (not all are avilable for all inst.)
;set to 0 to do nothing and set to 2 to overwrite Ouput (when relevant).
do_aia=1*[1,1,0,0]   
do_hmi=[0]
do_sharp=[0]   
do_euvia=1 ; not finished yet
do_euvib=1 ; not finished yet

;**********Main
;
;***********************SDO/AIA using VSO
if total(do_aia) gt 0 then begin
;download
 DATES=['2013/05/02','2013/05/02']
 TIMES=['05:00:00','05:25:00']
 TIMER=5. *60.   ; Time resolution [s]
 WAVE='193'  ; Wavelenght [nm]
 INSTR='AIA' ; Instrument
 SOURCE='SDO'; Observatory
 LEVEL='L1'  ; Data reduction level
 ;preprocessing
 PDIR='/preped' ; Preprocessed data output dir is ODIR+PDIR
 ;make movie
   MNAME='movie.mvi'
   LOG=0 ; Applies ALOG10() function to image before byte scaling  
   FOV=[1536.+[0, 1024], 768.+[0, 1024]] ; soze of the square FOV to use
   AVG=1 ; Num of frms to average 
   ;play move
   RDIFF=0 ; set to play a runnig difference movie

 ;
 ; odris 
 ODIR += '/'+SOURCE+'/'+INSTR+'/'+LEVEL+'/'+WAVE+'/'+strjoin(strsplit(DATES[0],'/',/extract))
 if file_test(ODIR) eq 0 then file_mkdir, ODIR   ; creates the output directory
 if file_test(ODIR+PDIR) eq 0 then file_mkdir, ODIR+PDIR   ; creates the output directory

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
   	if ~file_test(ODIR+'/*'+fname[i]+'*') or do_aia[0] eq 2 then $
        sopa = vso_get(bus[i], /RICE, out_dir=ODIR) $ ; Downloads the data
      else $
        print, 'ssw_data: Existing file, it will not be downloaded...' + fname[i]
 endif

 ; preprocesses files if they do not exist in ODIR+PDIR
 if do_aia[1] gt 0 then begin
   message, 'INFO: Preprocessing files.....' ,/informational 
    pfname=strjoin(strsplit(strjoin(strsplit(strjoin(strsplit(fname[0],'_',/extract)),'-',/extract)),'T',/extract),'_')
    for i=1, n_elements(bus.time.start)-1 do $ 
      pfname=[pfname,strjoin(strsplit(strjoin(strsplit(strjoin(strsplit(fname[i],'_',/extract)),'-',/extract)),'T',/extract),'_')]
    for i=0, n_elements(pfname)-1 do begin 
    	if ~file_test(ODIR+PDIR+'/*'+pfname[i]+'*') or do_aia[1] eq 2 then begin
        bla=file_Search(ODIR+'/*'+fname[i]+'*', count=count)
        if count eq 1 then $
         aia_prep, bla, -1, outdir=ODIR+PDIR, /do_write_fits ; preprocessing the data
      endif else begin
        print, 'ssw_data: Existing file, it will not be preprocessed...' + fname[i]
      endelse
    endfor
  endif
  
  ;makes movie
  if do_aia[2] gt 0 then begin
    f=file_search(ODIR+PDIR+'/*'+pfname[0]+'*')
    d=mrdfits(f[0])
    m=mean(d)
    sd=stddev(d)
    f=file_search(ODIR+PDIR+'/*.fits')
    if ~file_search(ODIR+PDIR+'/'+MNAME) or do_aia[2] eq 2 then $ 
      scc_mkmovie, f,m-1*sd,m+1*sd,'aia', save= ODIR+PDIR+'/'+MNAME, /times,/dorotate, log_scl=LOG ,average_frames=AVG,coords=FOV;,/getsubfield
  endif

  ; plays movie
  if do_aia[3] gt 0 then scc_playmovie, ODIR+PDIR+'/'+MNAME,running_diff=RDIFF
endif

;**********************SDO/HMI using VSO
if total(do_hmi) gt 0 then begin
 DATES=['2010/11/11','2010/11/11']
 TIMES=['06:20:00','06:50:00']
 TIMER=720.   ; Time resolution [s]. 45, 45, 45 and 90 s for Velocity, LOS B, Continuum and Full B respectivelly
 INSTR='hmi' ; Instrument
 SOURCE='SDO'; Observatory
 PHYSOBS='vector_magnetic_field'; VSO obserbable, for more options see: https://vso.nascom.nasa.gov/cgi-bin/show_details?instrument=HMI

 ;
 ; odris 
 ODIR += '/'+SOURCE+'/'+INSTR+'/'+strjoin(strsplit(DATES[0],'/',/extract))
 if file_test(ODIR) eq 0 then file_mkdir, ODIR   ; creates the output directory

 ; checks repository
 message, 'INFO: Checking repository.....' ,/informational 
  bus=vso_search(DATES[0]+' '+TIMES[0],DATES[1]+' '+TIMES[1], $ 
  source=SOURCE,instr=INSTR, sample=TIMER, physobs=PHYSOBS) ;Gets metadata from VSO
  fname=bus.time.start
  
 ; downloads files if they do not exist in ODIR
 if do_hmi[0] gt 0 then begin
  message, 'INFO: files found in repository.....' ,/informational 
   print, strjoin(fname," ; ")
   ef=file_search(ODIR+'/*.fits',count=enf)
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
        sopa = vso_get(bus[i], /RICE, out_dir=ODIR)  ; Downloads the data 
    endif else begin
        print, 'ssw_data: Existing file, it will not be downloaded...' + bus[i].info + ' ; ' + fname[i]
    endelse
  endfor
 endif

endif

;**********************SHARPS using JSOC
if total(do_sharp) gt 0 then begin
 DATES=['2010/11/11','2010/11/11']
 TIMES=['06:50:00','07:20:00']
 INSTR='hmi' ; Instrument
 SOURCE='SDO'; Observatory
 DATA_SERIES='hmi.sharp_cea_720s'; for more options see: http://jsoc.stanford.edu/doc/data/hmi/sharp/sharp.htm
 ;optional (comment out to avoid using) ; For more info see https://link.springer.com/content/pdf/10.1007%2Fs11207-014-0529-3.pdf
 SEGMENT=['Bp','Bt','Br','Bp_err','Bt_err','Br_err','conf_disambig','magnetogram'] ; data products
 ;KEWYWORDS=['totusjh','absnjzh','meanjzh'] ; header kewywords
 HARPNUM=245 ; NOAA AR number
 ;bla='hmi.sharp_cea_720s[245][2010.11.11_01:20_TAI-2010.11.11_14:59_TAI]{magnetogram,Bp,Bt,Br,Bp_err,Bt_err,Br_err,conf_disambig}'
 ; odris 
 ODIR += '/'+SOURCE+'/'+INSTR+'/'+strjoin(strsplit(DATES[0],'/',/extract))
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
      print, 'ssw_data: WARNING, all files for segment '+SEGMENT[i]+' already exist, skipping downloads'     
   endif else begin
     ssw_jsoc_time2data, ndates[0], ndates[1], ind, img, ds=DATA_SERIES, /copy_only, parent_out=ODIR $
                         ,keywords=KEWYWORDS, segment=SEGMENT[i], harpnum=HARPNUM, fnames_uncomp=ofname, /comp_delete,/silent
     ; reneames the uncompressed oputput files to match the official nameing and avoid overwriting.                          
     nofname = file_dirname(ofname)+'/'+DATA_SERIES+'.'+strtrim(ind.harpnum,2)+'.'+igl_dateconv(ind.t_obs,[4,5])+'.'+SEGMENT[i]+'.fits'
     file_move, ofname,nofname  
   endelse
 endfor
endif

;******************************STEREOA/EUVI
if do_euvia then begin
 DATE='20130502'
 FILES=['052030','051030','045530']+'_n4euA'
 EUVI_REPO='https://stereo-ssc.nascom.nasa.gov/data/ins_data/secchi/L0/a/img/euvi'
 EUVI_ODIR=ODIR+'/stereo/secchi/L0/a/img/euvi/'+DATE

 ;
 EUVI_SRC=EUVI_REPO+'/'+DATE+'/'+DATE+'_'+FILES+'.fts'
 for i=0, n_elements(EUVI_SRC)-1 do begin
  if file_test(EUVI_ODIR+'/'+DATE+'_'+FILES[i]+'.fts') eq 0 then $
    spawn,'wget '+string(34B)+EUVI_SRC[i]+string(34B)+' -P '+ EUVI_ODIR ; downloads file if it does not exist
 endfor
endif


;********************************STEREOB/EUVI
if do_euvib then begin
 DATE='20130502'
 FILES=['052030','051030','045530']+'_n4euB'
 EUVI_REPO='https://stereo-ssc.nascom.nasa.gov/data/ins_data/secchi/L0/b/img/euvi'
 EUVI_ODIR=ODIR+'/stereo/secchi/L0/b/img/euvi/'+DATE

 ;
 EUVI_SRC=EUVI_REPO+'/'+DATE+'/'+DATE+'_'+FILES+'.fts'
 for i=0, n_elements(EUVI_SRC)-1 do begin
  if file_test(EUVI_ODIR+'/'+DATE+'_'+FILES[i]+'.fts') eq 0 then $
    spawn,'wget '+string(34B)+EUVI_SRC[i]+string(34B)+' -P '+ EUVI_ODIR ; downloads file if it does not exist
 endfor
endif

print, 'ssw_data: DONE in '+strtrim(string(SYSTIME(/SECONDS)-time),2)+' s.... :-P'
end