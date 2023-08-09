;+
;*************************************************************************
; Francisco A. Iglesias - UTN-FRM/GEAA - franciscoaiglesias@frm.utn.edu.ar
;
; Reads all the .sav files and plots the timestamps of the images used 
; in the fits
; 
; History:
;	20180411 - First version
;*************************************************************************
;-
function timeline_all_search, path
; Reads the .pro files  in path directory, and extracts the dates of all the image files used.
; The dates are returned in Julian 
; See a file in the GCSs directory for a .pro example
@gcs_config
message, 'Reading file...' + path ,/informational
times= replicate({inst:' ', dat0:0d, dat1:0d}, 6)
nkw=n_elements(kw)
openr,lun,path,/get_lun      
rline=''
readf,lun,rline                       ; First line is a comment !!!!!
while not eof(lun) do begin           ; read line by line, only the valid are saved
readf,lun,rline
bla=strpos(rline,'=')
if bla gt -1 then begin 
 rkw=strsplit(rline,'=',/extract) ; get the variable name (up to the = sign)
 rkw=strtrim(rkw[0],2)
 case rkw of
	 ; list of variables to read
	 'ima': begin
	 			times[0].inst='sta/'+strmid(rline,strpos(rline,'cor'),4)
	 			bla=strsplit(rline,'/',/extract)
	 			times[0].dat1=date_conv(igl_dateconv(strmid(bla[n_elements(bla)-1],0,15),[6,0]),'J')
	 		end
	 'imaprev': begin
	 			bla=strsplit(rline,'/',/extract)
	 			times[0].dat0=date_conv(igl_dateconv(strmid(bla[n_elements(bla)-1],0,15),[6,0]),'J')
	 		end
	 'imb': begin
	 			times[1].inst='stb/'+strmid(rline,strpos(rline,'cor'),4)
	 			bla=strsplit(rline,'/',/extract)
	 			times[1].dat1=date_conv(igl_dateconv(strmid(bla[n_elements(bla)-1],0,15),[6,0]),'J')
	 		end
	 'imbprev': begin
	 			bla=strsplit(rline,'/',/extract)
	 	 		times[1].dat0=date_conv(igl_dateconv(strmid(bla[n_elements(bla)-1],0,15),[6,0]),'J')
	 	 	end
	 'lasco1': begin
	 	       if strpos(rline,'lascopath') gt -1 then begin
		 			times[2].inst='lasco/cor2'
		 			bla=LASCO_PATH+strmid(rline,strpos(rline,'/'),strpos(rline,',')-strpos(rline,'/')-1)
		 			lascohdr=headfits(bla)
		 			lascohdr=lasco_fitshdr2struct(lascohdr)
		 			times[2].dat1=date_conv(lascohdr.date_obs,'J')
		 		endif	
	 		end
	 'lasco0':	begin
	 	       if strpos(rline,'lascopath') gt -1 then begin	 	
		 			bla=strmid(rline,strpos(rline,'/'),strpos(rline,',')-strpos(rline,'/')-1)
		 			lascohdr=headfits(LASCO_PATH+bla)
		 			lascohdr=lasco_fitshdr2struct(lascohdr)
		 			times[2].dat0=date_conv(lascohdr.date_obs,'J')
		 		endif
	 		end
	 'imeuvib': begin
	 	        bla=strmid(rline,strpos(rline,'/')+1,3)
	 	        if bla eq 'sdo' then begin
	 	        	times[3].inst='sdo/aia'
	 			    times[3].dat0=date_conv(igl_dateconv(strmid(rline,strpos(rline,'preped/')+10,15),[6,0]),'J')
	 		    endif else begin
	 				times[3].inst='sdo/euvi'
	 				times[3].dat0=date_conv(igl_dateconv(strmid(rline,strpos(rline,'preped/')+7,15),[6,0]),'J')
	 		    endelse
	 		end
	 'imeuvia': begin
	 	        bla=strmid(rline,strpos(rline,'/')+1,3)
	 	        if bla eq 'sdo' then begin
	 	        	times[4].inst='sdo/aia'
	 			    times[4].dat0=date_conv(igl_dateconv(strmid(rline,strpos(rline,'preped/')+10,15),[6,0]),'J')
	 		    endif else begin
	 				times[4].inst='sdo/euvi'
	 				times[4].dat0=date_conv(igl_dateconv(strmid(rline,strpos(rline,'preped/')+7,15),[6,0]),'J')
	 		    endelse
	 		end
	 'SDATE': begin
	 			times[0].inst='sta/euvi'
	 			times[1].inst='stb/euvi'
	 			bla=strsplit(rline,"\'",/extract)
	 			times[0].dat1=bla[1]
	 		end
	 'STIME': begin
	 			bla=strsplit(rline,"\'",/extract)
	 			times[0].dat1=date_conv(igl_dateconv(string(times[0].dat1,format='(I8)')+'_'+bla[1],[6,0]),'J')
	 			times[1].dat1=times[0].dat1
	 		end
	 'SDATE0': begin
	 			bla=strsplit(rline,"\'",/extract)
	 			times[0].dat0=bla[1]
	 		end
	 'STIME0': begin
	 			bla=strsplit(rline,"\'",/extract)
	 			times[0].dat0=date_conv(igl_dateconv(string(times[0].dat0,format='(I8)')+'_'+bla[1],[6,0]),'J')
	 			times[1].dat0=times[0].dat0
	 		end
	 'DATE0': begin
	 			times[2].inst='sdo/aia'
	 			bla=strsplit(rline,"\'",/extract)
	 			times[2].dat0=bla[1]+bla[3]+bla[5]
	 		end
	 'TIME0': begin
	 			bla=strsplit(rline,"\'",/extract)
	 			times[2].dat0=date_conv(igl_dateconv(string(times[2].dat0,format='(I8)')+'_'+bla[1]+bla[3]+bla[5],[6,0]),'J')
	 		end
	 'DATE': begin
	 			bla=strsplit(rline,"\'",/extract)
	 			times[2].dat1=bla[1]+bla[3]+bla[5]
	 		end
	 'TIME': begin
	 			bla=strsplit(rline,"\'",/extract)
	 			times[2].dat1=date_conv(igl_dateconv(string(times[2].dat1,format='(I8)')+'_'+bla[1]+bla[3]+bla[5],[6,0]),'J')
	 		end
	 else:
  endcase
endif
endwhile
free_lun, lun
return, times
end

pro  timeline_all
@gcs_config
;OPATH=DATA_PATH+'/Polar_Observations/Polar_Documents/francisco/GCSs/timeline'
OPATH=GCSD_PATH + '/plot_time_stamps'

dirs=find_all_dir(GCSD_PATH+'/GCS*')
;output path
file_mkdir, file_dirname(dirs[0]) + '/timeline/'
file_chmod, file_dirname(dirs[0]) + '/timeline/' , CHMOD_DIR; prosper Linux permissions for directories

;plots settings
mydevice = !D.NAME
SET_PLOT, 'PS'
!P.Multi = 0
FONT=4.5
MARGIN=[1,1,1,1]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
SYM=[-2,-6,-4,-5]

;loops on all dirs
nd=n_elements(dirs)
for j=0, nd-1 do begin
 DEVICE, FILENAME=OPATH+'/'+file_basename(dirs[j])+'.ps', /LANDSCAPE ,/color
 LOADCT,38 ; rainbow + white color table
 f=file_search(dirs[j]+'/','FMw*.pro')
 nf=n_elements(f)
 times=timeline_all_search(f[0])
 colors=fltarr(n_elements(times))
 for k=1, nf-1 do begin
 	bla=timeline_all_search(f[k])
 	times=[times,bla]
 	if k mod 2 then col=256./17 * k/2. else col=128+256/17.* k/2.
    colors=[colors,fltarr(n_elements(bla)) + col]
 endfor
 ;plot the time line of the current directory
 ini=min(times[where(times.dat0 ne 0d)].dat0)
 for i=0, n_elements(times)-1 do begin
	INST=['sta/euvi','stb/euvi','sdo/aia','sdo/euvi','sta/cor1','stb/cor1','lasco/cor2', 'sta/cor2','stb/cor2']
	bla=where(strmatch(INST,times[i].inst) eq 1)
	if bla gt -1 then begin 
		if i eq 0 then begin
			if times[i].dat1 eq 0d then t1=times[i].dat0 else t1=times[i].dat1
			plot,([times[i].dat0,t1]-ini)*24d  , [bla,bla] , psym=-2, linestyle=0,title= 'Time line for '+file_basename(dirs[j]), $
	        xtitle='Time from '+date_conv(ini,'FITS')+' [h]', ytitle='Instrument', yrange=[-1,9],ytickname=[' ', INST, ' '], ytickinterval=1, $ 
	        xrange=[0,max(times.dat1)-ini]*24d, font=FONT, xmargin=MARGIN[0:1], ymargin=MARGIN[2:3], charthick=1.6, thick=3
	    endif else begin
	    	if times[i].dat1 eq 0d then t1=times[i].dat0 else t1=times[i].dat1
            oplot, ([times[i].dat0,t1]-ini)*24d  , [bla,bla], psym=-2, color= colors[i], thick=3
            oplot, [0,10]*24d, [2.5,2.5], linestyle=3, thick=4
            oplot, [0,10]*24d, [5.5,5.5], linestyle=3, thick=4
	    endelse
	endif
 endfor
endfor
DEVICE, /CLOSE
SET_PLOT, mydevice
end
