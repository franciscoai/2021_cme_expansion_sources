;+
;*************************************************************************
; Francisco A. Iglesias - UTN-FRM/GEAA - franciscoaiglesias@frm.utn.edu.ar
;
; Reads all the .sav files containing the parameters of the fitted CGS model
; ,following Thernisien et. al. 2009, and plots them vs time.
; 
; History:
;	20180411 - First version
;*************************************************************************
;-
pro  tevo_all
@gcs_config

DIR= GCSD_PATH + '/' + $
    ['GCS_20101212','GCS_20101214','GCS_20110317','GCS_20110605','GCS_20130129',$
    	'GCS_20130209','GCS_20130424(1)','GCS_20130608'] ;

SYM=reform(rebin([-1,-2,-4,-5,-6],5,ceil(n_elements(DIR)/5.)),ceil(n_elements(DIR)/5.)*5)
SYM=SYM[0:n_elements(DIR)-1]
LINE=reform(rebin([1,2,3,4],4,ceil(n_elements(DIR)/4.)),ceil(n_elements(DIR)/4.)*4)
LINE=shift(LINE[0:n_elements(DIR)-1],0)

	;reform(rebin([0,1,2,3],4,ceil(n_elements(DIR)/4.)),ceil(n_elements(DIR)/4.)*4)
;LINE=shift(LINE[0:n_elements(DIR)-1],-3)

;loops on all dirs and loads the data of each event in a struture
nd=n_elements(DIR)
mp=15 ; max number of time points per event
ini=replicate(!values.f_nan,MP)
for j=0, nd-1 do begin
 cd={np:0, awl:ini, awd:ini, hgt:ini, rat:ini, han:ini, rot:ini, lon:ini , lat:ini, date:replicate('YYYY-MM-DDTHH:mm:SS.SSS',MP),$
    dawl:ini, dawd:ini}  
 f=file_search(DIR[j]+'/','*.sav')
 nf=n_elements(f)
 ; stores relevant data from all .sav files in d struct
 for i=0, nf-1 do begin
  restore, f[i]
  awl=2.*(sgui.HAN + asin(sgui.RAT)) *!RADEG
  awd=2.*asin(sgui.RAT) *!RADEG
  cd.np=nf ; numper of valid data points
  cd.awl[i]=awl
  cd.awd[i]=awd
  cd.hgt[i]=sgui.HGT 
  cd.rat[i]=sgui.RAT 
  cd.han[i]=sgui.HAN *!RADEG 
  cd.rot[i]=sgui.ROT*!RADEG 
  cd.lon[i]=sgui.LON  
  cd.lat[i]=sgui.LAT
  cd.date[i]=sgui.ERUPTIONDATE
 endfor

  ;sorts according to date
  ind=sort(cd.date)
  cd.awl=cd.awl[ind]
  cd.awd=cd.awd[ind]
  cd.hgt=cd.hgt[ind]
  cd.rat=cd.rat[ind]
  cd.han=cd.han[ind]
  cd.rot=cd.rot[ind]
  cd.lon=cd.lon[ind] 
  cd.lat=cd.lat[ind]
  cd.date=cd.date[ind]

 ;computes some parameters
  time=str2utc(cd.date[0:nf-1])
  time=float(time.time)/1000./3600.
  cd.dawl[0:nf-1]=deriv(time, cd.awl[0:nf-1]) ; time derivatives of the ang
  cd.dawd[0:nf-1]=deriv(time, cd.awd[0:nf-1])
  
  ;cd.dawl[0:nf-1]=ts_diff(cd.awl[0:nf-1],1)/ts_diff(time,1) ; time deriv from  differences of the ang
  ;cd.dawd[0:nf-1]=ts_diff(cd.awd[0:nf-1],1)/ts_diff(time,1)

  if j eq 0 then d=cd else d=[d,cd] 
endfor

 ;delte last points
 ;d=d[0:4]
 ;if j eq 0 then alld=d else alld=[[alld],[d]] 

; norm. diff of the deriv 
awdiff=d.dawl-d.dawd

;***********Plots derv diff
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=file_dirname(DIR[0]) + '/aw_vel_diff.ps', /LANDSCAPE
!P.Multi = 0
FONT=4.5
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL

for j=0, nd-1 do begin
 if j eq 0 then begin
   plot, d[j].hgt[0:d[j].np-1], awdiff[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j],title= 'Difference of the expansion velocity (AWL-AWD)', $
      xtitle='Height [Rs]', ytitle='Diff. [Deg/h]', yrange=[min(awdiff,/nan),max(awdiff,/nan)],$
      font=FONT, YMARGIN=MARGIN, charthick=1.6, thick=1.2, xrange=[min(d[*].hgt[*],/nan),max(d[*].hgt[*],/nan)];,/ylog
 endif else begin
  oplot, d[j].hgt[0:d[j].np-1], awdiff[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j], thick=1.2;,/ylog
 endelse
endfor
;adds the exp fit 
legend,file_basename(DIR),psym=SYM, linestyle=LINE, /right,charsize=0.95
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;*********Plots norm derv diff
for j=0,nd-1 do awdiff[*,j]/=d[j].dawl[0]/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=file_dirname(DIR[0]) + '/aw_vel_diff_norm.ps', /LANDSCAPE
!P.Multi = 0
FONT=4.5
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL

for j=0, nd-1 do begin
 if j eq 0 then begin
   plot, d[j].hgt[0:d[j].np-1], awdiff[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j],title= 'Difference of the expansion velocity (AWL-AWD)', $
      xtitle='Height [Rs]', ytitle='Norm. Diff. [%]', yrange=[min(awdiff,/nan),max(awdiff,/nan)],$
      font=FONT, YMARGIN=MARGIN, charthick=1.6, thick=1.2, xrange=[min(d[*].hgt[*],/nan),max(d[*].hgt[*],/nan)];,/ylog
 endif else begin
  oplot, d[j].hgt[0:d[j].np-1], awdiff[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j], thick=1.2;,/ylog
 endelse
endfor
legend,file_basename(DIR),psym=SYM, linestyle=LINE, /right,charsize=0.95
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice


;***********Plots all AWL
toplot=d.awl
for j=0,nd-1 do toplot[*,j]/=mean(toplot[d[j].np-3:d[j].np-1,j])/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=file_dirname(DIR[0]) + '/awl_all.ps', /LANDSCAPE
!P.Multi = 0
FONT=4.5
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL

for j=0, nd-1 do begin
 ;d[j].np=4 ; use only first d[j].np time points 
 if j eq 0 then begin
   plot, d[j].hgt[0:d[j].np-1], toplot[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j],title= 'Norm. lateral angular width', $
      xtitle='Height [Rs]', ytitle='AWL/AWL_final [%]', yrange=[min(toplot[0:d[j].np-1,*],/nan),max(toplot[0:d[j].np-1,*],/nan)],$
      font=FONT, YMARGIN=MARGIN, charthick=1.6, thick=1.2, xrange=[min(d[*].hgt[*],/nan),max(d[*].hgt[*],/nan)];,/ylog
 endif else begin
  oplot, d[j].hgt[0:d[j].np-1], toplot[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j], thick=1.2;,/ylog
 endelse
endfor
legend,file_basename(DIR),psym=SYM, linestyle=LINE, position=[6,60],charsize=0.95
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice


;***********Plots all hgt vs time
toplot=d.hgt
;for j=0,nd-1 do toplot[*,j]/=mean(toplot[d[j].np-3:d[j].np-1,j])/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=file_dirname(DIR[0]) + '/hgt_vs_time.ps', /LANDSCAPE
!P.Multi = 0
FONT=4.5
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL

for j=0, nd-1 do begin
 ;d[j].np=4 ; use only first d[j].np time points 
 time=str2utc(d[j].date[0:d[j].np-1])
 time=float(time.time)/1000./3600.
 if file_basename(DIR[j]) eq 'GCS_20130608' then time[6]+=24.
 time-=time[0]
 if j eq 0 then begin
   plot, time, toplot[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j],title= 'Height vs time for all events', $
      xtitle='time from start [h]', ytitle='Height [Rs]', yrange=[min(d[*].hgt[*],/nan),max(d[*].hgt[*],/nan)],$
      font=FONT, YMARGIN=MARGIN, charthick=1.6, thick=1.2, xrange=[0,7];,/ylog
 endif else begin
  oplot, time, toplot[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j], thick=1.2;,/ylog
 endelse
endfor
legend,file_basename(DIR),psym=SYM, linestyle=LINE,charsize=0.95,position=[5,3.5]
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;***********Plots all derv_AWD vs derv_AWL
toplot=d.awd
for j=0,nd-1 do toplot[*,j]/=mean(toplot[d[j].np-3:d[j].np-1,j])/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=file_dirname(DIR[0]) + '/awd_all.ps', /LANDSCAPE
!P.Multi = 0
FONT=4.5
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL

for j=0, nd-1 do begin
 ;d[j].np=4 ; use only first d[j].np time points 
 if j eq 0 then begin
   plot, d[j].hgt[0:d[j].np-1], toplot[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j],title= 'Norm. axial angular width', $
      xtitle='Height [Rs]', ytitle='AWD/AWD_final [%]', yrange=[min(toplot[0:d[j].np-1,*],/nan),max(toplot[0:d[j].np-1,*],/nan)],$
      font=FONT, YMARGIN=MARGIN, charthick=1.6, thick=1.2, xrange=[min(d[*].hgt[*],/nan),max(d[*].hgt[*],/nan)];,/ylog
 endif else begin
  oplot, d[j].hgt[0:d[j].np-1], toplot[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j], thick=1.2;,/ylog
 endelse
endfor
legend,file_basename(DIR),psym=SYM, linestyle=LINE, position=[6,60],charsize=0.95
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;***********Plots
toploty=d.dawl
toplotx=d.dawd
x=0
y=0
for j=0,nd-1 do begin
  x=[x,toplotx[*,j]] ;0:max(where(d[j].hgt lt 10))
  y=[y,toploty[*,j]] ;0:max(where(d[j].hgt lt 10))
endfor
x=x[1:*]
y=y[1:*]
x=x[where(finite(x))]
y=y[where(finite(y))]
fit=linfit(x,y,CHISQR=chisqr)
chisqr=sqrt(chisqr/n_elements(x))

mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=file_dirname(DIR[0]) + '/deriv_awd_vs_deriv_awl.ps', /LANDSCAPE
!P.Multi = 0
FONT=4.5
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL

for j=0, nd-1 do begin
 ;d[j].np=4 ; use only first d[j].np time points 
 if j eq 0 then begin
   plot, toplotx[0:d[j].np-1,j], toploty[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j],title= 'Angular speeds for all events', $
      xtitle='deriv. AWD [Deg/h]', ytitle='deriv. AWL [Deg/h]', yrange=[min(toploty[0:d[j].np-1,*],/nan),1.2*max(toploty[0:d[j].np-1,*],/nan)],$
      font=FONT, YMARGIN=MARGIN, charthick=1.6, thick=1.2, xrange=[min(toplotx,/nan),max(toplotx,/nan)];,/ylog
 endif else begin
  oplot, toplotx[0:d[j].np-1,j], toploty[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j], thick=1.2;,/ylog
 endelse
endfor
xx=indgen(50)*(max(toplotx,/nan)-min(toplotx,/nan))/50+min(toplotx,/nan)
oplot, xx,fit[0]+fit[1]*xx, thick=7;,/ylog
;SYM=-SYM
legend,[file_basename(DIR),'FIT m = '+strtrim(fit[1],2)+'; rms_err='+strtrim(chisqr,2)],psym=[SYM,0], linestyle=[LINE,0], /left,charsize=0.95
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;***********Plots
toploty=d.awl/d.awd
toplotx=d.hgt
;for j=0,nd-1 do toploty[*,j]/=mean(toploty[d[j].np-3:d[j].np-1,j])/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=file_dirname(DIR[0]) + '/awl_over_awd_vs_height.ps', /LANDSCAPE
!P.Multi = 0
FONT=4.5
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL

for j=0, nd-1 do begin
 ;d[j].np=4 ; use only first d[j].np time points 
 if j eq 0 then begin
   plot, toplotx[0:d[j].np-1,j], toploty[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j],title= 'Angular width for all events', $
     xtitle='Height [Rs]', ytitle='AWL/AWD', yrange=[min(toploty[0:d[j].np-1,*],/nan),max(toploty[0:d[j].np-1,*],/nan)],$
      font=FONT, YMARGIN=MARGIN, charthick=1.6, thick=1.2, xrange=[min(toplotx,/nan),max(toplotx,/nan)];,/ylog
 endif else begin
  oplot, toplotx[0:d[j].np-1,j], toploty[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j], thick=1.2;,/ylog
 endelse
endfor
;SYM=-SYM
legend,file_basename(DIR),psym=SYM, linestyle=LINE,charsize=0.95
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;***********Plots
toploty=d.dawd
toplotx=d.hgt
for j=0,nd-1 do toploty[*,j]/=mean(d[j].awd[d[j].np-3:d[j].np-1])/100. ; nor to final value
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=file_dirname(DIR[0]) + '/deriv_awd_vs_height_norm.ps', /LANDSCAPE
!P.Multi = 0
FONT=4.5
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL

for j=0, nd-1 do begin
 ;d[j].np=4 ; use only first d[j].np time points 
 if j eq 0 then begin
   plot, toplotx[0:d[j].np-1,j], toploty[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j],title= 'Axial angular speed for all events', $
      xtitle='Height [Rs]', ytitle='Deriv. AWD/AWD_final [%/h]', yrange=[min(toploty[0:d[j].np-1,*],/nan),max(toploty[0:d[j].np-1,*],/nan)],$
      font=FONT, YMARGIN=MARGIN, charthick=1.6, thick=1.2, xrange=[min(toplotx,/nan),max(toplotx,/nan)];,/ylog
 endif else begin
  oplot, toplotx[0:d[j].np-1,j], toploty[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j], thick=1.2;,/ylog
 endelse
endfor
oplot, xx,1.+fltarr(n_elements(xx)), thick=5;,/ylog
;SYM=-SYM
legend,[file_basename(DIR),'1%'],psym=[SYM,0], linestyle=[LINE,1], /right,charsize=0.95
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;***********Plots
toploty=d.dawl
toplotx=d.hgt
for j=0,nd-1 do toploty[*,j]/=mean(d[j].awl[d[j].np-3:d[j].np-1])/100. ; norm to final value
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=file_dirname(DIR[0]) + '/deriv_awl_vs_height_norm.ps', /LANDSCAPE
!P.Multi = 0
FONT=4.5
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL

for j=0, nd-1 do begin
 ;d[j].np=4 ; use only first d[j].np time points 
 if j eq 0 then begin
   plot, toplotx[0:d[j].np-1,j], toploty[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j],title= 'Longitudinal angular speed for all events', $
      xtitle='Height [Rs]', ytitle='Deriv. AWL/AWL_final [%/h]', yrange=[min(toploty[0:d[j].np-1,*],/nan),max(toploty[0:d[j].np-1,*],/nan)],$
      font=FONT, YMARGIN=MARGIN, charthick=1.6, thick=1.2, xrange=[min(toplotx,/nan),max(toplotx,/nan)];,/ylog
 endif else begin
  oplot, toplotx[0:d[j].np-1,j], toploty[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j], thick=1.2;,/ylog
 endelse
endfor
oplot, xx,1.+fltarr(n_elements(xx)), thick=5;,/ylog
;SYM=-SYM
legend,[file_basename(DIR),'1%'],psym=[SYM,0], linestyle=[LINE,1], /right,charsize=0.95
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice
end


