
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
pro  oldtevo_FMwLASCO20110605
@gcs_config


DIR=GCSD_PATH+'/GCS_20110605'

f=file_search(DIR+'/','*.sav')
nf=n_elements(f)
; stores relevant data from all .sav files in d struct
restore, f[0]
awl=2.*(sgui.HAN + asin(sgui.RAT)) *!RADEG
awd=2.*asin(sgui.RAT)  *!RADEG

d={awl:awl, awd:awd, hgt:sgui.HGT, rat:sgui.RAT, han:sgui.HAN *!RADEG, rot:sgui.ROT, lon:sgui.LON , lat:sgui.LAT,date:sgui.ERUPTIONDATE}
for i=1, nf-1 do begin
  restore, f[i]
  awl=2.*(sgui.HAN + asin(sgui.RAT)) *!RADEG
  awd=2.*asin(sgui.RAT) *!RADEG
  d=[d, {awl:awl, awd:awd, hgt:sgui.HGT, rat:sgui.RAT, han:sgui.HAN *!RADEG, rot:sgui.ROT, lon:sgui.LON , lat:sgui.LAT,date:sgui.ERUPTIONDATE}]
endfor
d=d[sort(d.date)] ; sorts by date

;plots settings
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=DIR + '/tevo.ps', /LANDSCAPE
!P.Multi = [0,3,2,0,1]
FONT=4.5
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL

;Plots
time=str2utc(d.date)
time=float(time.time)/1000./3600.
;window, 2, xsize=1000, ysize=800
plot, time,  d.hgt, ylog=YLOG, psym=-2, ytitle='HEIGHT [Ro]', font=FONT, YMARGIN=MARGIN, $
      xtitle='time [h]', yrange=[min(d.hgt)*0.95,max(d.hgt)*1.05]
plot, time, d.han, ylog=YLOG, psym=-2, ytitle='HALF ANGLE [Deg]', font=FONT, YMARGIN=MARGIN, $
      xtitle='time [h]', yrange=[min(d.han)*0.95,max(d.han)*1.05]
plot, time, d.rat, ylog=YLOG, psym=-2, ytitle='RATIO', font=FONT, YMARGIN=MARGIN, $
      xtitle='time [h]', yrange=[min(d.rat)*0.95,max(d.rat)*1.05]
plot, time, d.rot, ylog=YLOG, psym=-2, ytitle='ROT [DEG?]', font=FONT, YMARGIN=MARGIN, $
      xtitle='time [h]', yrange=[min(d.rot)*0.95, max(d.rot)*1.05]
plot, time, d.awl, ylog=YLOG, psym=-2, ytitle='AWL [Deg]', font=FONT, YMARGIN=MARGIN, $
      xtitle='time [h]', yrange=[min(d.awl)*0.95,max(d.awl)*1.05]
plot, time, d.awd, psym=-2, ytitle='AWD [Deg]', font=FONT, YMARGIN=MARGIN, $
      xtitle='time [h]', yrange=[min(d.awd)*0.95,max(d.awd)*1.05]

DEVICE, FILENAME=DIR + '/tevo_der.ps', /LANDSCAPE
;window,3, xsize=900, ysize=500
!P.Multi = 0
plot, time, deriv(time,d.awl), psym=-4, linestyle=0,title= 'Derivatives of the Angular Widths', $
      xtitle='time [h]', ytitle='Derivative [Deg/h]', yrange=[-max(deriv(time,d.awl))*0.05,$
      max(deriv(time,d.awl))*1.05]
oplot, time, deriv(time,d.awd), psym=-2
legend,[' AWL  ',' AWD  '],psym=[4,2], /left
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

end
