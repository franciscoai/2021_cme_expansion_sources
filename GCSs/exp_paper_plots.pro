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

;Function to fit the AW's vs height profiles
function paper_plots_func2fit, x, p
  ;return, p[0]*(x-1.)^(p[1]) ; power law
  ;return, p[0]*(x-1.)^(p[1])*2.7182818284590452^(-p[2]*x) ; power law with exp cutoff
  ;return, p[0]*(1.-(p[1]/x)^p[2]);+p[2]*x ; what Bein_2013 uses for the mass vs height  !!!!!!
  return, p[0]+p[1]*x^p[2];+p[2]*x ; power law plus offset
  ;return, p[1]*x^p[2]-p[1]*1.^p[2] ; power law plus offset given by the initial AW (form WINI and HINI)
end 


;Derivative of paper_plots_func2fit
function paper_plots_func2fit_deriv, x, p
  return, p[1]*p[2]*x^(p[2]-1.);+p[2]*x 
end 

;Derivative of paper_plots_func2fit normalized
function paper_plots_func2fit_deriv_norm, x, p
  return, p[1]*p[2]*x^(p[2]-1.) / (p[0]+p[1]*x^p[2]);+p[2]*x 
end 

;Function to fit the height vs time profiles
function paper_plots_func2fit2  , x, p
    return, p[0]*x^2. + p[1]*x + p[2] ; regular cuadratic
end 

;Derivative of paper_plots_func2fit2,4 and 5
function paper_plots_func2fit2_deriv, x, p
  return, 2.*p[0]*x + p[1]
end 

; ;Function to fit the  height vs timeprofiles
; ;quadratic until the min, then constant
; function paper_plots_func2fit3, x, p
;   xmin = -0.5 * p[1]/p[0]
;   xl=where(x lt xmin,nxl)
;   xh=where(x ge xmin,nxh)

;   if nxl eq 0 and nxh gt 0 then $ ; all x > xmin
;     return, p[0]*x[xh]^2. + p[1]*x[xh] + p[2]

;   xl=fltarr(nxl)+xmin
;   if nxl gt 0 and nxh eq 0 then $  ; all x < xmin
;     return, p[0]*xl^2. + p[1]*xl + p[2] ; some x > xmin and some x < xmin

;   return, [p[0]*xl^2. + p[1]*xl + p[2],p[0]*x[xh]^2. + p[1]*x[xh] + p[2]]
; end 

;Function to fit the height [Rs] vs time [Jul] profiles
function paper_plots_func2fit4  , x, p
common paper_plots_func2fit_block, HINI, WINI
    return, p[0]*x^2. + p[1]*x + HINI - p[0]*P[2]^2. - p[1]*P[2] ; cuadratic plus forcing touching HINI at P[2]
end 

;Derivative wrt x of paper_plots_func2fit4
function deriv_paper_plots_func2fit4  , x, p
common paper_plots_func2fit_block, HINI, WINI
    return, 2.*p[0]*x+ p[1] ; cuadratic plus forcing touching HINI at P[2]
end 

;second derivative wrt x of paper_plots_func2fit4
function scnd_deriv_paper_plots_func2fit4  ,p
common paper_plots_func2fit_block, HINI, WINI
    return, 2.*p[0]; cuadratic plus forcing touching HINI at P[2]
end 

;largest root of paper_plots_func2fit4=y
function root_paper_plots_func2fit4,  y, p
  common paper_plots_func2fit_block, HINI, WINI
  bla=HINI - p[0]*P[2]^2. - p[1]*P[2] - y
  r1=(-p[1]+sqrt(p[1]^2.-4.*p[0]*bla))/(2*p[0])
  r2=(-p[1]-sqrt(p[1]^2.-4.*p[0]*bla))/(2*p[0])
  return, max(r1,r1) ; cuadratic plus forcing touching HINI at P[2]
end 

;Function to fit the widths vs time profiles
function paper_plots_func2fit5  , x, p
  common paper_plots_func2fit_block, HINI, WINI
    return, p[0]*x^2. + p[1]*x + WINI; cuadratic plus forcing touching HINI at x=0
end 

; To use in paper_plots_cgsWFO. Returns the negative of the BP_X vector 
; of the GCS. See EQ. 27 in Thernisien et. al. 2011
function paper_plots_cgsBPx, beta
  common paper_plots_cgsWFO_par, k, b, rho
  x=(rho + b * k^2. * sin(beta)) / (1. - k^2.)
  r=sqrt((x^2. + (b^2. * k^2. - rho^2.)) / (1. - k^2.)) 
  bp=(x + r)*cos(beta)
  ;print, beta, -bp
  ; stop
  return, -bp
end

; Computes the maximum Face-on (Lateral) width of the GCS based on the 3
; relevant  model parameters. This corresponds to the  maximum extension  of the CME
; in a line that is perpendicular to the propagation direction and lies in the lateral plane
; See EQ. 27 in Thernisien et. al. 2011
function  paper_plots_cgsWFO, alpha, delta, height
  common paper_plots_cgsWFO_par, k, b, rho
  k=sin(delta)
  b=height / cos(alpha)
  rho=height * tan(alpha)
  ; BLA=indgen(10)*0.01*3.1415/2.
  ; plot, BLA, paper_plots_cgsBPx(BLA)
  minb=AMOEBA(1.0e-5, function_name='paper_plots_cgsBPx', SCALE=3.16/2., P0 = [0.],FUNCTION_VALUE=fval)
  IF minb EQ -1 THEN MESSAGE, 'AMOEBA failed to converge'
  return, -fval[0]
end

; Computes the maximum Edge-on (axial) width of the GCS based on the 3
; relevant  model parameters. This corresponds to the  maximum extension  of the CME
; in a line that is perpendicular to the propagation direction and lies in the axial plane
; See EQ. 29 in Thernisien et. al. 2011
function  paper_plots_cgsWEO, alpha, delta, height
  common paper_plots_cgsWFO_par, k, b, rho
  k=sin(delta)
  b=height / cos(alpha)
  rho=height * tan(alpha)
  return, k*(b + rho) / (1. - k^2.)
end

;Derivative of paper_plots_func2fit normalized
function paper_plots_func2fit_deriv_norm_root, x
  common paper_plots_func2fit_block2, SETAW,pr
  return, (pr[1]*pr[2]*x^(pr[2]-1.) / (pr[0]+pr[1]*x^pr[2]))-SETAW;+p[2]*x 
end 


;Computes the Settling height of the angular widths, i.e. to reach AW=SETAW*final_AW 
function paper_plots_setHeight , p
  common paper_plots_func2fit_block2, SETAW,pr
    pr=p
    ;return, ((SETAW-1.)*p[0]/p[1])^(1./p[2])
    return, zbrent(2,30,func_name='paper_plots_func2fit_deriv_norm_root')
end 

pro  exp_paper_plots
@gcs_config
common paper_plots_func2fit_block, HINI, WINI
common paper_plots_func2fit_block2, SETAW,pr

DIR= $
    [ '20101212','20101214',$  ;not polar ,'20110215','20110213'
      '20110317','20110605','20130123',$
      '20130129','20130209','20130424',$
      '20130502','20130517','20130527',$
      '20130608'];,'20110214'] MAYBEEE
DIR1=GCSD_PATH +'/'+ DIR
DIR= GCSD_PATH + '/GCS_' + DIR
nd=n_elements(DIR)
OPATH=file_dirname(DIR[0])+'/plots' ; output path
AWRANGE=[15,145] ; AW range to plot
HRANGE=[1,8] ; height range to plot in Rs
VELRANGE=[0,1500] ; vel range to plot in km/s
ACCRANGE=[0,50] ; vel range to plot in km/s^2
func2fitINI=[90.,-90,-1.] ; initial conditions when fitting paper_plots_func2fit
func2fitINI4=[1.,1.,-10] ; initial conditions when fitting paper_plots_func2fit4
func2fitINI5=[1.,1.,-30] ; initial conditions when fitting paper_plots_func2fit5
HINI=1.25 ; Initial height of the CMEs, i.e. of their precursive quiscent cavity, see Gibson et. al. 2006
WINI=0.25 ; initial width of the CMEs, basically we take them as circles
SETAW=0.05 ;0.50 ; Used to compute the Settling height of the angular widths.
ERR=0.15 ; Max error in the fitted alpha (half angle) and kappa (aspect ratio)

SYM=reform(rebin([-1,-2,-4,-5,-6],5,ceil(nd/5.)),ceil(nd/5.)*5)
SYM=SYM[0:nd-1]
LINE=reform(rebin([1,2,3,4],4,ceil(nd/4.)),ceil(nd/4.)*4)
LINE=0*shift(LINE[0:nd-1],0) + 2
COLORS=fix(256./nd)*(indgen(nd))
COLORT=6 ; color table
RS2KM= 695700. ; Rs to km
XXMARGIN=[20,1]
FONTSZ=14 ; font size

;loops on all dirs and loads the data of each event in a struture
mp=15 ; max number of time points per event
ini=replicate(!values.f_nan,MP)
for j=0, nd-1 do begin
 cd={np:0, awl:ini, awd:ini, hgt:ini, dhgt:ini, rat:ini, han:ini, rot:ini, lon:ini , lat:ini, date:replicate('YYYY-MM-DDTHH:mm:SS.SSS',MP),$
    dawl:ini, dawd:ini, wl:ini, wd:ini, dwl:ini, dwd:ini, eawl_l:ini, eawl_h:ini, eawd_l:ini, eawd_h:ini}  
 f=file_search(DIR[j]+'/','*.sav')
 nf=n_elements(f)
 ; stores relevant data from all .sav files in d struct
 for i=0, nf-1 do begin
  restore, f[i]
  awl=2.*(sgui.HAN + asin(sgui.RAT)) *!RADEG ; lateral ang width
  awd=2.*asin(sgui.RAT) *!RADEG ; axial ang width
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
  cd.wl[i]=paper_plots_cgsWFO(sgui.HAN,sgui.RAT,sgui.HGT) ; lateral width in Rs
  cd.wd[i]=paper_plots_cgsWEO(sgui.HAN,sgui.RAT,sgui.HGT) ; axial width in Rs
  ; error estimation
  awl=2.*(sgui.HAN*(1+ERR) + asin(sgui.RAT*(1+ERR))) *!RADEG ; max lateral ang width
  awd=2.*asin(sgui.RAT*(1+ERR)) *!RADEG ; max axial ang width
  cd.eawl_h[i]=awl
  cd.eawd_h[i]=awd
  awl=2.*(sgui.HAN*(1-ERR) + asin(sgui.RAT*(1-ERR))) *!RADEG ; min lateral ang width
  awd=2.*asin(sgui.RAT*(1-ERR))*!RADEG ; min axial ang width
  cd.eawl_l[i]=awl
  cd.eawd_l[i]=awd  
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
  cd.wl=cd.wl[ind]
  cd.wd=cd.wd[ind]
  cd.eawl_h=cd.eawl_h[ind]
  cd.eawd_h=cd.eawd_h[ind]
  cd.eawl_l=cd.eawl_l[ind]
  cd.eawd_l=cd.eawd_l[ind]

 ;computes some parameters
  time=str2utc(cd.date[0:nf-1])
  time=double(time.time)/1000. + (time.mjd-time[0].mjd)*24.*3600. ; date to s
  cd.dawl[0:nf-1]=deriv(time, cd.awl[0:nf-1]) ; time derivatives of the ang widths
  cd.dawd[0:nf-1]=deriv(time, cd.awd[0:nf-1]) 
  ;cd.dawl[0:nf-1]=ts_diff(cd.awl[0:nf-1],1)/ts_diff(time,1) ; time deriv from  differences of the ang
  ;cd.dawd[0:nf-1]=ts_diff(cd.awd[0:nf-1],1)/ts_diff(time,1)

  if j eq nd-1 then begin ; for last case uses diff
    cd.dhgt[0:nf-1]=ts_diff(cd.hgt[0:nf-1],1)/ts_diff(time,1)
    cd.dwl[0:nf-1]=ts_diff(cd.wl[0:nf-1],1)/ts_diff(time,1)  
    cd.dwd[0:nf-1]=ts_diff(cd.wd[0:nf-1],1)/ts_diff(time,1)
  endif else begin 
    cd.dhgt[0:nf-1]=deriv(time, cd.hgt[0:nf-1])     
    cd.dwl[0:nf-1]=deriv(time, cd.wl[0:nf-1]) ; time derivatives of the widths
    cd.dwd[0:nf-1]=deriv(time, cd.wd[0:nf-1])    
  endelse

  if j eq 0 then d=cd else d=[d,cd] 
endfor

 ;delte last points
 ;d=d[0:4]
 ;if j eq 0 then alld=d else alld=[[alld],[d]] 

; norm. diff of the deriv 
awdiff=d.dawl-d.dawd

; creates OPATH
dummy=file_search(OPATH,count=count)
if count eq 0 then file_mkdir, OPATH   ; creates the output directory

;******00***Plots all hgt vs time
toplot=d.hgt
xx = indgen(80)*0.1 
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/height_vs_time.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

err=toplot-toplot+0.1
err[where(toplot lt 2.5)] = 1.55/mean([15., 23.])*0.1
for j=0, nd-1 do begin
 ;d[j].np=4 ; use only first d[j].np time points 
 time=str2utc(d[j].date[0:d[j].np-1])
 time=float(time.time)/1000./3600. + (time.mjd-time[0].mjd)*24. ; date to h of day
 ;if file_basename(DIR[j]) eq 'GCS_20130608' then time[6]+=24.
 ;time-=time[0]
 pi=replicate({limited:[0,0], limits:[0.D,0.D], tied:''},3)    ;creates structure for PARINFO
 ;pi[2].tied = strtrim(HINI,2)+'+ 0.25 *  P[1]^2. / P[0]'; eq. constraint to force the parabloa min to be at H0
 pi[2].limited[1]=1; enables upper limit
 pi[2].limits[1]=time[0]; upper limit
 func2fitINI4[2]=time[0] ; initial condition of starting time close to first meas
 pi[0].limited[0]=1 
 pi[0].limits[0]=0

 if j eq 0 then begin
  fit_ht=mpfitfun('paper_plots_func2fit4',time,toplot[0:d[j].np-1,j],err[0:d[j].np-1,j],func2fitINI4, PARINFO=pi)
  tini=fit_ht[2]
  time-=tini
  plot, time, toplot[0:d[j].np-1,j], psym=-SYM[j], linestyle=LINE[j],title= 'Height vs time for all events', $
      xtitle='time from start [h]', ytitle='Height [Rs]', yrange=[min(d[*].hgt[*],/nan),max(d[*].hgt[*],/nan)],$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=[0,8];,/ylog 
 endif else begin
  fit_ht=[fit_ht,mpfitfun('paper_plots_func2fit4',time,toplot[0:d[j].np-1,j],err[0:d[j].np-1,j],func2fitINI4, PARINFO=pi)]  
  tini=[tini,fit_ht[3*j+2]]; initial time, ie.e. height = HINI  
  time-=tini[j]
  oplot, time, toplot[0:d[j].np-1,j], psym=-SYM[j], linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
 endelse
  xx1 = xx+fit_ht[3*j+2];+time[0];
  oplot, xx, paper_plots_func2fit4(xx1,[fit_ht[3*j:3*j+2]]),linestyle=LINE[j],color=COLORS[j] 
endfor
oplot, [-1,6], replicate(HINI, 2),linestyle=1,color=COLORS[0] 
legend,file_basename(DIR1) ,psym=SYM, linestyle=LINE,charsize=0.8, color=COLORS,position=[6.2,4.0],$
  spacing=1, pspacing=2
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice


;***********Plots derv diff
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/dawd-dawl.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

for j=0, nd-1 do begin
 if j eq 0 then begin
   plot, d[j].hgt[0:d[j].np-1], awdiff[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j],title= 'Difference of the expansion velocity (AWL-AWD)', $
      xtitle='Height [Rs]', ytitle='Diff. [deg/s]', yrange=[min(awdiff,/nan),max(awdiff,/nan)],$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=[min(d[*].hgt[*],/nan),max(d[*].hgt[*],/nan)];,/ylog
 endif else begin
  oplot, d[j].hgt[0:d[j].np-1], awdiff[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
 endelse
endfor
;adds the exp fit 
legend,file_basename(DIR1),psym=SYM, linestyle=LINE, /right,charsize=0.8,color=COLORS
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

; ;*********Plots norm derv diff
; for j=0,nd-1 do awdiff[*,j]/=d[j].dawl[0]/100.
; mydevice = !D.NAME
; SET_PLOT, 'PS'
; DEVICE, FILENAME=OPATH + '/aw_vel_diff_norm.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
; !P.Multi = 0
; FONT=FONTSZ
; MARGIN=[5,5]
; YLOG=0
; !P.Color = '000000'xL
; !P.Background = 'FFFFFF'xL
; loadct, COLORT

; for j=0, nd-1 do begin
;  if j eq 0 then begin
;    plot, d[j].hgt[0:d[j].np-1], awdiff[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j],title= 'Difference of the expansion velocity (AWL-AWD)', $
;       xtitle='Height [Rs]', ytitle='Norm. Diff. [%]', yrange=[min(awdiff,/nan),max(awdiff,/nan)],$
;       font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=[min(d[*].hgt[*],/nan),max(d[*].hgt[*],/nan)];,/ylog
;  endif else begin
;   oplot, d[j].hgt[0:d[j].np-1], awdiff[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
;  endelse
; endfor
; legend,file_basename(DIR1),psym=SYM, linestyle=LINE, /right,charsize=0.8,color=COLORS
; ; Close the PostScript file:
; DEVICE, /CLOSE
; SET_PLOT, mydevice


;***********Plots all AWL
toplot=d.awl
;for j=0,nd-1 do toplot[*,j]/=mean(toplot[d[j].np-3:d[j].np-1,j])/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/awl_vs_height.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

xx = indgen(150)*0.1 
errawl=0
 for j=0, nd-1 do begin
  x=d[j].hgt[0:d[j].np-1]
  y=toplot[0:d[j].np-1,j]
 ; d[j].np=4 ; use only first d[j].np time points 
 if j eq 0 then begin
  ; fit_awlhs a power law (linear in log space)
   fit_awlh=mpfitfun('paper_plots_func2fit',x,y,5.+(y-y),func2fitINI)
   ;chi=sqrt(chi/n_elements(x))  
   plot, x, y, psym=-SYM[j], linestyle=LINE[j],title= 'Lateral angular width', $
      xtitle='Height [Rs]', ytitle='AWL', yrange=AWRANGE,$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=HRANGE;,/ylog
 endif else begin
  ;fits a power law
  fit_awlh=[fit_awlh,mpfitfun('paper_plots_func2fit',x,y,5.+(y-y),func2fitINI)]
  ;chi=[chi,sqrt(chi/n_elements(x))]
  oplot, x, y, psym=-SYM[j], linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
 endelse
 ;errplot, x,d[j].eawl_l[0:d[j].np-1],d[j].eawl_h[0:d[j].np-1], color=COLORS[j], thick=0.5
 oplot, xx+x[0], paper_plots_func2fit(xx+x[0],[fit_awlh[3*j:3*j+2]]),linestyle=LINE[j],color=COLORS[j] ;, psym=SYM[j], linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
 ; computes average relative error in awl
errawl=[errawl,mean((d[j].eawl_h[0:d[j].np-1]-d[j].eawl_l[0:d[j].np-1])),y[d[j].np-1]]
endfor
;print, fit;, chi
legend,file_basename(DIR1),psym=SYM, linestyle=LINE,charsize=0.8, color=COLORS
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;***********Plots all AWD
toplot=d.awd
;for j=0,nd-1 do toplot[*,j]/=mean(toplot[d[j].np-3:d[j].np-1,j])/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/awd_vs_height.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

xx = indgen(150)*0.1
errawd=0
 for j=0, nd-1 do begin
  x=d[j].hgt[0:d[j].np-1]
  y=toplot[0:d[j].np-1,j]
 ; d[j].np=4 ; use only first d[j].np time points 
 if j eq 0 then begin
  ; fit1s a power law (linear in log space)
   fit_awdh=mpfitfun('paper_plots_func2fit',x,y,5.+(y-y),func2fitINI)
   ;chi=sqrt(chi/n_elements(x))  
   ;yrange=[min(toplot[0:d[j].np-1,*],/nan),max(toplot[0:d[j].np-1,*],/nan)]
   ;xrange=[min(d[*].hgt[*],/nan),max(d[*].hgt[*],/nan)]
   plot, x, y, psym=-SYM[j], linestyle=LINE[j],title= 'Axial angular width', $
      xtitle='Height [Rs]', ytitle='AWD', yrange=[0,75],ystyle=1,$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j],$
       xrange=HRANGE;,/ylog
 endif else begin
  fit_awdh=[fit_awdh,mpfitfun('paper_plots_func2fit',x,y,5.+(y-y),func2fitINI)]
  ;chi=[chi,sqrt(chi/n_elements(x))]
  oplot, x, y, psym=-SYM[j], linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
 endelse
 ;errplot, x,d[j].eawd_l[0:d[j].np-1],d[j].eawd_h[0:d[j].np-1], color=COLORS[j], thick=0.5
 oplot, xx+x[0], paper_plots_func2fit(xx+x[0],[fit_awdh[3*j:3*j+2]]),linestyle=LINE[j],color=COLORS[j] ;, psym=SYM[j], linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
 errawd=[errawd,mean((d[j].eawd_h[0:d[j].np-1]-d[j].eawd_l[0:d[j].np-1])),y[d[j].np-1]]
 ;stop ; CHECK THAT THE TIME-DEPENDENT ERROR
endfor
;print, fit_awdh;, chi
;legend,file_basename(DIR1),psym=SYM, linestyle=LINE,charsize=0.8, color=COLORS
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice
;
;
;
;
;
;
;
print,"********************TABLE DATA************************"
TMFT="(I8)"
for j=0, nd-1 do begin
  t6=root_paper_plots_func2fit4(6., fit_ht[3*j:3*j+2]) ; time at 6 Rs
  ;print,fit_ht[3*j:3*j+2]
  ;print, t6
  vel=deriv_paper_plots_func2fit4(t6,fit_ht[3*j:3*j+2]) * RS2KM / 3600. ; vel at 6 Rs in km/s
  acc=scnd_deriv_paper_plots_func2fit4(fit_ht[3*j:3*j+2]) * RS2KM * 1000. / 3600.^2.; mean acc m/s^2.
  date=d[j].date[0] ; date of first point
  lat=d[j].lat[total(finite(d[j].lat[*]))-1]*!RADEG
  lon=d[j].lon[total(finite(d[j].lon[*]))-1]*!RADEG-tim2carr(d[j].date[total(finite(d[j].lat[*]))-1])
  awd=errawd[0:*:2]
  awd=awd[j+1]
  eawd=errawd[1:*:2]
  eawd=eawd[j]
  awl=errawl[0:*:2]
  awl=awl[j+1]
  eawl=errawl[1:*:2]
  eawl=eawl[j] 
  print, strmid(date,0,19)+' & '+string(lat,format=TMFT) $
        +' & '+string(lon,format=TMFT)+' & '+string(vel,format=TMFT) $
        +' & '+string(acc,format=TMFT)+' & '+string(awd,format=TMFT) $
        +'$\pm$'+string(eawd,format=TMFT)+' & '+string(awl,format=TMFT) $
        +'$\pm$'+string(eawl,format=TMFT)+'\\'
endfor
;
;
;
;
;
;

;***********Plots
;for j=0,nd-1 do toploty[*,j]/=mean(toploty[d[j].np-3:d[j].np-1,j])/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/dhawd_vs_height.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

for j=0, nd-1 do begin
toplotx=d[j].hgt[0:d[j].np-1]
toploty=paper_plots_func2fit_deriv_norm(toplotx,[fit_awdh[3*j:3*j+2]])

toplotxx=indgen(200)*0.1+toplotx[0]
toplotyy=paper_plots_func2fit_deriv_norm(toplotxx,[fit_awdh[3*j:3*j+2]])

 if j eq 0 then begin
   plot, toplotx, toploty, psym=-SYM[j], linestyle=LINE[j],title= 'Axial angular width change rate for all events', $
     xtitle='Height [Rs]', ytitle='Norm. change rate [1/Rs]', yrange=[0.001,6],$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=HRANGE,/ylog
   oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
 endif else begin
  oplot, toplotx, toploty, psym=-SYM[j], linestyle=0, thick=1.2, color=COLORS[j];,/ylog
  oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog  
 endelse
endfor
;SYM=-SYM
legend,file_basename(DIR1),psym=SYM, linestyle=LINE,charsize=0.8, color=COLORS
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice


;***********Plots
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/dhawl_minus_dhawd_norm_vs_height.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

for j=0, nd-1 do begin

time=str2utc(d[j].date[0:d[j].np-1])
time=float(time.time)/1000./3600. + (time.mjd-time[0].mjd)*24. ; date to h
toplotx=paper_plots_func2fit4(time,[fit_ht[3*j:3*j+2]])  

;toplotx=d[j].hgt[0:d[j].np-1]
toploty=paper_plots_func2fit_deriv_norm(toplotx,[fit_awlh[3*j:3*j+2]])-$
  paper_plots_func2fit_deriv_norm(toplotx,[fit_awdh[3*j:3*j+2]])

toplotxx=indgen(200)*0.1+toplotx[0]
toplotyy=paper_plots_func2fit_deriv_norm(toplotxx,[fit_awlh[3*j:3*j+2]])-$
  paper_plots_func2fit_deriv_norm(toplotxx,[fit_awdh[3*j:3*j+2]])

 if j eq 0 then begin
   plot, toplotx, toploty, psym=-SYM[j], linestyle=LINE[j], $
     xtitle='Height [Rs]', ytitle='Difference of norm. angular change rate [1/Rs]', yrange=[-0.7,0.7],$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=HRANGE;,/ylog
 oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
  myy=toplotyy
  mxx=toplotxx
 endif else begin
  oplot, toplotx, toploty, psym=-SYM[j], linestyle=0, thick=1.2, color=COLORS[j];,/ylog
  oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog  
  myy=[[myy],[toplotyy]]
  mxx=[[mxx],[toplotxx]]
 endelse
endfor
allx=indgen(100)*0.1
allm=allx
for r=0, n_elements(allx)-1 do begin
  ind=where(mxx[0,*] le allx[r], count)
  if count eq 12 then begin
    blam=interpol(myy[*,ind[0]],mxx[*,ind[0]],allx[r])
    for j=1, count-1 do blam=[blam,interpol(myy[*,ind[j]],mxx[*,ind[j]],allx[r])]
    allm[r]=mean(blam);/=count
  endif else begin
    allm[r]=!values.f_nan
  endelse
endfor
ind=where(~finite(allm))
;allm=smooth(allm,10,/nan)
;allm[ind]=!values.f_nan
oplot, allx, allm, psym=0, linestyle=0, thick=9, color=0;,/ylog 
legend,file_basename(DIR1),psym=SYM, linestyle=LINE,charsize=0.8, color=COLORS
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;***********Plots
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/dhawl_minus_dhawd_vs_height.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

for j=0, nd-1 do begin

time=str2utc(d[j].date[0:d[j].np-1])
time=float(time.time)/1000./3600. + (time.mjd-time[0].mjd)*24. ; date to h
toplotx=paper_plots_func2fit4(time,[fit_ht[3*j:3*j+2]])  

;toplotx=d[j].hgt[0:d[j].np-1]
toploty=paper_plots_func2fit_deriv(toplotx,[fit_awlh[3*j:3*j+2]])-$
  paper_plots_func2fit_deriv(toplotx,[fit_awdh[3*j:3*j+2]])

toplotxx=indgen(200)*0.1+toplotx[0]
toplotyy=paper_plots_func2fit_deriv(toplotxx,[fit_awlh[3*j:3*j+2]])-$
  paper_plots_func2fit_deriv(toplotxx,[fit_awdh[3*j:3*j+2]])

 if j eq 0 then begin
   plot, toplotx, toploty, psym=-SYM[j], linestyle=LINE[j], $
     xtitle='Height [Rs]', ytitle='Difference of angular change rate [deg/Rs]',ystyle=1, yrange=[-1,40],$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=HRANGE;,/ylog
 oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
  myy=toplotyy
  mxx=toplotxx
 endif else begin
  oplot, toplotx, toploty, psym=-SYM[j], linestyle=0, thick=1.2, color=COLORS[j];,/ylog
  oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog  
  myy=[[myy],[toplotyy]]
  mxx=[[mxx],[toplotxx]]
 endelse
endfor
allx=indgen(100)*0.1
allm=allx
for r=0, n_elements(allx)-1 do begin
  ind=where(mxx[0,*] le allx[r], count)
  if count eq 12 then begin
    blam=interpol(myy[*,ind[0]],mxx[*,ind[0]],allx[r])
    for j=1, count-1 do blam=[blam,interpol(myy[*,ind[j]],mxx[*,ind[j]],allx[r])]
    allm[r]=mean(blam);/=count
  endif else begin
    allm[r]=!values.f_nan
  endelse
endfor
ind=where(~finite(allm))
;allm=smooth(allm,10,/nan)
;allm[ind]=!values.f_nan
oplot, allx, allm, psym=0, linestyle=0, thick=9, color=0;,/ylog 
legend,file_basename(DIR1),psym=SYM, linestyle=LINE,charsize=0.8, color=COLORS,position=[5,37]
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;***********Plots
;for j=0,nd-1 do toploty[*,j]/=mean(toploty[d[j].np-3:d[j].np-1,j])/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/dhawl_vs_height.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

for j=0, nd-1 do begin
toplotx=d[j].hgt[0:d[j].np-1]
toploty=paper_plots_func2fit_deriv_norm(toplotx,[fit_awlh[3*j:3*j+2]])

toplotxx=indgen(200)*0.1+toplotx[0]
toplotyy=paper_plots_func2fit_deriv_norm(toplotxx,[fit_awlh[3*j:3*j+2]])

 if j eq 0 then begin
   plot, toplotx, toploty, psym=-SYM[j], linestyle=LINE[j],title= 'Lateral angular width change rate for all events', $
     xtitle='Height [Rs]', ytitle='Norm. change rate [1/Rs]', yrange=[0.001,6],$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=HRANGE ,/ylog
   oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
 endif else begin
  oplot, toplotx, toploty, psym=-SYM[j], linestyle=0, thick=1.2, color=COLORS[j];,/ylog
  oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog  
 endelse
endfor
;SYM=-SYM
legend,file_basename(DIR1),psym=SYM, linestyle=LINE,charsize=0.8, color=COLORS
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;***********Plots all AWD
toplot=d.awd
;for j=0,nd-1 do toplot[*,j]/=mean(toplot[d[j].np-3:d[j].np-1,j])/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/awd_vs_time.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

xx = indgen(150)*0.1

 for j=0, nd-1 do begin
  time=str2utc(d[j].date[0:d[j].np-1])
  time=float(time.time)/1000./3600. + (time.mjd-time[0].mjd)*24. ; date to h
  x=time-tini[j]
  y=toplot[0:d[j].np-1,j]
 ; d[j].np=4 ; use only first d[j].np time points 
 if j eq 0 then begin
  ; fit1s a power law (linear in log space)
   fit_awdt=mpfitfun('paper_plots_func2fit',x,y,5.+(y-y),func2fitINI)
   ;chi=sqrt(chi/n_elements(x))  
   ;yrange=[min(toplot[0:d[j].np-1,*],/nan),max(toplot[0:d[j].np-1,*],/nan)]
   ;xrange=[min(d[*].hgt[*],/nan),max(d[*].hgt[*],/nan)]
   plot, x, y, psym=-SYM[j], linestyle=LINE[j],title= 'Axial angular width', $
      xtitle='Time from start [h]', ytitle='AWD', yrange=AWRANGE,$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=HRANGE;,/ylog
 endif else begin
  fit_awdt=[fit_awdt,mpfitfun('paper_plots_func2fit',x,y,5.+(y-y),func2fitINI)]
  ;chi=[chi,sqrt(chi/n_elements(x))]
  oplot, x, y, psym=-SYM[j], linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
 endelse
 oplot, xx, paper_plots_func2fit(xx,[fit_awdt[3*j:3*j+2]]),linestyle=LINE[j],color=COLORS[j] ;, psym=SYM[j], linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
endfor
;print, fit_awdh;, chi
legend,file_basename(DIR1),psym=SYM, linestyle=LINE,charsize=0.8, color=COLORS
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice


;***********Plots all AWD
toplot=d.awl
;for j=0,nd-1 do toplot[*,j]/=mean(toplot[d[j].np-3:d[j].np-1,j])/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/awl_vs_time.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

xx = indgen(150)*0.1
 for j=0, nd-1 do begin
  time=str2utc(d[j].date[0:d[j].np-1])
  time=float(time.time)/1000./3600. + (time.mjd-time[0].mjd)*24. ; date to h
  x=time-tini[j]
  y=toplot[0:d[j].np-1,j]
 ; d[j].np=4 ; use only first d[j].np time points 
 if j eq 0 then begin
  ; fit1s a power law (linear in log space)
   fit_awlt=mpfitfun('paper_plots_func2fit',x,y,5.+(y-y),func2fitINI)
   ;chi=sqrt(chi/n_elements(x))  
   ;yrange=[min(toplot[0:d[j].np-1,*],/nan),max(toplot[0:d[j].np-1,*],/nan)]
   ;xrange=[min(d[*].hgt[*],/nan),max(d[*].hgt[*],/nan)]
   plot, x, y, psym=-SYM[j], linestyle=LINE[j],title= 'Lateral angular width', $
      xtitle='Time from start [h]', ytitle='AWL', yrange=AWRANGE,$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=HRANGE;,/ylog
 endif else begin
  fit_awlt=[fit_awlt,mpfitfun('paper_plots_func2fit',x,y,5.+(y-y),func2fitINI)]
  ;chi=[chi,sqrt(chi/n_elements(x))]
  oplot, x, y, psym=-SYM[j], linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
 endelse
 oplot, xx, paper_plots_func2fit(xx,[fit_awlt[3*j:3*j+2]]),linestyle=LINE[j],color=COLORS[j] ;, psym=SYM[j], linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
endfor
;print, fit_awdh;, chi
legend,file_basename(DIR1),psym=SYM, linestyle=LINE,charsize=0.8, color=COLORS
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;***********Plots
;for j=0,nd-1 do toploty[*,j]/=mean(toploty[d[j].np-3:d[j].np-1,j])/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/dawd_vs_height.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

xx=indgen(200)*0.1

for j=0, nd-1 do begin
time=str2utc(d[j].date[0:d[j].np-1])
time=float(time.time)/1000./3600. + (time.mjd-time[0].mjd)*24. ; date to h
;time-=tini[j]
toploty=paper_plots_func2fit_deriv(time-tini[j],[fit_awdt[3*j:3*j+2]])
toplotx=paper_plots_func2fit4(time,[fit_ht[3*j:3*j+2]])

toplotxx=paper_plots_func2fit4(xx+tini[j],[fit_ht[3*j:3*j+2]])
toplotyy=paper_plots_func2fit_deriv(xx,[fit_awdt[3*j:3*j+2]])

 if j eq 0 then begin
   plot, toplotx, toploty, psym=-SYM[j], linestyle=LINE[j],title= 'Axial angular width change rate for all events', $
     xtitle='Height [Rs]', ytitle='Change rate [deg/h]', yrange=[0.12,120], /ylog,$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=HRANGE ;,/ylog
   oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
 endif else begin
  oplot, toplotx, toploty, psym=-SYM[j], linestyle=0, thick=1.2, color=COLORS[j];,/ylog
  oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog  
 endelse
endfor
;SYM=-SYM
legend,file_basename(DIR1),psym=SYM, linestyle=LINE,charsize=0.8, color=COLORS, position=[6,400]
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;***********Plots
;for j=0,nd-1 do toploty[*,j]/=mean(toploty[d[j].np-3:d[j].np-1,j])/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/dawl_vs_height.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

xx=indgen(200)*0.1

for j=0, nd-1 do begin
time=str2utc(d[j].date[0:d[j].np-1])
time=float(time.time)/1000./3600. + (time.mjd-time[0].mjd)*24. ; date to h
;time-=tini[j]
toploty=paper_plots_func2fit_deriv(time-tini[j],[fit_awlt[3*j:3*j+2]])
toplotx=paper_plots_func2fit4(time,[fit_ht[3*j:3*j+2]])

toplotxx=paper_plots_func2fit4(xx+tini[j],[fit_ht[3*j:3*j+2]])
toplotyy=paper_plots_func2fit_deriv(xx,[fit_awlt[3*j:3*j+2]])

 if j eq 0 then begin
   plot, toplotx, toploty, psym=-SYM[j], linestyle=LINE[j],title= 'Lateral angular width change rate for all events', $
     xtitle='Height [Rs]', ytitle='Change rate [deg/h]', yrange=[0.12,120], /ylog,$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=HRANGE ;,/ylog
   oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
 endif else begin
  oplot, toplotx, toploty, psym=-SYM[j], linestyle=0, thick=1.2, color=COLORS[j];,/ylog
  oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog  
 endelse
endfor
;SYM=-SYM
legend,file_basename(DIR1),psym=SYM, linestyle=LINE,charsize=0.8, color=COLORS, position=[6,400]
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice


;************plots final AWs
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/finalAW_vs_radAcc.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

x=reform(fit_ht,[3,nd])
x=x[0,*]*2.*RS2KM/(3600.^2.)*1000. ; radial acceleration [m/s^2]
y=reform(fit_awdh,[3,nd])
y=y[0,*] ; final AWD
y1=reform(fit_awlh,[3,nd])
y1=y1[0,*]; final AWL

plot, x, y, psym=5, linestyle=0,title= 'Final angular widths vs radial acceleration for all events', $
     xtitle='Radial acc. [km/s^2]', ytitle='Final angular width [deg]', yrange=[min([y,y1],/nan),1.2*max([y,y1],/nan)],$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[0], xrange=[min(x,/nan),max(x,/nan)];,/ylog
oplot, x, y1, psym=2, linestyle=0, thick=1.2, color=COLORS[3];,/ylog
;SYM=-SYM
legend,['Final AWD','Final AWL'],psym=[5,2], linestyle=[0,0],charsize=0.8, color=[COLORS[0],COLORS[3]]
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;************plots Settling height for AWs
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/setHeight_vs_radAcc.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

x=reform(fit_ht,[3,nd])
x=x[0,*]*2.*RS2KM/(3600.^2.)*1000. ; radial acceleration [m/s^2]
y=paper_plots_setHeight(fit_awdh[0:2])
y1=paper_plots_setHeight(fit_awlh[0:2])
for j=1, nd-1 do begin ; Settling heights
  y =[y ,paper_plots_setHeight(fit_awdh[3*j:3*j+2])]
  y1=[y1,paper_plots_setHeight(fit_awlh[3*j:3*j+2])]
endfor
;y[where(y gt 50)] = 30
;y1[where(y1 gt 50)] = 30
plot, x, y, psym=5, linestyle=0,title= 'Settling heights vs radial acceleration for all events', $
     xtitle='Radial acc. [km/s^2]', ytitle='Settling height [Rs]', yrange=[min([y,y1],/nan),1.2*max([y,y1],/nan)],$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[0], xrange=[min(x,/nan),max(x,/nan)];,/ylog
oplot, x, y1, psym=2, linestyle=0, thick=1.2, color=COLORS[3];,/ylog
;SYM=-SYM
legend,['AWD','AWL'],psym=[5,2], linestyle=[0,0],charsize=0.8, color=[COLORS[0],COLORS[3]]
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;************plots Settling height for AWs
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/setHeight_vs_finalAW.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

x=reform(fit_awdh,[3,nd])
x=x[0,*] ; final AWD
x1=reform(fit_awlh,[3,nd])
x1=x1[0,*]; final AWL

y=paper_plots_setHeight(fit_awdh[0:2])
y1=paper_plots_setHeight(fit_awlh[0:2])
for j=1, nd-1 do begin ; Settling heights
  y =[y ,paper_plots_setHeight(fit_awdh[3*j:3*j+2])]
  y1=[y1,paper_plots_setHeight(fit_awlh[3*j:3*j+2])]
endfor
;y[where(y gt 50)] = 30
;y1[where(y1 gt 50)] = 30
plot, x, y, psym=5, linestyle=0,title= 'Settling heights vs final angular widths for all events', $
     xtitle='Final angular width [deg]', ytitle='Settling height [Rs]', yrange=[min([y,y1],/nan),1.2*max([y,y1],/nan)],$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[0], $
        xrange=[min([x,x1],/nan),max([x,x1],/nan)];,/ylog,/xlog
oplot, x1, y1, psym=2, linestyle=0, thick=1.2, color=COLORS[3]
;SYM=-SYM
legend,['AWD','AWL'],psym=[5,2], linestyle=[0,0],charsize=0.8, color=[COLORS[0],COLORS[3]]
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;************plots Settling height for AWs
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/awl_vs_awd_setHeight.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

; x=reform(fit_awdh,[3,nd])
; x=x[0,*] ; final AWD
; x1=reform(fit_awlh,[3,nd])
; x1=x1[0,*]; final AWL

y=paper_plots_setHeight(fit_awdh[0:2])
x=paper_plots_setHeight(fit_awlh[0:2])
for j=1, nd-1 do begin ; Settling heights
  y =[y ,paper_plots_setHeight(fit_awdh[3*j:3*j+2])]
  x=[x,paper_plots_setHeight(fit_awlh[3*j:3*j+2])]
endfor
;y[where(y gt 50)] = 30
;y1[where(y1 gt 50)] = 30
plot, [0,0], [0,0], psym=-1, linestyle=0,title= 'Settling height for all events', $
     xtitle='AWD height [Rs]', ytitle='AWL height [Rs]', yrange=[2,6],$
      font=16, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[0], $
        xrange=[2,6];,/ylog,/xlog
for j=0, nd-1 do oplot, [0,y[j]], [0,x[j]], psym=-SYM[j],color=COLORS[j], thick=5,symsize=1.4
oplot,[2,6],[2,6],linestyle=0
legend,file_basename(DIR1),psym=-SYM,charsize=1.2,color=COLORS, thick=5, position=[5,4.3]
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;************plots Settling height for AWs
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/awl_vs_awd_finalAW.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

x=reform(fit_awdh,[3,nd])
x=x[0,*] ; final AWD
x1=reform(fit_awlh,[3,nd])
x1=x1[0,*]; final AWL

;y[where(y gt 50)] = 30
;y1[where(y1 gt 50)] = 30
plot, [0,0], [0,0], psym=-1, linestyle=0,title= 'Final AW for all events', $
     xtitle='Final AWD [deg]', ytitle='Final AWL [deg]', yrange=[30,160],$
      font=16, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[0], $
        xrange=[30,60],xstyle=1,ystyle=1;,/ylog,/xlog
for j=0, nd-1 do oplot, [0,x[j]], [0,x1[j]], psym=-SYM[j],color=COLORS[j], thick=5,symsize=1.4
oplot,[0,160],[0,160],linestyle=0
oplot,[0,80],[0,160],linestyle=0
legend,file_basename(DIR1),psym=-SYM,charsize=1.2,color=COLORS, thick=5, position=[31,140]
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;************plots the fits parameters
; mydevice = !D.NAME
; SET_PLOT, 'PS'
; DEVICE, FILENAME=OPATH + '/pram_hist.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
; !P.Multi = 0
; FONT=FONTSZ
; MARGIN=[5,5]
; YLOG=0
; !P.Color = '000000'xL
; !P.Background = 'FFFFFF'xL
;loadct, COLORT

; fit=reform(fit,[3,nd])
; fit1=reform(fit1,[3,nd])
; ;h0=histogram(fit[0,*], locations=loc)
; plot, fit[0,*],fit1[0,*], xrange=[0,190],yrange=[0,190],psym=-SYM[0], xtitle='P0(AWL)', ytitle='P0(AWD)', title='Final angular widths'
; DEVICE, /CLOSE
; SET_PLOT, mydevice

; ;************************* individual AWL and AWD plots 
; toplot=d.awl
; toplot1=d.awd
; xx = indgen(150)*0.1
; ;for j=0,nd-1 do toplot[*,j]/=mean(toplot[d[j].np-3:d[j].np-1,j])/100.
; for j=0, nd-1 do begin
;  SET_PLOT, 'PS'
;  DEVICE, FILENAME=OPATH + '/ind_AW_'+file_basename(DIR[j])+'_noNorm.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ  
;  !P.Multi = 0
;  FONT=FONTSZ
;  MARGIN=[5,5]
;  YLOG=0
;  !P.Color = '000000'xL
;  !P.Background = 'FFFFFF'xL
;loadct, COLORT  
;  x=d[j].hgt[0:d[j].np-1]
;  y=toplot[0:d[j].np-1,j]
;  y1=toplot1[0:d[j].np-1,j]
;  fit=mpfitfun('paper_plots_func2fit',x,y,5.+(y-y),func2fitINI)
;  ;chi=sqrt(chi/n_elements(x))  
;  plot, x, y, psym=-SYM[0], linestyle=LINE[0],title= 'Angular widths evolution', $
;     xtitle='Height [Rs]', ytitle='Angular width [deg]', yrange=[min([y,y1],/nan),max([y,y1],/nan)],$
;     font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=[min(d[*].hgt[*],/nan),max(d[*].hgt[*],/nan)];,/ylog
;  oplot, xx, paper_plots_func2fit(xx,fit),linestyle=LINE[0] ;, psym=SYM[j], linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
;  xt=median(d[*].hgt[*])
;  yt=median([y,y1])
;  bla=strtrim(fit[0],2)+' ; '+strtrim(fit[1],2)+' ; '+strtrim(fit[2],2)

;  fit=mpfitfun('paper_plots_func2fit',x,y1,5.+(y-y),func2fitINI)
;  ;chi=sqrt(chi/n_elements(x))  
;  oplot, x, y1, psym=-SYM[1], linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
;  oplot, xx, paper_plots_func2fit(xx,fit),linestyle=LINE[0] ;, psym=SYM[j], linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
;  bla1=strtrim(fit[0],2)+' ; '+strtrim(fit[1],2)+' ; '+strtrim(fit[2],2)
;  legend,['AWL:'+bla, 'AWD:'+bla1], psym=SYM[0:1], linestyle=LINE[0:1],charsize=0.8
;  DEVICE, /CLOSE
;  SET_PLOT, mydevice
; endfor


;***********Plots all wl vs time
toplot=d.wl
;for j=0,nd-1 do toplot[*,j]/=mean(toplot[d[j].np-3:d[j].np-1,j])/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/wl_vs_time.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT
pi=replicate({limited:[0,0], limits:[0.D,0.D], tied:'', fixed:0},3)    ;creates structure for PARINFO
pi[0].limited[0]=1 ; enables lower limit
pi[0].limits[0]=0; lower limit
pi[1].limited[0]=1 
pi[1].limits[0]=0
pi[2].fixed=1

for j=0, nd-1 do begin
 ;d[j].np=4 ; use only first d[j].np time points 
 time=str2utc(d[j].date[0:d[j].np-1])
 time=float(time.time)/1000./3600. + (time.mjd-time[0].mjd)*24. ; date to h
 time-=tini[j]
 if j eq 0 then begin
   fit_wlt=mpfitfun('paper_plots_func2fit5',time,toplot[0:d[j].np-1,j],0.5,func2fitINI5,PARINFO=pi)
   plot, time, toplot[0:d[j].np-1,j], psym=-SYM[j], linestyle=LINE[j],title= 'Lateral width vs time for all events', $
      xtitle='time from start [h]', ytitle='WL [Rs]', yrange=[min(d[*].hgt[*],/nan),max(d[*].hgt[*],/nan)],$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=[0,7];,/ylog
 endif else begin
  fit_wlt=[fit_wlt,mpfitfun('paper_plots_func2fit5',time,toplot[0:d[j].np-1,j],0.5,func2fitINI5,PARINFO=pi)]  ;
  oplot, time, toplot[0:d[j].np-1,j], psym=-SYM[j], linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
 endelse
  oplot, xx, paper_plots_func2fit5(xx,[fit_wlt[3*j:3*j+2]]),linestyle=LINE[j],color=COLORS[j] 
endfor
legend,file_basename(DIR1),psym=SYM, linestyle=LINE,charsize=0.8,color=COLORS;,position=[5,3.5]
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;***********Plots all wd vs time
toplot=d.wd
;for j=0,nd-1 do toplot[*,j]/=mean(toplot[d[j].np-3:d[j].np-1,j])/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/wd_vs_time.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT
pi=replicate({limited:[0,0], limits:[0.D,0.D], tied:'', fixed:0},3)    ;creates structure for PARINFO
pi[0].limited[0]=1 ; enables lower limit
pi[0].limits[0]=0; lower limit
pi[1].limited[0]=1 
pi[1].limits[0]=0
pi[2].fixed=1

for j=0, nd-1 do begin
 ;d[j].np=4 ; use only first d[j].np time points 
 time=str2utc(d[j].date[0:d[j].np-1])
 time=float(time.time)/1000./3600. + (time.mjd-time[0].mjd)*24. ; date to h
 time-=tini[j]
 if j eq 0 then begin
   fit_wdt=mpfitfun('paper_plots_func2fit5',time,toplot[0:d[j].np-1,j],0.5,func2fitINI5,PARINFO=pi)
   plot, time, toplot[0:d[j].np-1,j], psym=-SYM[j], linestyle=LINE[j],title= 'Axial width vs time for all events', $
      xtitle='time from start [h]', ytitle='WD [Rs]', yrange=[min(d[*].hgt[*],/nan),max(d[*].hgt[*],/nan)],$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=[0,7];,/ylog
 endif else begin
  fit_wdt=[fit_wdt,mpfitfun('paper_plots_func2fit5',time,toplot[0:d[j].np-1,j],0.5,func2fitINI5,PARINFO=pi)]  
  oplot, time, toplot[0:d[j].np-1,j], psym=-SYM[j], linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
 endelse
  oplot, xx, paper_plots_func2fit5(xx,[fit_wdt[3*j:3*j+2]]),linestyle=LINE[j],color=COLORS[j] 
endfor
legend,file_basename(DIR1),psym=SYM, linestyle=LINE,charsize=0.8,color=COLORS;,position=[5,3.5]
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
DEVICE, FILENAME=OPATH + '/dawd_vs_dawl.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

for j=0, nd-1 do begin
 ;d[j].np=4 ; use only first d[j].np time points 
 if j eq 0 then begin
   plot, toplotx[0:d[j].np-1,j], toploty[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j],title= 'Angular speeds for all events', $
      xtitle='deriv. AWD [deg/s]', ytitle='deriv. AWL [deg/s]', yrange=[min(toploty[0:d[j].np-1,*],/nan),1.2*max(toploty[0:d[j].np-1,*],/nan)],$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=[min(toplotx,/nan),max(toplotx,/nan)];,/ylog
 endif else begin
  oplot, toplotx[0:d[j].np-1,j], toploty[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
 endelse
endfor
xx=indgen(50)*(max(toplotx,/nan)-min(toplotx,/nan))/50+min(toplotx,/nan)
oplot, xx,fit[0]+fit[1]*xx, thick=7;,/ylog
;SYM=-SYM
legend,[file_basename(DIR1),'FIT m = '+strtrim(fit[1],2)+'; rms_err='+strtrim(chisqr,2)],psym=[SYM,0], $
  linestyle=[LINE,0], /left,charsize=0.8,color=[COLORS,COLORS[1]]
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;***********Plots
;for j=0,nd-1 do toploty[*,j]/=mean(toploty[d[j].np-3:d[j].np-1,j])/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/awl_over_awd_vs_height.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

mxx=indgen(100)*0.1

for j=0, nd-1 do begin
 time=str2utc(d[j].date[0:d[j].np-1])
 time=float(time.time)/1000./3600. + (time.mjd-time[0].mjd)*24. ; date to h
 toplotx=paper_plots_func2fit4(time,[fit_ht[3*j:3*j+2]])  
 toploty=paper_plots_func2fit(toplotx,[fit_awlh[3*j:3*j+2]])/paper_plots_func2fit(toplotx,[fit_awdh[3*j:3*j+2]])

 toplotxx= mxx+ toplotx[0]
 toplotyy=paper_plots_func2fit(toplotxx,[fit_awlh[3*j:3*j+2]])/paper_plots_func2fit(toplotxx,[fit_awdh[3*j:3*j+2]])
 
 if j eq 0 then begin
  plot, toplotx, toploty, psym=-SYM[j], linestyle=0,title= 'Angular width for all events', $
     xtitle='Height [Rs]', ytitle='AWL/AWD', yrange=[1.1, 4.0],$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=HRANGE
  oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
  ; saves for the mean curve
  myy=paper_plots_func2fit(mxx,[fit_awlh[3*j:3*j+2]])/paper_plots_func2fit(mxx,[fit_awdh[3*j:3*j+2]])
  mini=toplotxx[0]
 endif else begin
  oplot, toplotx, toploty, psym=-SYM[j], linestyle=0, thick=1.2, color=COLORS[j];,/ylog
  oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog 
  myy=[[myy],[paper_plots_func2fit(mxx,[fit_awlh[3*j:3*j+2]])/paper_plots_func2fit(mxx,[fit_awdh[3*j:3*j+2]])]] 
  mini=[mini,toplotxx[0]]
 endelse
endfor
allm=mxx
allsd=mxx
for r=0, n_elements(allm)-1 do begin
  ind=where(mini lt mxx[r], count)
  if count eq 12 then begin
    allm[r]=median(myy[r,ind])
    allsd[r]=stddev(myy[r,ind])
  endif else begin
    allm[r]=!values.f_nan
    allsd[r]=!values.f_nan
  endelse
endfor
;oploterror, mxx, allm, allsd,ERRTHICK=0.5,Nskip=8, psym=0, linestyle=0, thick=7, color=0;,/ylog 
;oplot, mxx, allm+allsd, psym=0, linestyle=1, thick=9, color=0;,/ylog 
ind=where(~finite(allm))
;allm=smooth(allm,10,/nan)
;allm[ind]=!values.f_nan
oplot, mxx, allm, psym=0, linestyle=0, thick=9, color=0;,/ylog 
;oplot, mxx, allm-allsd, psym=0, linestyle=1, thick=9, color=0;,/ylog 
print, strtrim(mean(allm,/nan)) + '+-' + strtrim(mean(allsd,/nan))
;SYM=-SYM
legend,file_basename(DIR1),psym=SYM,/fill, linestyle=LINE,charsize=0.8, color=COLORS,$
  position=[5,3.95], spacing=1
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;***********Plots
toploty=d.wl
toplotx=d.hgt
;for j=0,nd-1 do toploty[*,j]/=mean(toploty[d[j].np-3:d[j].np-1,j])/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/wl_vs_height.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

for j=0, nd-1 do begin
 ;d[j].np=4 ; use only first d[j].np time points 
 if j eq 0 then begin
   plot, toplotx[0:d[j].np-1,j], toploty[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j],title= 'Lateral width vs Height for all events', $
     xtitle='Height [Rs]', ytitle='WL [Rs]', yrange=[0,15],$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=[0,15];,/ylog
 endif else begin
  oplot, toplotx[0:d[j].np-1,j], toploty[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
 endelse
endfor
;SYM=-SYM
legend,file_basename(DIR1),psym=SYM, linestyle=LINE,charsize=0.8, color=COLORS
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice


;***********Plots
toploty=d.wd
toplotx=d.hgt
;for j=0,nd-1 do toploty[*,j]/=mean(toploty[d[j].np-3:d[j].np-1,j])/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/wd_vs_height.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

for j=0, nd-1 do begin
 ;d[j].np=4 ; use only first d[j].np time points 
 if j eq 0 then begin
   plot, toplotx[0:d[j].np-1,j], toploty[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j],title= 'Axial width vs Height for all events', $
     xtitle='Height [Rs]', ytitle='WD [Rs]', yrange=[0,15],$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=[0,15];,/ylog
 endif else begin
  oplot, toplotx[0:d[j].np-1,j], toploty[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
 endelse
endfor
;SYM=-SYM
legend,file_basename(DIR1),psym=SYM, linestyle=LINE,charsize=0.8, color=COLORS
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice


;***********Plots
toploty=d.wl
toplotx=d.wd
;for j=0,nd-1 do toploty[*,j]/=mean(toploty[d[j].np-3:d[j].np-1,j])/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/wl_vs_wd.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL   
!P.Background = 'FFFFFF'xL
loadct, COLORT

for j=0, nd-1 do begin
 ;d[j].np=4 ; use only first d[j].np time points 
 if j eq 0 then begin
   plot, toplotx[0:d[j].np-1,j], toploty[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j],title= 'Axial vs Lateral widths for all events', $
     xtitle='WD [Rs]', ytitle='WL [Rs]', yrange=[0,15],$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=[0,15];,/ylog
 endif else begin
  oplot, toplotx[0:d[j].np-1,j], toploty[0:d[j].np-1,j], psym=SYM[j], linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
 endelse
endfor
;SYM=-SYM
legend,file_basename(DIR1),psym=SYM, linestyle=LINE,charsize=0.8, color=COLORS
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;***********Plots
xx = indgen(10)*4.
;for j=0,nd-1 do toploty[*,j]/=mean(toploty[d[j].np-3:d[j].np-1,j])/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/dwl_vs_dheight.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

for j=0, nd-1 do begin
time=str2utc(d[j].date[0:d[j].np-1])
time=float(time.time)/1000./3600. + (time.mjd-time[0].mjd)*24. ; date to h

toploty=paper_plots_func2fit2_deriv(time-tini[j],[fit_wlt[3*j:3*j+2]])*RS2KM/3600.
toplotx=paper_plots_func2fit2_deriv(time,[fit_ht[3*j:3*j+2]])*RS2KM/3600.  

xx=indgen(200)*0.1
toplotxx=paper_plots_func2fit2_deriv(xx,[fit_ht[3*j:3*j+2]])*RS2KM/3600.  
toplotyy=paper_plots_func2fit2_deriv(xx-tini[j],[fit_wlt[3*j:3*j+2]])*RS2KM/3600.

 if j eq 0 then begin
   plot, toplotx, toploty, psym=-SYM[j], linestyle=LINE[j],title= 'Lateral vs Radial velocity for all events', $
     xtitle='Radial vel. [km/s]', ytitle='Lateral expansion vel. [km/s]', yrange=VELRANGE,$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=VELRANGE;,/ylog
   oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
 endif else begin
  oplot, toplotx, toploty, psym=-SYM[j], linestyle=0, thick=1.2, color=COLORS[j];,/ylog
  oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog  
 endelse
 ;if j eq 7 then stop
endfor

;SYM=-SYM
legend,file_basename(DIR1),psym=SYM, linestyle=LINE,charsize=0.8, color=COLORS
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;***********Plots
xx = indgen(10)*4.
;for j=0,nd-1 do toploty[*,j]/=mean(toploty[d[j].np-3:d[j].np-1,j])/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/dwl_over_dheight_vs_height.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

for j=0, nd-1 do begin
time=str2utc(d[j].date[0:d[j].np-1])
time=float(time.time)/1000./3600. + (time.mjd-time[0].mjd)*24. ; date to h

toploty=paper_plots_func2fit2_deriv(time-tini[j],[fit_wlt[3*j:3*j+2]])/$
        paper_plots_func2fit2_deriv(time,[fit_ht[3*j:3*j+2]])
toplotx=paper_plots_func2fit4(time,[fit_ht[3*j:3*j+2]])  

xx=indgen(200)*0.1
toplotxx=paper_plots_func2fit4(xx+time[0],[fit_ht[3*j:3*j+2]]) 
toplotyy=paper_plots_func2fit2_deriv(xx+time[0]-tini[j],[fit_wlt[3*j:3*j+2]])/$
  paper_plots_func2fit2_deriv(xx+time[0],[fit_ht[3*j:3*j+2]])

 if j eq 0 then begin
   plot, toplotx, toploty, psym=-SYM[j], linestyle=LINE[j], $
     xtitle='Height[Rs]', ytitle='Lateral velocity / Radial velocity',ystyle=1, yrange=[0.3,3.1],$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=HRANGE;,/ylog
   oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
  myy=toplotyy
  mxx=toplotxx
 endif else begin
  oplot, toplotx, toploty, psym=-SYM[j], linestyle=0, thick=1.2, color=COLORS[j];,/ylog
  oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog  
  myy=[[myy],[toplotyy]]
  mxx=[[mxx],[toplotxx]]
 endelse
 ;print,file_basename(DIR1[j])
 ;stop
endfor
allx=indgen(100)*0.1
allm=allx
allsd=allx
for r=0, n_elements(allx)-1 do begin
  ind=where(mxx[0,*] le allx[r], count)
  if count eq 12 then begin
    blam=interpol(myy[*,ind[0]],mxx[*,ind[0]],allx[r])
    for j=1, count-1 do blam=[blam,interpol(myy[*,ind[j]],mxx[*,ind[j]],allx[r])]
    allm[r]=mean(blam);/=count
    allsd[r]=stddev(blam)
  endif else begin
    allm[r]=!values.f_nan
    allsd[r]=!values.f_nan
  endelse
endfor

print,'wl_over_dheight_vs_height:'
print, strtrim(mean(allm,/nan)) + '+-' + strtrim(mean(allsd,/nan))

ind=where(~finite(allm))
;allm=smooth(allm,10,/nan)
;allm[ind]=!values.f_nan
oplot, allx, allm, psym=0, linestyle=0, thick=9, color=0;,/ylog 
;SYM=-SYM
legend,file_basename(DIR1),psym=SYM, linestyle=LINE,charsize=0.8, color=COLORS,position=[4.5,2.7]
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;***********Plots
;for j=0,nd-1 do toploty[*,j]/=mean(toploty[d[j].np-3:d[j].np-1,j])/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/dwd_vs_dheight.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

for j=0, nd-1 do begin
time=str2utc(d[j].date[0:d[j].np-1])
time=float(time.time)/1000./3600. + (time.mjd-time[0].mjd)*24. ; date to h

toploty=paper_plots_func2fit2_deriv(time-tini[j],[fit_wdt[3*j:3*j+2]])*RS2KM/3600.
toplotx=paper_plots_func2fit2_deriv(time,[fit_ht[3*j:3*j+2]])*RS2KM/3600.  

xx=indgen(200)*0.1
toplotxx=paper_plots_func2fit2_deriv(xx,[fit_ht[3*j:3*j+2]])*RS2KM/3600.  
toplotyy=paper_plots_func2fit2_deriv(xx-tini[j],[fit_wdt[3*j:3*j+2]])*RS2KM/3600.

t6=root_paper_plots_func2fit4(6., fit_ht[3*j:3*j+2]) ; time at 6 Rs
print, 'DWD at 6 Rs',paper_plots_func2fit2_deriv(t6-tini[j],[fit_wdt[3*j:3*j+2]])*RS2KM/3600.
 if j eq 0 then begin
   plot, toplotx, toploty, psym=-SYM[j], linestyle=LINE[j],title= 'Axial vs Radial velocity for all events', $
     xtitle='Radial vel. [km/s]', ytitle='Axial expansion vel. [km/s]', yrange=VELRANGE,$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=VELRANGE;,/ylog
   oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
 endif else begin
  oplot, toplotx, toploty, psym=-SYM[j], linestyle=0, thick=1.2, color=COLORS[j];,/ylog
  oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog  
 endelse
endfor
;SYM=-SYM
legend,file_basename(DIR1),psym=SYM, linestyle=LINE,charsize=0.8, color=COLORS
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;***********Plots
xx = indgen(10)*4.
;for j=0,nd-1 do toploty[*,j]/=mean(toploty[d[j].np-3:d[j].np-1,j])/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/dwl_vs_height.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

for j=0, nd-1 do begin
time=str2utc(d[j].date[0:d[j].np-1])
time=float(time.time)/1000./3600. + (time.mjd-time[0].mjd)*24. ; date to h

toploty=paper_plots_func2fit2_deriv(time-tini[j],[fit_wlt[3*j:3*j+2]])*RS2KM/3600.
toplotx=paper_plots_func2fit4(time,[fit_ht[3*j:3*j+2]])  

xx=indgen(200)*0.1+tini[j]
toplotxx=paper_plots_func2fit4(xx,[fit_ht[3*j:3*j+2]]) 
toplotyy=paper_plots_func2fit2_deriv(xx-tini[j],[fit_wlt[3*j:3*j+2]])*RS2KM/3600.

t6=root_paper_plots_func2fit4(6., fit_ht[3*j:3*j+2]) ; time at 6 Rs
print, 'DWL at 6 Rs',paper_plots_func2fit2_deriv(t6-tini[j],[fit_wlt[3*j:3*j+2]])*RS2KM/3600. 

 if j eq 0 then begin
   plot, toplotx, toploty, psym=-SYM[j], linestyle=LINE[j],title= 'Lateral vs Radial velocity for all events', $
     xtitle='Height[Sr]', ytitle='Lateral expansion vel. [km/s]', yrange=VELRANGE,$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=HRANGE;,/ylog
   oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
 endif else begin
  oplot, toplotx, toploty, psym=-SYM[j], linestyle=0, thick=1.2, color=COLORS[j];,/ylog
  oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog  
 endelse
 ;if j eq 7 then stop
endfor

;SYM=-SYM
legend,file_basename(DIR1),psym=SYM, linestyle=LINE,charsize=0.8, color=COLORS
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice


;***********Plots
;for j=0,nd-1 do toploty[*,j]/=mean(toploty[d[j].np-3:d[j].np-1,j])/100.
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/dwd_vs_height.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

for j=0, nd-1 do begin
time=str2utc(d[j].date[0:d[j].np-1])
time=float(time.time)/1000./3600. + (time.mjd-time[0].mjd)*24. ; date to h

toploty=paper_plots_func2fit2_deriv(time-tini[j],[fit_wdt[3*j:3*j+2]])*RS2KM/3600.
toplotx=paper_plots_func2fit4(time,[fit_ht[3*j:3*j+2]])  

xx=indgen(200)*0.1+tini[j]
toplotxx=paper_plots_func2fit4(xx,[fit_ht[3*j:3*j+2]]) 
toplotyy=paper_plots_func2fit2_deriv(xx-tini[j],[fit_wdt[3*j:3*j+2]])*RS2KM/3600.

 if j eq 0 then begin
   plot, toplotx, toploty, psym=-SYM[j], linestyle=LINE[j],title= 'Axial vs Radial velocity for all events', $
     xtitle='Height [Rs]', ytitle='Axial expansion vel. [km/s]', yrange=VELRANGE,$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=HRANGE ;,/ylog
   oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
 endif else begin
  oplot, toplotx, toploty, psym=-SYM[j], linestyle=0, thick=1.2, color=COLORS[j];,/ylog
  oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog  
 endelse
endfor
;SYM=-SYM
legend,file_basename(DIR1),psym=SYM, linestyle=LINE,charsize=0.8, color=COLORS
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;***********Plots
; mean values of the fitted paramerters P[0]
mPod=reform(fit_wdt,3,nd)
mPod=mean(mPod[0,*])
mPol=reform(fit_wlt,3,nd)
mPol=mean(mPol[0,*])

mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/dwl_over_dwd_vs_height.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

for j=0, nd-1 do begin
time=str2utc(d[j].date[0:d[j].np-1])
time=float(time.time)/1000./3600. + (time.mjd-time[0].mjd)*24. ; date to h

toplotx=paper_plots_func2fit4(time,[fit_ht[3*j:3*j+2]])  
toploty=(paper_plots_func2fit2_deriv(time-tini[j],[fit_wlt[3*j:3*j+2]])/$
  paper_plots_func2fit2_deriv(time-tini[j],[fit_wdt[3*j:3*j+2]]))

xx=indgen(200)*0.1+time[0]
toplotxx=paper_plots_func2fit4(xx,[fit_ht[3*j:3*j+2]]) 
toplotyy=(paper_plots_func2fit2_deriv(xx-tini[j],[fit_wlt[3*j:3*j+2]])/$
  paper_plots_func2fit2_deriv(xx-tini[j],[fit_wdt[3*j:3*j+2]]))

 if j eq 0 then begin
   plot, toplotx, toploty, psym=-SYM[j], linestyle=2, $
     xtitle='Height [Rs]', ytitle='Lateral velocity / Axial velocity',ystyle=1,yrange=[1,2.0],$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j], xrange=HRANGE
   oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
  myy=toplotyy
  mxx=toplotxx
 endif else begin
  oplot, toplotx, toploty, psym=-SYM[j], linestyle=0, thick=1.2, color=COLORS[j];,/ylog
  oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog  
  myy=[[myy],[toplotyy]]
  mxx=[[mxx],[toplotxx]]
 endelse
endfor
allx=indgen(100)*0.1
allm=allx
allsd=allx
for r=0, n_elements(allx)-1 do begin
  ind=where(mxx[0,*] le allx[r], count)
  if count eq 12 then begin
    blam=interpol(myy[*,ind[0]],mxx[*,ind[0]],allx[r])
    for j=1, count-1 do blam=[blam,interpol(myy[*,ind[j]],mxx[*,ind[j]],allx[r])]
    allm[r]=mean(blam);/=count
    allsd[r]=stddev(blam)
  endif else begin
    allm[r]=!values.f_nan
    allsd[r]=!values.f_nan
  endelse
endfor

print,'dwl_over_dwd_vs_height:'
print, strtrim(mean(allm,/nan)) + '+-' + strtrim(mean(allsd,/nan))

ind=where(~finite(allm))
;allm=smooth(allm,10,/nan)
;allm[ind]=!values.f_nan
oplot, allx, allm, psym=0, linestyle=0, thick=9, color=0;,/ylog 
;legend,file_basename(DIR1),psym=[SYM], $
;  linestyle=[LINE], /left,charsize=0.8,color=[COLORS];,position=[2,2.45]
;legend,file_basename(DIR1),psym=SYM, linestyle=LINE,charsize=0.8, color=COLORS
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice

;***********Plots
; mean values of the fitted paramerters P[0]
mPod=reform(fit_wdt,3,nd)
mPod=mean(mPod[0,*])
mPol=reform(fit_wlt,3,nd)
mPol=mean(mPol[0,*])

mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/dwl_vs_dwd.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, FONT_SIZE=FONTSZ
!P.Multi = 0
FONT=FONTSZ
MARGIN=[5,5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT

for j=0, nd-1 do begin
time=str2utc(d[j].date[0:d[j].np-1])
time=float(time.time)/1000./3600. + (time.mjd-time[0].mjd)*24. ; date to h
time-=tini[j]
toploty=paper_plots_func2fit2_deriv(time,[fit_wlt[3*j:3*j+2]])*RS2KM/3600.
toplotx=paper_plots_func2fit2_deriv(time,[fit_wdt[3*j:3*j+2]])*RS2KM/3600.  


xx=indgen(200)*0.1
toplotxx=paper_plots_func2fit2_deriv(xx,[fit_wdt[3*j:3*j+2]])*RS2KM/3600.  
toplotyy=paper_plots_func2fit2_deriv(xx,[fit_wlt[3*j:3*j+2]])*RS2KM/3600.

 if j eq 0 then begin
   plot, toplotx, toploty, psym=SYM[j], linestyle=2,title= 'Lateral vs axial velocity for all events', $
     xtitle='Axial expansion vel. [km/s]', ytitle='Lateral expansion vel. [km/s]',ystyle=1, yrange=[20,2000],$
      font=FONT, YMARGIN=MARGIN, XMARGIN=XXMARGIN, charthick=1.6, thick=1.2, color=COLORS[j],xstyle=1, xrange=[20,1000],/ylog,/xlog
   ;oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog
 endif else begin
  oplot, toplotx, toploty, psym=SYM[j], linestyle=2, thick=1.2, color=COLORS[j]
  ;oplot, toplotxx, toplotyy, psym=0, linestyle=LINE[j], thick=1.2, color=COLORS[j];,/ylog  
 endelse
endfor
xx=indgen(100)*15
oplot, xx,mPol / mPod*xx, thick=7;,/ylog
legend,[file_basename(DIR1),'Slope = '+strtrim(mPol / mPod,2)],psym=[SYM,0], $
  linestyle=[LINE,0], /left,charsize=0.8,color=[COLORS,COLORS[1]]
;legend,file_basename(DIR1),psym=SYM, linestyle=LINE,charsize=0.8, color=COLORS
; Close the PostScript file:
DEVICE, /CLOSE
SET_PLOT, mydevice


; *****************histogram of dwl/dwd
mydevice = !D.NAME
SET_PLOT, 'PS'
DEVICE, FILENAME=OPATH + '/dwl_over_dwd_hist.eps', /LANDSCAPE, DECOMPOSED=0,/ENCAPSULATED, COLOR=1, $
XSIZE=11, YSIZE=11, FONT_SIZE=10
!P.Multi = 0
FONT=1.5
MARGIN=[3.5,2.5]
YLOG=0
!P.Color = '000000'xL
!P.Background = 'FFFFFF'xL
loadct, COLORT
j=0

for j=0, nd-1 do begin
time=str2utc(d[j].date[0:d[j].np-1])
time=float(time.time)/1000./3600. + (time.mjd-time[0].mjd)*24. ; date to h
time-=tini[j]
y=paper_plots_func2fit2_deriv(time,[fit_wlt[3*j:3*j+2]])*RS2KM/3600.
x=paper_plots_func2fit2_deriv(time,[fit_wdt[3*j:3*j+2]])*RS2KM/3600.  

fit=linfit(x,y)

if j eq 8 then fit[1]=mean(y)/mean(x)

 if j eq 0 then begin
   dhist=fit[1]
 endif else begin
   dhist=[dhist,fit[1]]
 endelse
endfor

h0=histogram(dhist, locations=loc, nbins=4)
plot,loc, h0, psym=10,linestyle=LINE[0], $
   xtitle='Velocities ratio',ytitle='Number of CMEs',$
   charsize=FONT, YMARGIN=MARGIN,XMARGIN=[6,2.5], xstyle=1,charthick=1.6,symsize=0.5, thick=1.6,$
   yrange=[min(h0),max(h0)+1], color=COLORS[0];,/ylog
;print, dhist, mean(dhist), stddev(dhist), stddev(dhist)/ mean(dhist)
DEVICE, /CLOSE
SET_PLOT, mydevice


end


