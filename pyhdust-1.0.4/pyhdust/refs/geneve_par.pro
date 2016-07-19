function wrot,par,oblat=oblat

;------------------------------------------------
;
; Converts w=Omega/Omega_c into W=vrot/vorb
;
; CALL:
; w=WROT(par[,/OBLAT])
;
;------------------------------------------------

on_error,1
compile_opt idl2

if n_elements(oblat) ne 0 then begin
  w=(1.5d0^1.5d0)*sqrt(2d0*(par-1d0)/par^3)	; Ekstrom et al. 2008, Eq. (9)
endif else w=par

pi=!dpi
gam=2d0*cos((pi+acos(w))/3d0)
ww=sqrt(gam^3/w)

return,ww
end


;::::::::::::::::::::::::::::::::::::::::::::::::

pro beta_hdust,par,beta,oblat=oblat

;------------------------------------------------
;
; Beta Espinosa-Lara (VLTI School lecture)
;
;------------------------------------------------

on_error,1
compile_opt idl2

; Ekstrom et al. 2008, Eq. 9
if n_elements(oblat) ne 0 then begin
  omega_c=(1.5d0^1.5d0)*sqrt(2d0*(par-1d0)/par^3)
endif else omega_c=par

; Espinosa-Lara lecture, slide 18
delt=1.
omega1=0d0
omega=omega_c
while delt ge 1d-5 do begin
  f=(3d0/(2d0+omega^2))^3*omega^2-omega_c^2
  df=-108d0*omega*(omega^2-1d0)/(omega^2+2d0)^4
  omega1=omega-f/df
  delt=abs(omega1-omega)/omega
  omega=omega1
endwhile

nom=1
omega=[omega]

pi=!dpi
nthe=100
theta=(pi/2d0)*(findgen(nthe)+1d0)/nthe
grav=fltarr(nom,nthe)
teff=fltarr(nom,nthe)
corr=fltarr(nom,nthe)
beta=fltarr(nom)

for iom=0,nom-1 do begin
  for ithe=0,nthe-1 do begin

    delt=1.
    r1=0d0
    r=1d0
    while delt ge 1d-5 do begin
      f=omega[iom]^2*r^3*sin(theta[ithe])^2-(2d0+omega[iom]^2)*r+2d0
      df=3d0*omega[iom]^2*r^2*sin(theta[ithe])^2-(2d0+omega[iom]^2)
      r1=r-f/df
      delt=abs(r1-r)/r
      r=r1
    endwhile

    delt=1.
    n1=0d0
    ftheta=1d0/3d0*omega[iom]^2*r^3*cos(theta[ithe])^3$
           +cos(theta[ithe])+alog(tan(theta[ithe]/2d0))
    n=theta[ithe]
    while delt ge 1d-5 do begin
      f=cos(n)+alog(tan(n/2d0))-ftheta
      df=-sin(n)+1d0/sin(n)
      n1=n-f/df
      delt=abs(n1-n)/n
      n=n1
    endwhile

    grav[iom,ithe]=sqrt(1d0/r^4+omega[iom]^4*r^2*sin(theta[ithe])^2 $
                   -2d0*omega[iom]^2*sin(theta[ithe])^2/r)

    corr[iom,ithe]=sqrt(tan(n)/tan(theta[ithe]))

    teff[iom,ithe]=corr[iom,ithe]*grav[iom,ithe]^0.25

  endfor

  u=where(finite(teff[iom,*]) ne 0)
  coef=poly_fit(alog(grav[iom,u]),alog(teff[iom,u]),1)
  beta[iom]=coef[1]

endfor

beta=beta[0]

end


;::::::::::::::::::::::::::::::::::::::::::::::::

pro geneve_par,par,Hfrac,oblat=oblat,makeeps=makeeps,zams=zams

;------------------------------------------------
;
; HDUST source parameters
;
;------------------------------------------------

on_error,1
compile_opt idl2

if n_elements(oblat) ne 0 then begin
  omega_interp=(1.5d0^1.5d0)*sqrt(2d0*(par-1d0)/par^3)	; Ekstrom et al. 2008, Eq. (9)
endif else omega_interp=par

beta_hdust,omega_interp,beta
ww=wrot(omega_interp)

mass=[$
14.6d0,$
12.5d0,$
10.8d0,$
9.6d0,$
8.6d0,$
7.7d0,$
6.4d0,$
5.5d0,$
4.8d0,$
4.2d0,$
3.8d0,$
3.4d0 $
]
nm=n_elements(mass)

str_mass=[$
'M14p60',$
'M12p50',$
'M10p80',$
'M9p600',$
'M8p600',$
'M7p700',$
'M6p400',$
'M5p500',$
'M4p800',$
'M4p200',$
'M3p800',$
'M3p400' $
]

st=['B0.5','B1','B1.5','B2','B2.5','B3','B4','B5','B6','B7','B8','B9']
zsun='Z01400'

vel=[$
0.6d0,$
0.7d0,$
0.8d0,$
0.9d0,$
0.95d0 $
]
nv=n_elements(vel)

str_vel=[$
'V60000',$
'V70000',$
'V80000',$
'V90000',$
'V95000' $
]

rp=dblarr(nm,nv)
omega=dblarr(nm,nv)
logL=dblarr(nm,nv)
logTeff=dblarr(nm,nv)
rp_interp=dblarr(nm)
logL_interp=dblarr(nm)
logTeff_interp=dblarr(nm)

;Hfrac=0.5d0
openw,lun1,'geneve_par.txt',/get_lun

printf,lun1,'      Omega           W               beta'
printf,lun1,[omega_interp,ww,beta]
printf,lun1

printf,lun1,'                       Mass            Rp              L'

for im=0,nm-1 do begin
  for iv=0,nv-1 do begin

    file='./stmodels/'+str_mass[im]+zsun+str_vel[iv]+'.dat'
    openr,lun,file,/get_lun
    header=strarr(2)
    data=dblarr(55,400)
    readf,lun,header,data
    close,lun
    free_lun,lun

    if n_elements(zams) eq 0 or hfrac gt 0.7 then $
      iMS=where(abs(data[21,*]-Hfrac) eq min(abs(data[21,*]-Hfrac))) $
    else begin
      iMS=0
      print,'ZAMS PARAMETERS WERE CHOSEN'
    endelse

    rp[im,iv]=data[44,iMS]
    omega[im,iv]=data[39,iMS]
    logL[im,iv]=data[3,iMS]
    logTeff[im,iv]=data[4,iMS]

  endfor

  u=reform(where(omega[im,*] eq omega_interp))
  if u[0] ne -1 then begin
    iomega=where(omega[im,*] eq omega_interp)
    rp_interp=rp[im,iomega[0]]
    logL_interp[im]=logL[im,iomega[0]]
    logL_interp[im]=logL[im,iomega[0]]
    logTeff_interp[im]=logTeff[im,iomega[0]]
    printf,lun1,'INTERPOLATION:  ',[mass[im],rp_interp,10.^logL_interp[im]]
  endif else $
  if omega_interp lt min(omega[im,*]) then begin
    printf,lun1,'OMEGA TOO SMALL FOR INTERPOLATION. PLEASE INCLUDE NEW GENEVE MODELS.'
  endif else $
  if omega_interp gt max(omega[im,*]) then $
  begin
    coef=poly_fit(omega[im,*],rp[im,*],3)
    rp_interp[im]=coef[0]+coef[1]*omega_interp+coef[2]*omega_interp^2 $
                    +coef[3]*omega_interp^3

    coef=poly_fit(omega[im,*],logL[im,*],3)
    logL_interp[im]=coef[0]+coef[1]*omega_interp+coef[2]*omega_interp^2 $
                    +coef[3]*omega_interp^3

    coef=poly_fit(omega[im,*],logTeff[im,*],3)
    logTeff_interp[im]=coef[0]+coef[1]*omega_interp+coef[2]*omega_interp^2 $
                    +coef[3]*omega_interp^3

    printf,lun1,'EXTRAPOLATION:  ',[mass[im],rp_interp[im],10.^logL_interp[im]]
  endif else $
  begin
    iomega=sort(abs(omega[im,*]-omega_interp))


; to make sure to get the two closest points AROUND omega_interp (and not at the same side)
    if (omega_interp-omega[im,iomega[0]])*(omega_interp-omega[im,iomega[1]]) lt 0 then begin
      i_closer=1
    endif else i_closer=2


    rp_interp[im]=rp[im,iomega[0]]+(omega_interp-omega[im,iomega[0]])* $
    (rp[im,iomega[1]]-rp[im,iomega[0]])/(omega[im,iomega[i_closer]]-omega[im,iomega[0]])

    logL_interp[im]=logL[im,iomega[0]]+(omega_interp-omega[im,iomega[0]])* $
    (logL[im,iomega[1]]-logL[im,iomega[0]])/(omega[im,iomega[i_closer]]-omega[im,iomega[0]])

    logTeff_interp[im]=logTeff[im,iomega[0]]+(omega_interp-omega[im,iomega[0]])* $
    (logTeff[im,iomega[1]]-logTeff[im,iomega[0]])/(omega[im,iomega[i_closer]]-omega[im,iomega[0]])

    printf,lun1,'INTERPOLATION:  ',[mass[im],rp_interp[im],10.^logL_interp[im]]
  endelse

endfor

close,lun1
free_lun,lun1

print
print,'STELLAR PARAMETERS WRITTEN AT geneve_par.txt'
print


; PLOT ------------------------------------------

if n_elements(makeeps) ne 0 then begin

set_plot,'ps'

device,filename='geneve_rp.eps',xsize=16,ysize=10,/color,/encapsulated
!p.multi=[0,4,3]
loadct,39,/silent

omtitle=textoidl('\Omega/\Omega_{crit}')
rptitle=textoidl('R_p/R_{sun}')
rptitle=textoidl('R_p [R')+sunsymbol()+']'
plotsym,8,/fill,color=250

for im=0,nm-1 do begin

  title=st[im]+'    '+string(float(mass[im]),format='(f4.1)')+' M'+sunsymbol()
  xrange=[0.5,1.1]
  yrange=[0.9*min(rp[im,*]),1.05*max(rp[im,*])]
  plot,xrange,yrange,$
       xrange=xrange,yrange=yrange,$
       xtitle=omtitle,ytitle=rptitle,$
       xstyle=1,ystyle=1,/nodata,$
       xticks=3,yticks=3,title=title

  coef=poly_fit(omega[im,*],rp[im,*],3)
  x_fit=xrange[0]+(xrange[1]-xrange[0])*findgen(100)/99d0
  y_fit=coef[0]+coef[1]*x_fit+coef[2]*x_fit^2+coef[3]*x_fit^3

  oplot,x_fit,y_fit,line=1,color=50
  oplot,omega[im,*],rp[im,*],color=250
  oplot,omega[im,*],rp[im,*],psym=6,symsize=0.4
  oplot,[omega_interp],[rp_interp[im]],psym=8,symsize=0.5

endfor
device,/close


device,filename='geneve_lum.eps',xsize=16,ysize=10,/color,/encapsulated
!p.multi=[0,4,3]
loadct,39,/silent

omtitle=textoidl('\Omega/\Omega_{crit}')
lgLtitle=textoidl('log(L/L')+sunsymbol()+')'
plotsym,8,/fill,color=250

for im=0,nm-1 do begin

  title=st[im]+'    '+string(float(mass[im]),format='(f4.1)')+' M'+sunsymbol()
  xrange=[0.5,1.1]
  yrange=[0.99*min(logL[im,*]),1.01*max(logL[im,*])]
  plot,xrange,yrange,$
       xrange=xrange,yrange=yrange,$
       xtitle=omtitle,ytitle=lgLtitle,$
       xstyle=1,ystyle=1,/nodata,$
       xticks=3,yticks=3,title=title

  coef=poly_fit(omega[im,*],logL[im,*],3)
  x_fit=xrange[0]+(xrange[1]-xrange[0])*findgen(100)/99d0
  y_fit=coef[0]+coef[1]*x_fit+coef[2]*x_fit^2+coef[3]*x_fit^3

  oplot,x_fit,y_fit,line=1,color=50
  oplot,omega[im,*],logL[im,*],color=250
  oplot,omega[im,*],logL[im,*],psym=6,symsize=0.4
  oplot,[omega_interp],[logL_interp[im]],psym=8,symsize=0.5

endfor
device,/close


set_plot,'x'
!p.multi=0
loadct,0,/silent

endif

end
