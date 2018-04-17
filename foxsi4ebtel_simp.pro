pro foxsi4ebtel_simp

  ; Code to take an EBTEL sim DEM and fold through FOXSI response
  ;   Case 1:  Base model
  ;   Case 2:  10 times more frequent nanoflares with 1/10 the energy
  ;   Case 3:  2 times more frequent nanoflares with 1/2 the energy
  ;   Case 4:  10 times shorter duration nanoflares with 10 times greater peak heating rate
  ;
  ; Need foxsi-smex stuff from github in your IDL path
  ;
  ;-------------------------------------------------
  ;
  ; 16-Aug-2016 IGH   New EBTEL sims
  ; 17-Aug-2016 IGH   Add in some RHESSI
  ; 16-Sep-2016 IGH   "Final" version using sensitivity limit of 2\sigma=2*sqrt(background)
  ;
  ; 23-Sep-2016 IGH   New simplier version based on foxsi_resp_ebtel4_hsi.pro
  ;    Sent to Amir et al. who then generated the final fig (include SXR telescope)
  ;    Looks like different EBTEL were added/used for the final fig and not just the ones here???
  ;
  ; 12-Apr-2018 IGH   Updated to use new background and response file
  ;    No 0.17% in the background now?
  ;    Plotting 2\sigma background threshold not the background (as before)
  ;
  ;-------------------------------------------------
  ; Jim says use "pixel area" in cm^2 so multiple by the size of one FOXSI pixel?
  ; Surely emission would be bigger than that?
  ; 9'x9' fov but how many pixels ?
  ; Just use the 8" FWHM instead, then?
  aras=8.0*8.0
  arcm=aras*7.25d7^2

  ;-------------------------------------------------
  ; Case 1
  restore,file='case1.sav';,/ver
  ; Just using the coronal DEM
  dem1=dem_cor_avg
  logt=logtdem
  dlogt=logt[1]-logt[0]
  dt=10d^(logt+0.5*dlogt)-10d^(logt-0.5*dlogt)
  ;-------------------------------------------------
  ; Case 2
  restore,file='case2.sav';,/ver
  dem2=dem_cor_avg
  ;-------------------------------------------------
  ; Case 3
  restore,file='case3.sav';,/ver
  dem3=dem_cor_avg
  ;-------------------------------------------------
  ; Case 2
  restore,file='case4.sav';,/ver
  dem4=dem_cor_avg
  ;-------------------------------------------------
  ; To get from cm^-5 K^-1 to cm^-3
  em1=dem1*dt*arcm
  em2=dem2*dt*arcm
  em3=dem3*dt*arcm
  em4=dem4*dt*arcm
  ;-------------------------------------------------
  gd=where(logt gt 6.0,ngd)
  logtx=logt[gd]
  mk2kev=0.08617
  tkevx=10d^logtx*1d-6*mk2kev

  emx1=em1[gd]
  emx2=em2[gd]
  emx3=em3[gd]
  emx4=em4[gd]

  ;-------------------------------------------------
  ; Energy binning of the FOXSI spectrum
  mine=5
  maxe=25
  ebin=1/3.
  ns=(maxe-mine)/ebin

  edges=dblarr(2,ns)
  edges[0,*]=mine+findgen(ns)*ebin
  edges[1,*]=mine+(1.+findgen(ns))*ebin
  de=edges[1,0]-edges[0,0]
  engs=get_edges(edges,/mean)
  
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Load in the foxsi stuff (bloody common blocks!)
  ; This assumes the foxsi-smex stuff is in your IDL path
  @foxsi-smex-setup-script

  ; FOXSI background
  ; One factor of two as two telescopes
  ; And another factor of two for higher background due to orbit (Steve's email)
  ; And using the 0.17% HPD Albert suggested
  ; Old version
  ;  fxb=2*2*engs^(-0.8) * exp(-engs/30.)*ebin*0.17/100.

  ; Albert's new version
  ; fxb is in count/s

  fxb=foxsi_get_instrument_background(energy_arr=engs,/in_hpd)*ebin
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  phfx1=dblarr(ngd,ns)
  phfx2=dblarr(ngd,ns)
  phfx3=dblarr(ngd,ns)
  phfx4=dblarr(ngd,ns)

  for i=0,ngd-1 do begin
    ; f_vth.pro only works >1MK
    ; Not an issue for us as FOXSI is only above 5 keV
    ; Probably needs actual CHIANTI thermal spectrum for SXR instrument
    phfx1[i,*]=f_vth( edges, [emx1[i]*1d-49, tkevx[i], 1.] )
    phfx2[i,*]=f_vth( edges, [emx2[i]*1d-49, tkevx[i], 1.] )
    phfx3[i,*]=f_vth( edges, [emx3[i]*1d-49, tkevx[i], 1.] )
    phfx4[i,*]=f_vth( edges, [emx4[i]*1d-49, tkevx[i], 1.] )
  endfor

  phfxt1=total(phfx1,1,/nan)
  phfxt2=total(phfx2,1,/nan)
  phfxt3=total(phfx3,1,/nan)
  phfxt4=total(phfx4,1,/nan)

  ;---------------------------------
  ; DEM total eng is this want the non-thermal power to match it
  dem_th_eng=3.75d10*arcm
  ; Then divide by how often repeat
  nn_pow=dem_th_eng/3000.

  ; Just using some typical non-thermal values
  del=7
  ec=10.
  gam=del-1
  fe=nn_pow*(gam-1)/(gam*ec*1.6d-9)
  thck_p=[fe*1d-35,del,3000,3,ec,3000.]
  phnn=f_thick2(edges,thck_p)

  del=7
  ec2=5.
  gam=del-1
  fe=nn_pow*(gam-1)/(gam*ec2*1.6d-9)
  thck_p=[fe*1d-35,del,3000,3,ec2,3000.]
  phnn2=f_thick2(edges,thck_p)

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  foxsi_root_path = '~/github/foxsi-smex/idl'
  foxsi_data_path = strmid(foxsi_root_path, 0, strlen(foxsi_root_path)-3) + 'data/'
  fea = foxsi_get_effective_area( energy_arr=engs)

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  fxr1 = fea.eff_area_cm2 * phfxt1 * ebin
  fxr2 = fea.eff_area_cm2 * phfxt2 * ebin
  fxr3 = fea.eff_area_cm2 * phfxt3 * ebin
  fxr4 = fea.eff_area_cm2 * phfxt4 * ebin
  fxrn = fea.eff_area_cm2 * phnn * ebin
  fxrn2 = fea.eff_area_cm2 * phnn2 * ebin

  ; Doing this for an observation of dur
  dur=60.*60.
  fxc1=fxr1*dur
  fxc2=fxr2*dur
  fxc3=fxr3*dur
  fxc4=fxr4*dur
  fxcn=fxrn*dur
  fxcn2=fxrn2*dur

  !p.multi=[0,3,1]
  !x.style=17
  !y.style=17
  linecolors
  !p.thick=2
  !p.charsize=3

  ;Plot the EBTEL DEMs
  plot,logtx,emx1,xtitle='Log!D10!NT',ytitle='EBTEL Emission Measure [1/cm3]',$
    /ylog,xrange=[6.35,7.25],yrange=[2d41,3d45],/nodata
  oplot,logtx,emx1,color=2
  oplot,logtx,emx2,color=4
  oplot,logtx,emx3,color=7
  oplot,logtx,emx4,color=5,lines=2

  ; Plot the X-ray spectrum
  plot,engs,phfxt1,xtitle='Energy [keV]', ytitle='X-ray Spectrum [ph/s/keV/cm2]',/nodata,$
    /ylog,/xlog,xrange=[5,25],yrange=[2e-7,1e1]

  oplot,engs,phfxt1,color=2
  oplot,engs,phfxt2,color=4
  oplot,engs,phfxt3,color=7
  oplot,engs,phfxt4,color=5
  oplot,engs,phnn,color=9
  oplot,engs,phnn2,color=9

  ; Plot the FOXSI Spectrum
  plot,engs,fxc1,xtitle='Energy [keV]',ytitle='FOXSI Spectrum [count]',$
    /nodata,/ylog,xrange=[5,19],yrange=[1e0,2e4]
  oplot,engs,fxc1,color=2,psym=10
  oplot,engs,fxc2,color=4,psym=10
  oplot,engs,fxc3,color=7,psym=10
  oplot,engs,fxc4,color=5,psym=10
  oplot,engs,fxcn,color=9,psym=10
  oplot,engs,fxcn2,color=9,psym=10

  loadct,0,/silent
  oplot,engs,3*sqrt(fxb*dur),color=150,psym=10

  stop
end
