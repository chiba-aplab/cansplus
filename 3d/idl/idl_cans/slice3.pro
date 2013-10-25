pro slice3,out,dout,iout,jout,data,x,y,z,x0,x1,y0,y1,z0,z1

szx=size(x) & ix=szx(1)
szy=size(y) & jx=szy(1)
szz=size(z) & kx=szz(1)

xsml=min([x0,x1]) & xbig=max([x0,x1])
ysml=min([y0,y1]) & ybig=max([y0,y1])
zsml=min([z0,z1]) & zbig=max([z0,z1])

isml=0 & ibig=ix-1
jsml=0 & jbig=jx-1
ksml=0 & kbig=jx-1

for i=0,ix-1 do begin
  if (x(i) ge xsml) then begin
    isml=i
    goto,L101
  endif
endfor
L101 :

for i=ix-1,0,-1 do begin
  if (x(i) le xbig) then begin
    ibig=i
    goto,L111
  endif
endfor
L111 :

for j=0,jx-1 do begin
  if (y(j) ge zsml) then begin
    jsml=j
    goto,L201
  endif
endfor
L201 :

for j=jx-1,0,-1 do begin
  if (y(j) le zbig) then begin
    jbig=j
    goto,L211
  endif
endfor
L211 :

for k=0,kx-1 do begin
  if (z(k) ge zsml) then begin
    ksml=k
    goto,L301
  endif
endfor
L301 :

for k=kx-1,0,-1 do begin
  if (z(k) le zbig) then begin
    kbig=k
    goto,L311
  endif
endfor
L311 :

isz=abs(ibig-isml)+1
jsz=abs(jbig-jsml)+1
ksz=abs(kbig-ksml)+1
msz=max([isz,jsz,ksz],mmx)

out=fltarr(msz)
iout=lonarr(msz)
jout=lonarr(msz)
kout=lonarr(msz)

if (mmx eq 0) then begin
  inc=(z1-z0)/(x1-x0)
  for m=0,msz-1 do begin
    i=isml+m
    zon=inc*(x(i)-x0)+z0
  
    for j=jsml,jbig do begin
      if (y(j) ge zon) then goto,L301
    endfor
    L301 :
      zon1=y(j)   & dzon1=abs(zon1-zon)
      zon0=y(j-1) & dzon0=abs(zon0-zon)
      tmp=min([dzon0,dzon1],mx)
      if (mx eq 0) then jon=j-1 else jon=j
      out(m)=data(i,jon)
      iout(m)=i
      jout(m)=jon
  
  endfor
endif else begin
  inc=(x1-x0)/(z1-z0)
  for m=0,msz-1 do begin
    j=jsml+m
    xon=inc*(y(j)-z0)+x0
  
    for i=isml,ibig do begin
      if (x(i) ge xon) then goto,L401
    endfor
    L401 :
      xon1=x(i)   & dxon1=abs(xon1-xon)
      xon0=x(i-1) & dxon0=abs(xon0-xon)
      tmp=min([dxon0,dxon1],mx)
      if (mx eq 0) then ion=i-1 else ion=i
      out(m)=data(ion,j)
      iout(m)=ion
      jout(m)=j
  
  endfor
endelse
dout=sqrt((x(iout)-x0)^2+(y(jout)-z0)^2)
end
