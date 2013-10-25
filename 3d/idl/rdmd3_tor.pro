narg=-1

;read switch, 0:OFF, 1:ON
filereads=1
rros=1
rprs=1
rvs=0
rbs=0
retas=0
s3d=0 ; store the data in 3D array

;fs:start file number, fe:end file number
if filereads eq 0 then begin
fs=1
fe=8
endif else begin
read,'startdata = ',fs
read,'enddata = ',fe
endelse
nx=fe-fs+1

dacgetparam,'etc/params.txt.0000','gm',gm
dacgetparam,'etc/params.txt.0000','ix',ix
dacgetparam,'etc/params.txt.0000','jx',jx
dacgetparam,'etc/params.txt.0000','kx',kx
dacgetparam,'etc/params.txt.0000','mpisize',npe
dacgetparam,'etc/params.txt.0000','margin',margin
dacgetparam,'etc/params.txt.0000','dtout',dt
dacgetparam,'etc/params.txt.0000','mpix',mpix
dacgetparam,'etc/params.txt.0000','mpiy',mpiy
dacgetparam,'etc/params.txt.0000','mpiz',mpiz
dacgetparam,'etc/params.txt.0000','nrmlt',nrmlt
;; read,'mpix=',mpix
;; read,'mpiy=',mpiy
;; read,'mpiz=',mpiz

igx=ix*mpix-(mpix-1)*margin*2
jgx=jx*mpiy-(mpiy-1)*margin*2
kgx=kx*mpiz-(mpiz-1)*margin*2
hjx=jgx/2

x=fltarr(igx)
y=fltarr(jgx)
z=fltarr(kgx)



;read x,y,z
mpirank=0
for mz=1,mpiz do begin
   for my=1,mpiy do begin
      for mx=1,mpix do begin
         xd=fltarr(ix)
         zd=fltarr(kx)
         files='etc/x.dac.'+string(mpirank,form='(i4.4)')
         dacget1d,files,xd,margin=margin
         files='etc/y.dac.'+string(mpirank,form='(i4.4)')
         dacget1d,files,yd,margin=margin
         files='etc/z.dac.'+string(mpirank,form='(i4.4)')
         dacget1d,files,zd,margin=margin
         ixs=(mx-1)*(ix-2*margin)+1-1
         ixe=ixs+ix-1
         jxs=(my-1)*(jx-2*margin)+1-1
         jxe=jxs+jx-1
         kxs=(mz-1)*(kx-2*margin)+1-1
         kxe=kxs+kx-1
         x[ixs:ixe]=xd
         y[jxs:jxe]=yd
         z[kxs:kxe]=zd
         mpirank=mpirank+1
      endfor
   endfor
endfor

read,' input z : ',temp1
o=kgx-1
for k=0,kgx-1 do begin
  if z(k) le temp1 then begin
    o=k
  endif
endfor

;read ro
data=fltarr(igx,jgx,kgx)
if rros eq 1 then begin
   ro=fltarr(igx,jgx,nx)
   if s3d eq 1 then ro_3d=fltarr(igx,jgx,kgx,nx)
   for i=fs,fe do begin
      mpirank=0
      for mz=1,mpiz do begin
         for my=1,mpiy do begin
            for mx=1,mpix do begin
               files='ro/ro.dac.'+string(i,form='(i4.4)')+'.'+string(mpirank,form='(i4.4)')
               dacget3s,files,temp,narg=narg,margin=margin
               ixs=(mx-1)*(ix-2*margin)+1-1
               ixe=ixs+ix-1
               jxs=(my-1)*(jx-2*margin)+1-1
               jxe=jxs+jx-1
               kxs=(mz-1)*(kx-2*margin)+1-1
               kxe=kxs+kx-1
               data[ixs:ixe,jxs:jxe,kxs:kxe]=temp[*,*,*]
               mpirank=mpirank+1
            endfor
         endfor
      endfor
      ro[*,*,i-fs]=data[*,*,o]
      if s3d eq 1 then ro_3d[*,*,*,i-fs]=data[*,*,*]
      print,'ro data',i
   endfor
endif

;read pr
if rprs eq 1 then begin
   pr=fltarr(igx,jgx,nx)
   if s3d eq 1 then  pr_3d=fltarr(igx,jgx,kgx,nx)
   for i=fs,fe do begin
      mpirank=0
      for mz=1,mpiz do begin
         for my=1,mpiy do begin
            for mx=1,mpix do begin
               files='pr/pr.dac.'+string(i,form='(i4.4)')+'.'+string(mpirank,form='(i4.4)')
               dacget3s,files,temp,narg=narg,margin=margin
               ixs=(mx-1)*(ix-2*margin)+1-1
               ixe=ixs+ix-1
               jxs=(my-1)*(jx-2*margin)+1-1
               jxe=jxs+jx-1
               kxs=(mz-1)*(kx-2*margin)+1-1
               kxe=kxs+kx-1
               data[ixs:ixe,jxs:jxe,kxs:kxe]=temp[*,*,*]
               mpirank=mpirank+1
            endfor
         endfor
      endfor
      pr[*,*,i-fs]=data[*,*,o]
      if s3d eq 1 then pr_3d[*,*,*,i-fs]=data[*,*,*]
      print,'pr data',i
   endfor
endif

;read vx
if rvs eq 1 then begin
   vx=fltarr(igx,jgx,nx)
   vy=fltarr(igx,jgx,nx)
   vz=fltarr(igx,jgx,nx)
   if s3d eq 1 then begin
      vx_3d=fltarr(igx,jgx,kgx,nx)
      vy_3d=fltarr(igx,jgx,kgx,nx)
      vz_3d=fltarr(igx,jgx,kgx,nx)
   endif
   for i=fs,fe do begin
      mpirank=0
      for mz=1,mpiz do begin
         for my=1,mpiy do begin
            for mx=1,mpix do begin
               files='v/vx.dac.'+string(i,form='(i4.4)')+'.'+string(mpirank,form='(i4.4)')
               dacget3s,files,temp,narg=narg,margin=margin
               ixs=(mx-1)*(ix-2*margin)+1-1
               ixe=ixs+ix-1
               jxs=(my-1)*(jx-2*margin)+1-1
               jxe=jxs+jx-1
               kxs=(mz-1)*(kx-2*margin)+1-1
               kxe=kxs+kx-1
               data[ixs:ixe,jxs:jxe,kxs:kxe]=temp[*,*,*]
               mpirank=mpirank+1
            endfor
         endfor
      endfor
      vx[*,*,i-fs]=data[*,*,o]
      if s3d eq 1 then vx_3d[*,*,*,i-fs]=data[*,*,*]
      print,'vx data',i
;read vy
      mpirank=0
      for mz=1,mpiz do begin
         for my=1,mpiy do begin
            for mx=1,mpix do begin
               files='v/vy.dac.'+string(i,form='(i4.4)')+'.'+string(mpirank,form='(i4.4)')
               dacget3s,files,temp,narg=narg,margin=margin
               ixs=(mx-1)*(ix-2*margin)+1-1
               ixe=ixs+ix-1
               jxs=(my-1)*(jx-2*margin)+1-1
               jxe=jxs+jx-1
               kxs=(mz-1)*(kx-2*margin)+1-1
               kxe=kxs+kx-1
               data[ixs:ixe,jxs:jxe,kxs:kxe]=temp[*,*,*]
               mpirank=mpirank+1
            endfor
         endfor
      endfor
      vy[*,*,i-fs]=data[*,*,o]
      if s3d eq 1 then vy_3d[*,*,*,i-fs]=data[*,*,*]
      print,'vy data',i
;read vz
      mpirank=0
      for mz=1,mpiz do begin
         for my=1,mpiy do begin
            for mx=1,mpix do begin
               files='v/vz.dac.'+string(i,form='(i4.4)')+'.'+string(mpirank,form='(i4.4)')
               dacget3s,files,temp,narg=narg,margin=margin
               ixs=(mx-1)*(ix-2*margin)+1-1
               ixe=ixs+ix-1
               jxs=(my-1)*(jx-2*margin)+1-1
               jxe=jxs+jx-1
               kxs=(mz-1)*(kx-2*margin)+1-1
               kxe=kxs+kx-1
               data[ixs:ixe,jxs:jxe,kxs:kxe]=temp[*,*,*]
               mpirank=mpirank+1
            endfor
         endfor
      endfor
      vz[*,*,i-fs]=data[*,*,o]
      if s3d eq 1 then vz_3d[*,*,*,i-fs]=data[*,*,*]
      print,'vz data',i
   endfor
endif

;read bx
if rbs eq 1 then begin
   bx=fltarr(igx,jgx,nx)
   by=fltarr(igx,jgx,nx)
   bz=fltarr(igx,jgx,nx)
   if s3d eq 1 then begin
      bx_3d=fltarr(igx,jgx,kgx,nx)
      by_3d=fltarr(igx,jgx,kgx,nx)
      bz_3d=fltarr(igx,jgx,kgx,nx)
   endif
   for i=fs,fe do begin
      mpirank=0
      for mz=1,mpiz do begin
         for my=1,mpiy do begin
            for mx=1,mpix do begin
               files='b/bx.dac.'+string(i,form='(i4.4)')+'.'+string(mpirank,form='(i4.4)')
               dacget3s,files,temp,narg=narg,margin=margin
               ixs=(mx-1)*(ix-2*margin)+1-1
               ixe=ixs+ix-1
               jxs=(my-1)*(jx-2*margin)+1-1
               jxe=jxs+jx-1
               kxs=(mz-1)*(kx-2*margin)+1-1
               kxe=kxs+kx-1
               data[ixs:ixe,jxs:jxe,kxs:kxe]=temp[*,*,*]
               mpirank=mpirank+1
            endfor
         endfor
      endfor
      bx[*,*,i-fs]=data[*,*,o]
      if s3d eq 1 then bx_3d[*,*,*,i-fs]=data[*,*,*]
      print,'bx data',i
;read by
      mpirank=0
      for mz=1,mpiz do begin
         for my=1,mpiy do begin
            for mx=1,mpix do begin
               files='b/by.dac.'+string(i,form='(i4.4)')+'.'+string(mpirank,form='(i4.4)')
               dacget3s,files,temp,narg=narg,margin=margin
               ixs=(mx-1)*(ix-2*margin)+1-1
               ixe=ixs+ix-1
               jxs=(my-1)*(jx-2*margin)+1-1
               jxe=jxs+jx-1
               kxs=(mz-1)*(kx-2*margin)+1-1
               kxe=kxs+kx-1
               data[ixs:ixe,jxs:jxe,kxs:kxe]=temp[*,*,*]
               mpirank=mpirank+1
            endfor
         endfor
      endfor
      by[*,*,i-fs]=data[*,*,o]
      if s3d eq 1 then by_3d[*,*,*,i-fs]=data[*,*,*]
      print,'by data',i
;read bz
      mpirank=0
      for mz=1,mpiz do begin
         for my=1,mpiy do begin
            for mx=1,mpix do begin
               files='b/bz.dac.'+string(i,form='(i4.4)')+'.'+string(mpirank,form='(i4.4)')
               dacget3s,files,temp,narg=narg,margin=margin
               ixs=(mx-1)*(ix-2*margin)+1-1
               ixe=ixs+ix-1
               jxs=(my-1)*(jx-2*margin)+1-1
               jxe=jxs+jx-1
               kxs=(mz-1)*(kx-2*margin)+1-1
               kxe=kxs+kx-1
               data[ixs:ixe,jxs:jxe,kxs:kxe]=temp[*,*,*]
               mpirank=mpirank+1
            endfor
         endfor
      endfor
      bz[*,*,i-fs]=data[*,*,o]
      if s3d eq 1 then bz_3d[*,*,*,i-fs]=data[*,*,*]
      print,'bz data',i
   endfor
endif

;read eta
data=fltarr(igx,jgx,kgx)
if retas eq 1 then begin
   eta=fltarr(igx,jgx,nx)
   if s3d eq 1 then eta_3d=fltarr(igx,jgx,kgx,nx)
   for i=fs,fe do begin
      mpirank=0
      for mz=1,mpiz do begin
         for my=1,mpiy do begin
            for mx=1,mpix do begin
               files='eta/eta.dac.'+string(i,form='(i4.4)')+'.'+string(mpirank,form='(i4.4)')
               dacget3s,files,temp,narg=narg,margin=margin
               ixs=(mx-1)*(ix-2*margin)+1-1
               ixe=ixs+ix-1
               jxs=(my-1)*(jx-2*margin)+1-1
               jxe=jxs+jx-1
               kxs=(mz-1)*(kx-2*margin)+1-1
               kxe=kxs+kx-1
               data[ixs:ixe,jxs:jxe,kxs:kxe]=temp[*,*,*]
               mpirank=mpirank+1
            endfor
         endfor
      endfor
      eta[*,*,i-fs]=data[*,*,o]
      if s3d eq 1 then eta_3d[*,*,*,i-fs]=data[*,*,*]
      print,'eta data',i
   endfor
endif

if rros*rprs eq 1 then begin
    te=fltarr(igx,jgx,nx)
    if s3d eq 1 then te_3d=fltarr(igx,jgx,kgx,nx)
    for i=fs,fe do begin
        te(*,*,i-fs)=pr(*,*,i-fs)*1.6d-24/(ro(*,*,i-fs)*1.38d-16)
        if s3d eq 1 then te_3d[*,*,*,i-fs]=pr_3d(*,*,*,i-fs)*1.6d-24/(ro_3d(*,*,*,i-fs)*1.38d-16)
    endfor
endif

end
