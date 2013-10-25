pro dacget4s,files,data,narg=narg $
   ,ix=ix,jx=jx,kx=kx,mx,nx=nx,narr=narr,ndx=ndx,margin=margin $
   ,pardir=pardir,ipear=ipear,jpear=jpear,kpear=kpear,mpear=mpear


if n_elements(narg) eq 0 then narg=-1

sz=size(narg)
allreadflag=0
if (sz[0] ge 1) then  begin
  narr=narg
endif else begin
  if (narg ge 0) then begin
    narr=[narg]
  endif else begin
    allreadflag=1
  endelse
endelse

mfx=n_elements(files)

if (mfx eq 1) then begin

file=files[0]

; check header
;
openr,unit,file,/get_lun

unity=0L
pos=4
point_lun,unit,pos
readu,unit,unity
if (unity gt 1000) then endian_flag=1 else endian_flag=0
; if (endian_flag = 1) then machine-arch and file-endian is different.

mver=0L
pos=pos+12
point_lun,unit,pos
readu,unit,mver
if endian_flag then byteorder,mver,/ntohl

;   for data type: integer=4, real=5, double=6
mtype=0L
pos=pos+12
point_lun,unit,pos
readu,unit,mtype
if endian_flag then byteorder,mtype,/ntohl

mndim=0L
pos=pos+12
point_lun,unit,pos
readu,unit,mndim
if endian_flag then byteorder,mndim,/ntohl

mdim=lonarr(mndim)
pos=pos+12
point_lun,unit,pos
readu,unit,mdim
if endian_flag then byteorder,mdim,/ntohl

ix=mdim[0] & jx=mdim[1] & kx=mdim[2] & mx=mdim[3]
if (mdim[mndim-1] eq -1) then begin
headersize=12L*4+mndim*4+8
tmp=FSTAT(unit) & filesize=tmp.size
datasize=filesize-headersize
if (mtype eq 6) then dataunitsize=8L*ix*jx*kx*mx $
                else dataunitsize=4L*ix*jx*kx*mx
dataperstep=dataunitsize+8
ndx=datasize/dataperstep
endif else begin
  ndx=mdim[mndim-1]
endelse

if (allreadflag eq 1) then narr=indgen(ndx)

sz=size(narr)
ndmax=sz[1]

case mtype of
  4: begin
    data0=lonarr(ix,jx,kx,mx)
    data=lonarr(ix,jx,kx,mx,ndmax)
  end
  5: begin
    data0=fltarr(ix,jx,kx,mx)
    data=fltarr(ix,jx,kx,mx,ndmax)
  end
  6: begin
    data0=dblarr(ix,jx,kx,mx)
    data=dblarr(ix,jx,kx,mx,ndmax)
  end
endcase

for nd=0,ndmax-1 do begin
   nn=narr[nd]
   pos=headersize+nn*dataperstep      ; position in file

   if (pos lt filesize) then begin
     point_lun,unit,pos+4
     readu,unit,data0
     data[*,*,*,*,nd]=data0
   endif else begin
     goto,label0
   endelse

endfor
label0:

ndx=nd
data=data[*,*,*,*,0:ndx-1]
nx=ndx

close,unit
free_lun,unit

if endian_flag then begin
case mtype of
  4: byteorder,data,/ntohl
  5: byteorder,data,/xdrtof
  6: byteorder,data,/xdrtod
endcase
endif

endif else begin ; (mfx eq 1)

for mf=0,mfx-1 do begin

file=files[mf]

dacget4s,file,data0,narg=narg $
   ,ix=ix,jx=jx,kx=kx,mx=mx,nx=nx,narr=narr,ndx=ndx

if (n_elements(pardir) ne 0) then begin
if (mf eq 0) then begin

case pardir of
 0: begin
  mgx=mx
  kgx=kx
  jgx=jx
  igx=ix*mfx-2*margin*(mfx-1)
  sz=size(data0)
  mtype=sz[sz[0]+1]+1
  case mtype of
    4: begin
      data=lonarr(igx,jgx,kgx,mgx,ndx)
    end
    5: begin
      data=fltarr(igx,jgx,kgx,mgx,ndx)
    end
    6: begin
      data=dblarr(igx,jgx,kgx,mgx,ndx)
    end
  endcase
 end
 1: begin
  igx=ix
  jgx=jx*mfx-2*margin*(mfx-1)
  kgx=kx
  mgx=mx
  sz=size(data0)
  mtype=sz[sz[0]+1]+1
  case mtype of
    4: begin
      data=lonarr(igx,jgx,kgx,mgx,ndx)
    end
    5: begin
      data=fltarr(igx,jgx,kgx,mgx,ndx)
    end
    6: begin
      data=dblarr(igx,jgx,kgx,mgx,ndx)
    end
  endcase
 end
 2: begin
  igx=ix
  jgx=jx
  kgx=kx*mfx-2*margin*(mfx-1)
  mgx=mx
  sz=size(data0)
  mtype=sz[sz[0]+1]+1
  case mtype of
    4: begin
      data=lonarr(igx,jgx,kgx,mgx,ndx)
    end
    5: begin
      data=fltarr(igx,jgx,kgx,mgx,ndx)
    end
    6: begin
      data=dblarr(igx,jgx,kgx,mgx,ndx)
    end
  endcase
 end
endcase
endif ; (mf eq 0)

case pardir of
 0: begin
   moff=mf*(ix-2*margin)+margin
   data[moff:moff+ix-1-margin,*,*,*,*]=data0[margin:ix-1,*,*,*,*]
   if (margin ge 1) then $
   if (mf eq 0) then data[0:margin-1,*,*,*,*]=data0[0:margin-1,*,*,*,*]
 end
 1: begin
   moff=mf*(jx-2*margin)+margin
   data[*,moff:moff+jx-1-margin,*,*,*]=data0[*,margin:jx-1,*,*,*]
   if (margin ge 1) then $
   if (mf eq 0) then data[*,0:margin-1,*,*,*]=data0[*,0:margin-1,*,*,*]
 end
 2: begin
   moff=mf*(kx-2*margin)+margin
   data[*,*,moff:moff+kx-1-margin,*,*]=data0[*,*,margin:kx-1,*,*]
   if (margin ge 1) then $
   if (mf eq 0) then data[*,*,0:margin-1,*,*]=data0[*,*,0:margin-1,*,*]
 end
endcase


ix=igx
jx=jgx
kx=kgx
mx=mgx

endif ; (n_elements(pardir) ne 0) 

if (mf eq 0) then begin

  ipex=max(ipear)+1
  jpex=max(jpear)+1
  kpex=max(kpear)+1
  mpex=max(mpear)+1

  igx=ix*ipex-2*margin*(ipex-1)
  jgx=jx*jpex-2*margin*(jpex-1)
  kgx=kx*kpex-2*margin*(kpex-1)
  mgx=mx*mpex-2*margin*(mpex-1)

  sz=size(data0)
  mtype=sz[sz[0]+1]+1

  case mtype of
    4: begin
      data=lonarr(igx,jgx,kgx,mgx,nx)
    end
    5: begin
      data=fltarr(igx,jgx,kgx,mgx,nx)
    end
    6: begin
      data=dblarr(igx,jgx,kgx,mgx,nx)
    end
  endcase

endif ; (mf eq 0)

ipe=ipear[mf]
jpe=jpear[mf]
kpe=kpear[mf]
mpe=mpear[mf]

ig0=ipe*(ix-2*margin)
jg0=jpe*(jx-2*margin)
kg0=kpe*(kx-2*margin)
mg0=mpe*(mx-2*margin)
ig1=ig0+ix-1
jg1=jg0+jx-1
kg1=kg0+kx-1
mg1=mg0+mx-1

data[ig0:ig1,jg0:jg1,kg0:kg1,mg0:mg1,*]=data0

ix=igx
jx=jgx
kx=kgx
mx=mgx

endfor ; mf=0,mfx-1

endelse ; (mfx eq 1)



return
end
