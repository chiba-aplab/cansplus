pro dacgetparam,file,name,value,swap_endian=swap_endian

; check header
;
openr,unit,file,/get_lun,/swap_endian
tmp=FSTAT(unit) & filesize=tmp.size

str=''
m=0

while (not eof(unit)) do begin

  readf,unit,str
  if (strlen(str) ge 4) then begin
  if (strmid(str,0,4) eq '#dac') then begin
    strs=strsplit(str,' ',/extract)
    if (name eq strs[2]) then begin
      mtype=fix(strs[1])
      case mtype of
        2 : begin
          pos=strsplit(str,' ')
          value=strmid(str,pos[3])
        end
        4 : value=long(strs[3])
        5 : value=float(strs[3])
      endcase
      goto,label0
    endif
  endif
  endif

endwhile

label0:


close,unit
free_lun,unit

return
end


pro dacget1d,files,data $
   ,ix=ix,margin=margin,swap_endian=swap_endian

if n_elements(msode0) eq 0 then msode0=4
msode=msode0

mfx=n_elements(files)

if (mfx eq 1) then begin ; single file

file=files[0]

; check header
;
openr,unit,file,/get_lun,/swap_endian

L300:
unity=0L
pos=msode
point_lun,unit,pos
readu,unit,unity
if (unity eq 0) then begin
  msode=8
  goto,L300
endif
if (unity gt 1000) then endian_flag=1 else endian_flag=0
; if (endian_flag = 1) then machine-arch and file-endian is different.

mver=0L
pos=pos+4+msode*2
point_lun,unit,pos
readu,unit,mver
if endian_flag then byteorder,mver,/ntohl

;   for data type: integer=4, real=5, double=6
mtype=0L
pos=pos+4+msode*2
point_lun,unit,pos
readu,unit,mtype
if endian_flag then byteorder,mtype,/ntohl

mndim=0L
pos=pos+4+msode*2
point_lun,unit,pos
readu,unit,mndim
if endian_flag then byteorder,mndim,/ntohl

mdim=lonarr(mndim)
pos=pos+4+msode*2
point_lun,unit,pos
readu,unit,mdim
if endian_flag then byteorder,mdim,/ntohl

ix=mdim[mndim-1]

case mtype of
  4: begin
    data=lonarr(ix)
  end
  5: begin
    data=fltarr(ix)
  end
  6: begin
    data=dblarr(ix)
  end
endcase

headersize=(4L+msode*2)*4+(4L*mndim+msode*2)
tmp=FSTAT(unit) & filesize=tmp.size

pos=headersize      ; position in file

   if (pos lt filesize) then begin
     point_lun,unit,pos+msode
     readu,unit,data
   endif else begin
     goto,label0
   endelse

label0:

close,unit
free_lun,unit

if endian_flag then begin
case mtype of
  4: byteorder,data,/ntohl
  5: byteorder,data,/xdrtof
  6: byteorder,data,/xdrtod
endcase
endif

endif else begin ; multi files

for mf=0,mfx-1 do begin

file=files[mf]

dacget1d,file,data0,/swap_endian

if (mf eq 0) then begin
  sz=size(data0)
  ix=sz[1]
  igx=ix*mfx-2*margin*(mfx-1)
  mtype=sz[2]+1
  case mtype of
    4: begin
      data=lonarr(igx)
    end
    5: begin
      data=fltarr(igx)
    end
    6: begin
      data=dblarr(igx)
    end
  endcase
endif

moff=mf*(ix-2*margin)+margin
data[moff:moff+ix-1-margin]=data0[margin:ix-1]

if (margin ge 1) then $
  if (mf eq 0) then data[0:margin-1]=data0[0:margin-1]

endfor


endelse

ix=n_elements(data)


return
end


pro dacget3s,files,data,narg=narg $
   ,ix=ix,jx=jx,kx=kx,nx=nx,narr=narr,ndx=ndx,margin=margin,swap_endian=swap_endian


if n_elements(narg) eq 0 then narg=-1
if n_elements(msode0) eq 0 then msode0=4
msode=msode0

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
openr,unit,file,/get_lun,/swap_endian

L300:
unity=0L
pos=msode
point_lun,unit,pos
readu,unit,unity
if (unity eq 0) then begin
  msode=8
  goto,L300
endif
if (unity gt 1000) then endian_flag=1 else endian_flag=0
; if (endian_flag = 1) then machine-arch and file-endian is different.

mver=0L
pos=pos+4+msode*2
point_lun,unit,pos
readu,unit,mver
if endian_flag then byteorder,mver,/ntohl

;   for data type: integer=4, real=5, double=6
mtype=0L
pos=pos+4+msode*2
point_lun,unit,pos
readu,unit,mtype
if endian_flag then byteorder,mtype,/ntohl

mndim=0L
pos=pos+4+msode*2
point_lun,unit,pos
readu,unit,mndim
if endian_flag then byteorder,mndim,/ntohl

mdim=lonarr(mndim)
pos=pos+4+msode*2
point_lun,unit,pos
readu,unit,mdim
if endian_flag then byteorder,mdim,/ntohl

ix=mdim[0] & jx=mdim[1] & kx=mdim[2]
if (mdim[mndim-1] eq -1) then begin
headersize=(4L+msode*2)*4+(4L*mndim+msode*2)
tmp=FSTAT(unit) & filesize=tmp.size
datasize=filesize-headersize
if (mtype eq 6) then dataunitsize=8L*ix*jx*kx else dataunitsize=4L*ix*jx*kx
dataperstep=dataunitsize+msode*2
ndx=datasize/dataperstep
endif else begin
  ndx=mdim[mndim-1]
endelse

if (allreadflag eq 1) then narr=indgen(ndx)

sz=size(narr)
ndmax=sz[1]

case mtype of
  4: begin
    data0=lonarr(ix,jx,kx)
    data=lonarr(ix,jx,kx,ndmax)
  end
  5: begin
    data0=fltarr(ix,jx,kx)
    data=fltarr(ix,jx,kx,ndmax)
  end
  6: begin
    data0=dblarr(ix,jx,kx)
    data=dblarr(ix,jx,kx,ndmax)
  end
endcase

for nd=0,ndmax-1 do begin
   nn=narr[nd]
   pos=headersize+nn*dataperstep      ; position in file

   if (pos lt filesize) then begin
     point_lun,unit,pos+msode
     readu,unit,data0
     data[*,*,*,nd]=data0
   endif else begin
     goto,label0
   endelse

endfor
label0:

ndx=nd
data=data[*,*,*,0:ndx-1]
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

endif else begin

for mf=0,mfx-1 do begin

file=files[mf]

dacget3s,file,data0,narg=narg $
   ,ix=ix,jx=jx,kx=kx,nx=nx,narr=narr,ndx=ndx,/swap_endian

if (mf eq 0) then begin
  igx=ix
  jgx=jx
  kgx=kx*mfx-2*margin*(mfx-1)
  sz=size(data0)
  mtype=sz[sz[0]+1]+1
  case mtype of
    4: begin
      data=lonarr(ix,jx,kgx,ndx)
    end
    5: begin
      data=fltarr(ix,jx,kgx,ndx)
    end
    6: begin
      data=dblarr(ix,jx,kgx,ndx)
    end
  endcase
endif

moff=mf*(kx-2*margin)+margin
data[*,*,moff:moff+kx-1-margin,*]=data0[*,*,margin:kx-1,*]

if (margin ge 1) then $
if (mf eq 0) then data[*,*,0:margin-1,*]=data0[*,*,0:margin-1,*]

endfor

ix=igx
jx=jgx
kx=kgx

endelse

return
end


pro dac_read, data, x, y, z, filename, swap_endian=swap_endian

if(not(keyword_set(filename)))then retall

flist = file_search(filename,count=count)
if(count eq 0) then begin
    print,'No such file'
    retall
endif

dir = file_dirname(flist[0],/mark_directory)
dacgetparam,dir+'params_rank=0000.txt','gm',gm,/swap_endian
dacgetparam,dir+'params_rank=0000.txt','ix',ix,/swap_endian
dacgetparam,dir+'params_rank=0000.txt','jx',jx,/swap_endian
dacgetparam,dir+'params_rank=0000.txt','kx',kx,/swap_endian
dacgetparam,dir+'params_rank=0000.txt','mpisize',npe,/swap_endian
dacgetparam,dir+'params_rank=0000.txt','margin',margin,/swap_endian
dacgetparam,dir+'params_rank=0000.txt','dtout',dt,/swap_endian
dacgetparam,dir+'params_rank=0000.txt','mpix',mpix,/swap_endian
dacgetparam,dir+'params_rank=0000.txt','mpiy',mpiy,/swap_endian
dacgetparam,dir+'params_rank=0000.txt','mpiz',mpiz,/swap_endian

igx=ix*mpix-(mpix-1)*margin*2
jgx=jx*mpiy-(mpiy-1)*margin*2
kgx=kx*mpiz-(mpiz-1)*margin*2

x=fltarr(igx)
y=fltarr(jgx)
z=fltarr(kgx)

;read x,y,z
mpirank=0
for mz=1,mpiz do begin
   for my=1,mpiy do begin
      for mx=1,mpix do begin
         files=dir+'x_rank='+string(mpirank,form='(i4.4)')+'.dac'
         dacget1d,files,xd,margin=margin,/swap_endian
         files=dir+'y_rank='+string(mpirank,form='(i4.4)')+'.dac'
         dacget1d,files,yd,margin=margin,/swap_endian
         files=dir+'z_rank='+string(mpirank,form='(i4.4)')+'.dac'
         dacget1d,files,zd,margin=margin,/swap_endian
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

if( (count mod mpirank) ne 0)then begin
  print, 'Number of files and procs. is mismatch.'
  retall
endif
nproc = mpirank
nd = count/nproc
data=fltarr(igx,jgx,kgx,nd)

for l=0L,nd-1 do begin
   mpirank=0
   for mz=1,mpiz do begin
      for my=1,mpiy do begin
         for mx=1,mpix do begin
            print, 'Reading...',flist[l*nproc+mpirank]
            dacget3s,flist[l*nproc+mpirank],temp,margin=margin,/swap_endian
            ixs=(mx-1)*(ix-2*margin)+1-1
            ixe=ixs+ix-1
            jxs=(my-1)*(jx-2*margin)+1-1
            jxe=jxs+jx-1
            kxs=(mz-1)*(kx-2*margin)+1-1
            kxe=kxs+kx-1
            data[ixs:ixe,jxs:jxe,kxs:kxe,l]=temp[*,*,*]
            mpirank=mpirank+1
         endfor
      endfor
   endfor
endfor

data = reform(data[margin:igx-1-margin,margin:jgx-1-margin,margin:kgx-1-margin,0:nd-1])
x = reform(x[margin:igx-1-margin])
y = reform(y[margin:jgx-1-margin])
z = reform(z[margin:kgx-1-margin])

end
