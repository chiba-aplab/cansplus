function field_lines_2d, ux,uy,nsteps=nsteps,len=len,r0=r0,$
                      thick=thick,npos=npos,dir=dir

info = size(ux)
nx = info(1)
ny = info(2)

nmax = max([nx,ny])

if(not(keyword_set(nsteps))) then nsteps = 1e4

;INITIAL POSITION
if(not(keyword_set(r0)))then begin
    if(not(keyword_set(npos)))then begin
        npos = 100
    endif
    r0 = fltarr(2,npos)
    r0(0,*) = randomu(seed, npos) * (nx-1)
    r0(1,*) = randomu(seed, npos) * (ny-1)
endif else begin
    info=size(r0)
    npos = info(2)
endelse

r = dblarr(2,nsteps,npos)

for k=0,npos-1 do begin
    r(0,0,k) = r0(0,k)
    r(1,0,k) = r0(1,k)
endfor

;Default length
if(not(keyword_set(len)))then len = 1.0

;Default direction
if(not(keyword_set(dir)))then dir = 1.0

;;TRACE FIELD LINES
umax = max(sqrt(double(ux)^2+double(uy)^2),/NAN)
;; zscale = len/umax

for k=0,npos-1 do begin
for i=1L,nsteps-1 do begin	
    x0 = r(0,i-1,k)
    y0 = r(1,i-1,k)

    utx = interpolate(ux,x0,y0,cubic=-0.5)
    uty = interpolate(uy,x0,y0,cubic=-0.5)
    uabs = sqrt(utx^2+uty^2)
    
    if(uabs eq 0.0)then uabs = uabs+1.D-10
    r(0,i,k) = len*dir*utx/uabs + x0
    r(1,i,k) = len*dir*uty/uabs + y0
;;     r(0,i,k) = dir*utx*zscale + x0
;;     r(1,i,k) = dir*uty*zscale + y0

    flags = (r(0,i,k) gt 0) and (r(0,i,k) lt nx-1) and $
	    (r(1,i,k) gt 0) and (r(1,i,k) lt ny-1)

    if(not(flags))then begin
        istop=i
        r(0,istop+1:*,k) = r(0,istop,k)
        r(1,istop+1:*,k) = r(1,istop,k)
        break
    endif
endfor
endfor

return, r


end
