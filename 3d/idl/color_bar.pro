;###########################################################################
; AUTHOR
;
;   Yosuke Matsumoto
;   E-mail: ymatumot@ybb.ne.jp
;
;
; LICENSE
;
; This software is OSI Certified Open Source Software.
; OSI Certified is a certification mark of the Open Source Initiative.
;
; Copyright 2003 Yosuke Matsumoto
;
; This software is provided "as-is", without any express or
; implied warranty. In no event will the authors be held liable
; for any damages arising from the use of this software.
;
; Permission is granted to anyone to use this software for any
; purpose, including commercial applications, and to alter it and
; redistribute it freely, subject to the following restrictions:
;
; 1. The origin of this software must not be misrepresented; you must
;    not claim you wrote the original software. If you use this software
;    in a product, an acknowledgment in the product documentation
;    would be appreciated, but is not required.
;
; 2. Altered source versions must be plainly marked as such, and must
;    not be misrepresented as being the original software.
;
; 3. This notice may not be removed or altered from any source distribution.
;
; For more information on Open Source Software, visit the Open Source
; web site: http://www.opensource.org.
;
;###########################################################################

pro color_bar, z,ct=ct,xs,ys,xe,ye,sym=sym,maxz=maxz,minz=minz,charsize=charsize,$
               vertical=vertical
common wdw, px,py,sx,sy

;; CHECK COLOR TABLE
if(not(keyword_set(ct)))then ct=0
if(ct gt 40) then begin
    print,'No such color table.'
    return
endif

;; MAX AND MIN VALUE OF Z
if( not(keyword_set(maxz)) )then maxz = max(z)
if( not(keyword_set(minz)) )then minz = min(z)
if(keyword_set(sym) and min(z) lt 0.0) then begin
    if (maxz gt abs(minz)) then minz = -maxz
    if (abs(minz) gt maxz) then maxz = -minz
endif

;; DRAW A COLOR BAR
set_window,xs,xe,ys,ye

if(not(keyword_set(vertical)))then begin
   bar = findgen(128,2)
   bar_x = minz+findgen(128)/127.0*(maxz-minz)
   bar_y = fltarr(2)
   bar(*,1) = bar(*,0)
endif else begin
   bar = findgen(2,128)
   bar_x = fltarr(2)
   bar_y = minz+findgen(128)/127.0*(maxz-minz)
   bar(1,*) = bar(0,*)
endelse

loadct,ct,/silent
;sx = sx+0.5
;sy = sy+0.5
if(!d.name eq 'PS') then begin
   tv, bytscl(bar),px(0)+1,py(0)+1,xsize=sx,ysize=sy
endif else begin
   tv, congrid(bytscl(bar),sx, sy), px(0)+1, py(0)+1
endelse

if( not(keyword_set(charsize)) )then charsize = 1.0
loadct,12,/silent
if(not(keyword_set(vertical)))then begin
   if(0.5*abs(maxz-minz) lt 1.0e-2) then begin
      contour, bar, bar_x, bar_y, /noerase, /nodata, $
               position=[px(0),py(0),px(1),py(1)], yminor=1, yticks=1,$
               ytickname=replicate(' ',2), xticks=4, xs=1, ys=1,/device,$
               xtickformat='(e10.2)',charsize=charsize
   endif else begin
      contour, bar, bar_x, bar_y, /noerase, /nodata, $
               position=[px(0),py(0),px(1),py(1)], yminor=1, yticks=1,$
               ytickname=replicate(' ',2), xticks=4, xs=1, ys=1,/device,$
               xtickformat='(f10.2)',charsize=charsize
   endelse
endif else begin
   if(0.5*abs(maxz-minz) lt 1.0e-2) then begin
      contour, bar, bar_x, bar_y, /noerase, /nodata, $
               position=[px(0),py(0),px(1),py(1)], xminor=1, xticks=1,$
               xtickname=replicate(' ',2), ytickname=replicate(' ',2), $
               yticks=1, xs=1, ys=8,/device
      axis,yaxis=1,yrange=[min(bar_y),max(bar_y)],ytickformat='(e10.2)',$
           charsize=charsize,yticks=4,ys=1,/device
   endif else begin
      contour, bar, bar_x, bar_y, /noerase, /nodata, $
               position=[px(0),py(0),px(1),py(1)], xminor=1, xticks=1,$
               xtickname=replicate(' ',2), ytickname=replicate(' ',2), $
               yticks=1, xs=1, ys=8,/device
      axis,yaxis=1,yrange=[min(bar_y),max(bar_y)],ytickformat='(f10.2)',$
           charsize=charsize,yticks=4,ys=1,/device
   endelse
endelse
  

end

