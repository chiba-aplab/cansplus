;###########################################################################
; AUTHOR
;
;   Yosuke Matsumoto
;
;###########################################################################


pro plot_clcnt, z,xax=xax,yax=yax,ct=ct,sym=sym,maxz=maxz,minz=minz,$
                keep_aspect=keep_aspect,pos=pos,nocolbar=nocolbar,$
                noerase=noerase,_Extra=_Extra,charsize=charsize,$
                ver_colbar=ver_colbar

common wdw, px,py,sx,sy

;; CHECK INPUT VARIABLE
if(not(keyword_set(z)))then begin
    z = file_read()
    if(n_elements(z) eq 1)then begin
        if(z eq 0)then return
    endif
endif

;; CHECK COLOR TABLE
if(not(keyword_set(ct)))then ct=0
if(ct gt 40) then begin
    print,'No such color table.'
    return
endif

;; READ THE SIZE OF THE MATRIX
z = reform(z)
info_z = size(z)
n_col = info_z(1)
n_line = info_z(2)
if(keyword_set(xax) and keyword_set(yax))then begin
   xsz = (max(xax)-min(xax))/(max(xax)+max(yax)-min(xax)-min(yax))*1280
   ysz = (max(yax)-min(yax))/(max(xax)+max(yax)-min(xax)-min(yax))*1280
endif else begin
   xsz = float(n_col)/float(n_col+n_line)*1280
   ysz = float(n_line)/float(n_col+n_line)*1280
endelse
if(!d.window eq -1 and !d.name eq 'X')then window
if(!d.name ne 'PS'  and !d.name ne 'Z' and not keyword_set(keep_aspect))then $
   window,!d.window,xsize=xsz,ysize=ysz


;; SET X AND Y AXIS IF NOT
if(not(keyword_set(xax)))then xax = findgen(n_col)
if(not(keyword_set(yax)))then yax = findgen(n_line)
;; ADJUST X and Y AXIS TO PIXELS
xax2 = fltarr(n_col+1)
xax2(0) = 2.*xax(0)-xax(1)
xax2(1:n_col) = xax
yax2 = fltarr(n_line+1)
yax2(0) = 2.*yax(0)-yax(1)
yax2(1:n_line) = yax
xax2=xax2+(max(xax2)-min(xax2))/n_col*0.5
yax2=yax2+(max(yax2)-min(yax2))/n_line*0.5

if(not(keyword_set(noerase)))then erase

;; Drawing area
if(not(keyword_set(ver_colbar)))then begin
   xs = 0.15
   ys = 0.18
   xe = 0.92
   ye = 0.95
   ;; for color bar
   xsc = 0.15
   ysc = 0.06
   xec = 0.92
   yec = 0.08
endif else begin
   xs = 0.15
   ys = 0.15
   xe = 0.87
   ye = 0.87
   ;; for color bar
   xsc = 0.89
   ysc = 0.15
   xec = 0.91
   yec = 0.87
endelse
  
if(keyword_set(pos) and n_elements(pos) eq 4)then begin
   xs=pos(0)
   ys=pos(1)
   xe=pos(2)
   ye=pos(3)
endif

;; ADD COLOR BAR
if(not(keyword_set(nocolbar)) and not(keyword_set(pos)))then begin

   if(not(keyword_set(ver_colbar)))then begin
      if(keyword_set(sym) and min(z) lt 0.0) then begin
         color_bar,z,ct=ct,xsc,ysc,xec,yec,maxz=maxz,minz=minz,/sym,charsize=charsize
      endif else begin
         color_bar,z,ct=ct,xsc,ysc,xec,yec,maxz=maxz,minz=minz,charsize=charsize
      endelse
   endif else begin
      if(keyword_set(sym) and min(z) lt 0.0) then begin
         color_bar,z,ct=ct,xsc,ysc,xec,yec,maxz=maxz,minz=minz,/sym,charsize=charsize,/ver
      endif else begin
         color_bar,z,ct=ct,xsc,ysc,xec,yec,maxz=maxz,minz=minz,charsize=charsize,/ver
      endelse
   endelse

endif

set_window, xs, xe, ys, ye

;; DRAW
if( not(keyword_set(maxz)) )then maxz = max(z)
if( not(keyword_set(minz)) )then minz = min(z)
if (keyword_set(sym) and min(z) lt 0.0) then begin
    if (maxz gt abs(minz)) then minz = -maxz
    if (abs(minz) gt maxz) then maxz = -minz
endif

loadct, ct, /silent
if(!d.name eq 'PS') then begin
    tv,bytscl(z,max=maxz,min=minz),px(0)+1,py(0)+1,xsize=sx,ysize=sy
endif else begin
    tv, congrid(bytscl(z,max=maxz, min=minz),sx,sy), px(0)+1, py(0)+1
endelse

if( not(keyword_set(charsize)) )then charsize = 1.0
loadct,12,/silent
contour,congrid(z,n_col+1,n_line+1),xax2,yax2,$
  /nodata,position=[px(0),py(0),px(1),py(1)],$
  xs=1,ys=1,/noerase,/device,xticklen=-0.02,yticklen=-0.02,charsize=charsize,$
  _Extra=_Extra


;; RESET
loadct, 12, /silent


end




