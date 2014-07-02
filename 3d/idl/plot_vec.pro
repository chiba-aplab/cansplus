;###########################################################################
; AUTHOR
;
;   Yosuke Matsumoto
;
;###########################################################################

pro plot_vec, ux, uy, facx, facy, xax=xax, yax=yax, gui=gui,color=color
common wdw, px,py,sx,sy

;; SET SIZE OF THE WINDOW
set_window,!x.window(0),!x.window(1),!y.window(0),!y.window(1)


;; CHECK INPUT VARIABLE
; if(not(keyword_set(ux)) or not(keyword_set(uy)) or keyword_set(gui))then begin
;     z = file_read(/gui)
;     if(n_elements(z) eq 1)then begin
;         if(z eq 0)then return
;     endif
; endif


;; READ SIZE OF THE MATRIX
ux = reform(ux)
uy = reform(uy)
info_ux =  size(ux)
info_uy =  size(uy)


;; CHECK
if (info_ux(0) ne info_uy(0)) then begin
    print, 'size unmatch!' 
    return
endif else begin
    n_col = info_ux(1)
    n_line = info_ux(2)
endelse


;; set x and y axis
if(not(keyword_set(xax)))then begin
    xax = findgen(n_col)
endif
if(not(keyword_set(yax)))then begin
    yax = findgen(n_line)
endif


;; shrink
n_col = n_col/facx
n_line = n_line/facy


;; draw
;;arrow_len = sqrt(max(ux)^2+max(uy)^2)
hux = congrid(ux,n_col,n_line)
huy = congrid(uy,n_col,n_line)
hx = congrid(xax,n_col)
hy = congrid(yax,n_line)

velovect,hux, huy, hx, hy, /noerase,$
  position=[px(0)/!d.x_vsize,py(0)/!d.y_vsize,px(1)/!d.x_vsize,py(1)/!d.y_vsize],$
  xrange=[min(xax),max(xax)], yrange=[min(yax),max(yax)],$
	length=2.5, xstyle=5, ystyle=5, color=color


end


