;###########################################################################
; AUTHOR
;
;   Yosuke Matsumoto
;
;###########################################################################

pro plot_cnt, z,xax=xax,yax=yax,lvl,maxd=maxd,mind=mind,_Extra=_Extra
common wdw, px,py,sx,sy


if(not(keyword_set(lvl)))then lvl = 15

;; READ THE SIZE OF THE MATRIX
info_z = size(reform(z))
n_col = info_z(1)
n_line = info_z(2)

;; set the size of the window
set_window,!x.window(0),!x.window(1),!y.window(0),!y.window(1)

;; draw
if(not(keyword_set(xax)))then xax = findgen(n_col)
if(not(keyword_set(yax)))then yax = findgen(n_line)
if(not(keyword_set(maxd)))then maxd = max(z)
if(not(keyword_set(mind)))then mind = min(z)

contour, reform(z), xax, yax, /noerase, /follow, xstyle=5, ystyle=5,$
         position=[px(0),py(0),px(1),py(1)], nlevels=lvl, $
         max_value=maxd, min_value=mind, /device,         $
         _Extra=_Extra

end



