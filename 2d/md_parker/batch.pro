;; SELECT IMFORMATION
print, '0:bx 1:by 2:bz 3:eta  4:phi'
print, '5:pr 6:ro 7:vx 8:vy  10:vz'

;; SIMULATION PARAMETER
beta = 2.0
gam = 5.d0/3.d0
gr = 0.077

;; READING DATA
pushd,'data'
dac_read,data,x,y,'0030_*.dac'
popd
dx = abs(x[1]-x[0])
dy = abs(y[1]-y[0])

;; DRAW 2D PROFILE WITH IMAGE, VECTOR, COLORBAR IDL NATIVE ROUTINES
img = image(alog10(data[*,*,5]),x-dx/2.,y-dy/2., $
             axis_style=2,xtickdir=1,ytickdir=1,rgb_table=33, $
             xr=[-20,20],yr=[-20,20], $
             dimension=[650,550],position=[0.15,0.1,0.8,0.9])
cb  = colorbar(target=img,orientation=1,tickdir=1,textpos=-1, $
             position=[0.85,0.1,0.9,0.9])
vct = vector(data[*,*,7],data[*,*,8],x-dx/2.,y-dy/2.,transparency=50,$
             x_subsample=10,y_subsample=10,/overplot)

;;ADJUSTING FOTMATS
img.font_name='Times'
img.font_size='16'
img.xtitle='$X $' ;;CAN USE LATEX COMMANDS FOR TITLE
img.ytitle='$Y $'
;img.scale,2.0,2.0 ;; ENLARGE IF IMAGE IS SMALL
cb.font_name='Times'
cb.font_size='14'
vct.x_subsample=4
vct.y_subsample=4*dx/dy

img.save,'output.png',resolution=200 ;; FILE FORMAT IS AUTOMATICALLY DETECTED

end
