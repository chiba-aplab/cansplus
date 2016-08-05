;; SIMULATION PARAMETER
beta = 2.0
gam = 5.d0/3.d0
gr = 0.077

;; READING DATA
pushd,'data'
dac_read,data,x,y,'0013*.dac'
popd
dx = abs(x[1]-x[0])
dy = abs(y[1]-y[0])

;; DRAW 2D PROFILE WITH IMAGE, VECTOR, COLORBAR IDL NATIVE ROUTINES
img = image(data[*,*,5],x-dx/2.,y-dy/2.,axis_style=2,xtickdir=1,ytickdir=1,rgb_table=33)
cb  = colorbar(target=img,orientation=1,tickdir=1,position=[0.85,0.1,0.9,0.9])
vct = vector(data[*,*,7],data[*,*,8],x-dx/2.,y-dy/2.,transparency=50,$
             x_subsample=10,y_subsample=10,/overplot)

;;ADJUSTING FOTMATS
img.font_name='Times'
img.font_size='16'
img.xtitle='$X (\lambda)$' ;;CAN USE LATEX COMMANDS FOR TITLE
img.ytitle='$Y (\lambda)$'
img.scale,2.0,2.0 ;; ENLARGE IF IMAGE IS SMALL
cb.font_name='Times'
cb.font_size='14'
vct.y_subsample=7

img.save,'output.png',resolution=200 ;; FILE FORMAT IS AUTOMATICALLY DETECTED


;; PLOT LINEAR GROWTH OF THE FASTEST GROWING MODE
pushd,'data'
dac_read,vy,x,y,'*_vy*.dac'
popd

info=size(vy)
nx=info[1]
ny=info[2]
nt=info[3]

time=findgen(nt)*5.0*sqrt(1.+beta*gam/2)
tmp = fltarr(nx,ny)
fgm = fltarr(nt)

for l=0,nt-1 do begin
 for j=0,ny-1 do begin
   tmp[*,j] = fft(vy[*,j,l],-1)
 endfor

 fgm[l] = total(abs(tmp[1,*]),2)/ny

endfor

pl1 = plot(time[1:*],fgm[1:*])
pl2 = plot(time[1:*],2d-4*exp(gr*time[1:*]),/over)
label = text(0.2,0.8,'$exp(\gamma_{fgm} t)$',/normal)

;;FOR ADJUSTING FOTMATS
pl1.xstyle=1
pl1.ystyle=1
pl1.font_size=16
pl1.font_name='Time'
pl1.ylog=1
pl1.xtitle='$time (\lambda/V_0)$
pl1.ytitle='$V_{y}/V_0 (m=1)$
pl1.yrange=[2d-4,0.02]
pl1.thick=2
pl2.linestyle=2
pl2.thick=2
pl1.scale,0.9
pl1.color='red'

pl1.save,'kh_linear_growth.png',resolution=200

end
