;; SIMULATION PARAMETER

;; READING DATA
pushd,'data'
dac_read,data,x,y,'0026*.dac'
popd
info=size(data)
nx = info[1]
ny = info[2]
np = info[3]
dx = abs(x[1]-x[0])
dy = abs(y[1]-y[0])

;; COPY ORIGINAL DATA TO Y<0 REGION
data2 = findgen(nx,ny*2,np)
data2[*,ny:2*ny-1,*]=data
data2[*,0:ny-1,*]=reverse(data,2)

data2[*,0:ny-1,0] = -data2[*,0:ny-1,0] ;; Bx
data2[*,0:ny-1,8] = -data2[*,0:ny-1,8] ;; Vy
data2[*,0:ny-1,9] = -data2[*,0:ny-1,9] ;; Vz

yy = findgen(ny*2)
yy[ny:2*ny-1]=y
yy[0:ny-1]=-reverse(y)

;; STARTING POINTS FOR FIELD LINES
npos=20

r0=fltarr(2,npos)
r1=fltarr(2,npos)
r2=fltarr(2,npos)
r0[0,*]=(findgen(npos)+1)*nx/(npos+1.)
r0[1,*]=ny-1+0.5
r1[0,*]=0.0
r1[1,*]=(findgen(npos)+1)*ny/(npos+1.)
r2[0,*]=0.0
r2[1,*]=(findgen(npos)+1)*ny/(npos+1.)+ny
fl1 = field_lines_2d(data2[*,*,0],data2[*,*,1],r0=r0,npos=npos,len=[0.5,0.5],dir=+1)
fl2 = field_lines_2d(data2[*,*,0],data2[*,*,1],r0=r0,npos=npos,len=[0.5,0.5],dir=-1)
fl3 = field_lines_2d(data2[*,*,0],data2[*,*,1],r0=r1,npos=npos,len=[0.5,0.5],dir=-1)
fl4 = field_lines_2d(data2[*,*,0],data2[*,*,1],r0=r2,npos=npos,len=[0.5,0.5],dir=+1)

;; DRAW 2D PROFILE WITH IMAGE, VECTOR, COLORBAR IDL NATIVE ROUTINES
img = image(data2[*,*,7],x-dx/2.,yy-dy/2,axis_style=2,xtickdir=1,ytickdir=1,rgb_table=33,yrange=[-15,15])
cb  = colorbar(target=img,orientation=0,tickdir=0,position=[0.2,0.1,0.8,0.15])
for k=0,npos-1 do begin
  pl=plot(min(x)+fl1[0,*,k]*max(x),min(yy)+fl1[1,*,k]*(max(yy)-min(yy)), /overplot, color='black',thick=1.5)
  pl=plot(min(x)+fl2[0,*,k]*max(x),min(yy)+fl2[1,*,k]*(max(yy)-min(yy)), /overplot, color='black',thick=1.5)
  pl=plot(min(x)+fl3[0,*,k]*max(x),min(yy)+fl3[1,*,k]*(max(yy)-min(yy)),/overplot, color='gray',thick=1.5)
  pl=plot(min(x)+fl4[0,*,k]*max(x),min(yy)+fl4[1,*,k]*(max(yy)-min(yy)),/overplot, color='gray',thick=1.5)
endfor

;;ADJUSTING FOTMATS
img.font_name='Times'
img.font_size='16'
img.xtitle='$X (\lambda)$' ;;CAN USE LATEX COMMANDS FOR TITLE
img.ytitle='$Y (\lambda)$'
img.title='$V_x$'
img.scale,1.0,2.0 ;; ARRANGE ASPECT RATIO
img.position=[0.15,0.25,0.95,0.85]
cb.font_name='Times'
cb.font_size='14'

img.save,'output.png',resolution=200 ;; FILE FORMAT IS AUTOMATICALLY DETECTED

end
