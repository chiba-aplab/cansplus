pro plot_tube,xline,yline,zline,t3d=t3d,n_vertices=n_vertices,radius=radius

if (n_elements(n_vertices) eq 0) then n_vertices=4
if (n_elements(radius) eq 0) then radius=1.

; make a circle
;
mcircx=n_vertices & r0=radius
s0=2.*!pi/mcircx*findgen(mcircx)
xcirc=r0*cos(s0)
ycirc=r0*sin(s0)
zcirc=fltarr(mcircx)

; cylinder polygon indices
;
p=intarr(5,mcircx)
for m=0,mcircx-1 do begin
  p[0,m]=4
  p[1,m]=m
  p[4,m]=mcircx+m
  if (m ne mcircx-1) then begin
    p[2,m]=m+1
    p[3,m]=mcircx+m+1
  endif else begin
    p[2,m]=0
    p[3,m]=mcircx
  endelse
endfor

; define verces & polygons
;
mlx=n_elements(xline)
pol=intarr(5,mcircx*(mlx-1))
xvertar=fltarr(mcircx*mlx)
yvertar=fltarr(mcircx*mlx)
zvertar=fltarr(mcircx*mlx)

for ml=0,mlx-1 do begin

  if (ml ne mlx-1) then begin
    vect=[xline[ml+1]-xline[ml] $
       ,yline[ml+1]-yline[ml] $
       ,zline[ml+1]-zline[ml]]
  endif
  cosp=vect[0]/sqrt(vect[0]^2+vect[1]^2)
  sinp=vect[1]/sqrt(vect[0]^2+vect[1]^2)
  cost=vect[2]/sqrt(total(vect^2))
  sint=sqrt(vect[0]^2+vect[1]^2)/sqrt(total(vect^2))
  matrot=[[ cosp*cost,-sinp,cosp*sint] $
       ,[ sinp*cost, cosp,sinp*sint] $
       ,[     -sint,    0,     cost]]
  xyzcirc=[[xcirc],[ycirc],[zcirc]]#matrot

  xvert=reform(xyzcirc[*,0])+xline[ml]
  yvert=reform(xyzcirc[*,1])+yline[ml]
  zvert=reform(xyzcirc[*,2])+zline[ml]

  xvertar[mcircx*ml:mcircx*(ml+1)-1]=xvert
  yvertar[mcircx*ml:mcircx*(ml+1)-1]=yvert
  zvertar[mcircx*ml:mcircx*(ml+1)-1]=zvert

endfor

for ml=0,mlx-2 do pol[*,mcircx*ml:mcircx*(ml+1)-1]=p+mcircx*ml
pol[0,*]=4


res=polyshade(xvertar,yvertar,zvertar,pol,t3d=t3d)
if (!d.name ne 'Z') then tv,res


return
end
