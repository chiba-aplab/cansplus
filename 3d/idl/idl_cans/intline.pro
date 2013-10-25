pro intline,Vx,Vy,Vz,X,Y,Z ,Xplot, Yplot, Zplot $
       ,index,xline,yline,zline $
       ,ds=ds,n_points=n_points
;+
; NAME:
;   intline
; PURPOSE:
;   make a field line vertices
; CATEGORY:
;   General graphics.
; CALLING SEQUENCE:
;   intline, vx, vy, vz, x, y, z, xplot,yplot,zplot
; INPUTS:
;  v[xyz] = 3 dimensional array indicating the [xyz]-component of a vector filed
;   [xyz] = 1 dimensional array indicating the [xyz]-coordinate
; OPTIONAL INPUT PARAMETERS:
;
; KEYWORD PARAMETERS:
;   ds = distance between the points in each line
;   n_points = number of points in each line
;
; OUTPUTS:
;   xline
; COMMON BLOCKS:
;       none.
; SIDE EFFECTS:
;       The currently selected display is affected.
; RESTRICTIONS:
;       None that are obvious.
; PROCEDURE:
;
; MODIFICATION HISTORY:
;    ver 0.0  T. yokoyama (Ntl Astron. Obs. Japan) 1997. 7.
;-


;;;;;;;;; Error checking ;;;;;;;;;;;;;;;
on_error,2                      ;Return to caller if an error occurs

;;;;;;;;; default parameters ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
szv = size(vx)             ;Size of image

if (n_elements(x) eq 0) then x=lindgen(szv(1))
if (n_elements(y) eq 0) then y=lindgen(szv(2))
if (n_elements(z) eq 0) then z=lindgen(szv(3))
if (n_elements(vlimit) eq 0) then vlimit=1.d-2

xrange=[min(x),max(x)]
yrange=[min(y),max(y)]
zrange=[min(z),max(z)]

if n_elements(n_points) eq 0 then n_points=2*max([szv(1),szv(2),szv(3)])
if n_elements(ds) eq 0 then $
    ds=0.5*min([float(xrange(1)-xrange(0))/szv(1) $
               ,float(yrange(1)-yrange(0))/szv(2) $
               ,float(zrange(1)-zrange(0))/szv(3)])

;;;;;;;;; Initialize ;;;;;;;;;;;;;;;

mlflg=0

;;;;;;;;; Line Start Position ;;;;;;;;;;;;;;;

      marx=n_elements(xplot)

;;;;;;;;; Lines ;;;;;;;;;;;;;;;

      index=lonarr(marx)
      mpoints=0L
      mline=0

        for mar=0,marx-1 do begin     ; lineID
           xp0=xplot(mar)
           yp0=yplot(mar)
           zp0=zplot(mar)

        for mbf=-1,1,2 do begin       ; direction

           xline0=fltarr(n_points)
           yline0=fltarr(n_points)
           zline0=fltarr(n_points)
           xline0(0)=xp0
           yline0(0)=yp0
           zline0(0)=zp0

        for m=1,n_points-1 do begin

           xp=xline0(m-1)
           yp=yline0(m-1)
           zp=zline0(m-1)

           if ((xp gt max(x)) or (xp lt min(x)) $    ; inside?
            or (yp gt max(y)) or (yp lt min(y)) $
            or (zp gt max(z)) or (zp lt min(z))) then goto,out

             whr=where((x le xp),count) & if (count ne 0) then ip0=max(whr)
             whr=where((y le yp),count) & if (count ne 0) then jp0=max(whr)
             whr=where((z le zp),count) & if (count ne 0) then kp0=max(whr)
             ip1=ip0+1 & jp1=jp0+1 & kp1=kp0+1
             ax0=(xp-x(ip0))/(x(ip1)-x(ip0)) & ax1=1.-ax0
             ay0=(yp-y(jp0))/(y(jp1)-y(jp0)) & ay1=1.-ay0
             az0=(zp-z(kp0))/(z(kp1)-z(kp0)) & az1=1.-az0

             vecx= ax1*ay1*az1*vx(ip0,jp0,kp0) $
                  +ax1*ay1*az0*vx(ip0,jp0,kp1) $
                  +ax1*ay0*az1*vx(ip0,jp1,kp0) $
                  +ax1*ay0*az0*vx(ip0,jp1,kp1) $
                  +ax0*ay1*az1*vx(ip1,jp0,kp0) $
                  +ax0*ay1*az0*vx(ip1,jp0,kp1) $
                  +ax0*ay0*az1*vx(ip1,jp1,kp0) $
                  +ax0*ay0*az0*vx(ip1,jp1,kp1)
             vecy= ax1*ay1*az1*vy(ip0,jp0,kp0) $
                  +ax1*ay1*az0*vy(ip0,jp0,kp1) $
                  +ax1*ay0*az1*vy(ip0,jp1,kp0) $
                  +ax1*ay0*az0*vy(ip0,jp1,kp1) $
                  +ax0*ay1*az1*vy(ip1,jp0,kp0) $
                  +ax0*ay1*az0*vy(ip1,jp0,kp1) $
                  +ax0*ay0*az1*vy(ip1,jp1,kp0) $
                  +ax0*ay0*az0*vy(ip1,jp1,kp1)
             vecz= ax1*ay1*az1*vz(ip0,jp0,kp0) $
                  +ax1*ay1*az0*vz(ip0,jp0,kp1) $
                  +ax1*ay0*az1*vz(ip0,jp1,kp0) $
                  +ax1*ay0*az0*vz(ip0,jp1,kp1) $
                  +ax0*ay1*az1*vz(ip1,jp0,kp0) $
                  +ax0*ay1*az0*vz(ip1,jp0,kp1) $
                  +ax0*ay0*az1*vz(ip1,jp1,kp0) $
                  +ax0*ay0*az0*vz(ip1,jp1,kp1)
             vec=sqrt(vecx^2+vecy^2+vecz^2)

           if (vec le vlimit) then begin
             goto,out
           endif

             dx=mbf*vecx/vec*ds
             dy=mbf*vecy/vec*ds
             dz=mbf*vecz/vec*ds

             xline0(m)=xp+dx
             yline0(m)=yp+dy
             zline0(m)=zp+dz


        endfor
out :
        mend=m
        if (mend le 1) then goto,skip
        mpoints=mpoints+mend
        xline0=xline0[0:mend-1]
        yline0=yline0[0:mend-1]
        zline0=zline0[0:mend-1]
        if (mbf eq -1) then begin
          xline0=reverse(xline0[1:*])
          yline0=reverse(yline0[1:*])
          zline0=reverse(zline0[1:*])
          mpoints=mpoints-1
        endif

;       if ((mbf eq -1) and (mar eq 0)) then begin
        if (mlflg eq 0) then begin
          xline=xline0
          yline=yline0
          zline=zline0
          mlflg=1
        endif else begin
          xline=[xline,xline0]
          yline=[yline,yline0]
          zline=[zline,zline0]
        endelse

skip :

;       stop
        endfor ; mbf
         index[mline]=mpoints-1
         mline=mline+1
        endfor ; mar

index=index[uniq(index)]

return
end
