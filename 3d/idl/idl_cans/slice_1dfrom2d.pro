pro slice_1dfrom2d,data1d,s,data2d,xy0,xy1,mx=mx,x=x,y=y

sz=size(data2d)

case sz(0) of
  2: begin ; single time step
    nx=1
  end
  3: begin ; multiple time steps
    nx=sz(3)
  end
  else: begin ; illeagal
    print,' Error: The input data array is not 2-dimensional'
    return
  end
endcase

ix=sz(1) & jx=sz(2)

mx0=max(sz(1:2))
if (n_elements(mx) eq 0) then mx=mx0

if (n_elements(x) eq 0) then x=fltarr(ix)
if (n_elements(y) eq 0) then y=fltarr(jx)

x0=xy0(0) & y0=xy0(1)
x1=xy1(0) & y1=xy1(1)


xi=x0+(x1-x0)/float(mx-1)*findgen(mx)
yi=y0+(y1-y0)/float(mx-1)*findgen(mx)


s=sqrt((xi-x0)^2+(yi-y0)^2)

ires=interpol(findgen(ix),x,xi)
jres=interpol(findgen(jx),y,yi)

data1d=fltarr(mx,nx)

for n=0,nx-1 do begin
  data1d0=interpolate(data2d(*,*,n),ires,jres)
  data1d(*,n)=data1d0
endfor

data1d=reform(data1d)

return
end
