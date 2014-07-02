;-----------------------------------------------------------------------
;  #  GAUSSIAN ELIMINATION
;-----------------------------------------------------------------------
pro gssel, dm,j
common pssn, p, x
common numarr, nx, ny

w = dcomplexarr(nx)
g = dcomplexarr(nx)

w[0] = 1.0/(-2.0*dm-1.0)
g[0] = p[0,j]*w[0]

for i=1,nx-2 do begin
   w[i] = 1.0/(-2.0*dm-w[i-1])
   g[i] = (p[i,j]-g[i-1])*w[i]
endfor
w[nx-1] = 1./(-2.0*dm-1.0-w[nx-2])
g[nx-1] = (p[nx-1,j]-g[nx-2])*w[nx-1]

x[nx-1,j] = g[nx-1]

for i=nx-2,0,-1 do begin
   x[i,j] = g[i]-w[i]*x[i+1,j]
endfor

end


;-----------------------------------------------------------------------
;     POISSON EQ. SOLVER FOR 2D BOUNDED-PERIODIC SIMULATION BOX
;-----------------------------------------------------------------------
pro poisson_bp, b, phi
common pssn, p,x
common numarr, nx,ny

;; SETING
info = size(b)
nx = info[1]
ny = info[2]
phi = dcomplexarr(nx,ny)

;-----------------------------------------------------------------------
;DIRECT-in-X, FFT-in-Y
;-----------------------------------------------------------------------
x = dcomplexarr(nx,ny)
aa = dcomplexarr(ny)
p  = dcomplexarr(nx,ny)

for i=0,nx-1 do begin
    for j=0,ny-1 do begin
       aa[j] = b[i,j]
    endfor

    bb = fft(aa,-1)

    for j=0,ny-1 do begin
       p[i,j] = bb[j]
    endfor
endfor

j=0
dm = 1.0
gssel,dm,j

j=ny/2
dm = 1.0+2.0*sin(!pi/2.0)^2
gssel,dm,j

for j=1,ny/2-1 do begin
    jj = ny-j
    argj = float(j)*!pi/float(ny)
    dm = 1.0+2.0*sin(argj)^2
    gssel,dm,j
    gssel,dm,jj
endfor

for i=0,nx-1 do begin
    for j=0,ny-1 do begin
        bb[j] = x[i,j]
    endfor

    aa = fft(bb,1)

    for j=0,ny-1 do begin
        phi[i,j] = aa[j]
    endfor
endfor


end

