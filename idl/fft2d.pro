function fft2d,data,x,t,direction,_extra=extra,wnum=wnum,freq=freq

;+
; FFT in 2D with frequency and wavenumber arrays
; Arrays are shifted in ascending order
; output=fft2d(input,x,t,direction[,wnum=wnum][,freq=freq][,_extra=_extra])
;-

if n_elements(direction) eq 0 then direction=-1

ans=fft(data,direction,_extra=extra)

nx=n_elements(x)
nt=n_elements(t)
sized=size(data,/dim)

nx21=0
nt21=0

if nx eq sized[0] then begin
   dx=x[1]-x[0]
   nx21=nx/2+1
   wnum=indgen(nx)
   wnum[nx21]=nx21-nx+findgen(nx21-2)
   wnum=wnum/(nx*dx)
   wnum=shift(wnum,-nx21)
endif

if nt eq sized[1] then begin
   dt=t[1]-t[0]
   nt21=nt/2+1
   freq=indgen(nt)
   freq[nt21]=nt21-nt+findgen(nt21-2)
   freq=freq/(nt*dt)
   freq=shift(freq,-nt21)
endif

ans=reverse(shift(ans,-nx21+1,-nt21))

; factor (max(x)-min(x))*(max(t)-min(t)) is multiplied for
; normalization to meet parseval's law
if nx eq sized[0] then begin
   ans=ans*(max(x)-min(x))
   if direction eq +1 then ans=ans/nx
endif
if nt eq sized[1] then begin
   ans=ans*(max(t)-min(t))
   if direction eq +1 then ans=ans/nt
endif

return,ans
end

