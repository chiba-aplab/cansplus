function fftf,data,time,direction,_extra=extra,freq=freq

;+
; FFT with frequency array
; Arrays are shifted in ascending order
; output=fftf(input,time,direction[,freq=freq][,_extra=_extra])
;-

if n_elements(direction) eq 0 then direction=-1

ans=fft(data,direction,_extra=extra)

nn=n_elements(time)

nn21=0

if nn eq n_elements(data) then begin
   dt=time[1]-time[0]
   nn21=nn/2+1
   freq=findgen(nn)
   freq[nn21]=nn21-nn+findgen(nn21-2)
   freq=freq/(nn*dt)
   freq=shift(freq,-nn21)
endif

ans=shift(ans,-nn21)

; factor max(time)-min(time) is multiplied for normalization to meet parseval's law
if nn eq n_elements(data) then begin
   ans=ans*(max(time)-min(time))
   if direction eq +1 then  ans=ans/nn
endif

return,ans
end



