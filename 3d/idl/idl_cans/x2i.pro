function x2i,xi,x0,Inverse=inverse

sz=size(xi) & ix=sz(1)

if (keyword_set(inverse)) then begin
  i0=ix-1
  for i=ix-1,0,-1 do begin
    if (xi(i) le x0) then begin
      i0=i
      goto,L111
    endif
  endfor
  L111 :
endif else begin
  i0=0
  for i=0,ix-1 do begin
    if (xi(i) ge x0) then begin
      i0=i
      goto,L101
    endif
  endfor
  L101 :
endelse

return,i0
end
