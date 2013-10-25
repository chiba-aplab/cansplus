pro put_time,time,charsize=charsize,xposition=xp,yposition=yp,color=color

if (n_elements(xp) eq 0) then xp=!x.window(0) 
if (n_elements(yp) eq 0) then begin
  dyw=!y.window(1)-!y.window(0)
  yp=!y.window(0)+1.01*dyw
endif
if (n_elements(color) eq 0) then begin
  if (!d.flags and 512) ne 0 then color=0 else color=!d.table_size-1
endif
if (n_elements(charsize) eq 0) then charsize=1

if time lt 100. then $
time_st=strmid(strcompress(string(time),/remove_all),0,4)
if time ge 100. then $
time_st=strmid(strcompress(string(time),/remove_all),0,5)

xyouts,/normal,xp,yp,'t = '+time_st,charsize=charsize,color=color

end
