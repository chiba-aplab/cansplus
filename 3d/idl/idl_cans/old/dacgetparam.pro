pro dacgetparam,file,name,value

; check header
;
openr,unit,file,/get_lun
tmp=FSTAT(unit) & filesize=tmp.size

str=''
m=0

while (not eof(unit)) do begin

  readf,unit,str
  if (strlen(str) ge 4) then begin
  if (strmid(str,0,4) eq '#dac') then begin
    strs=strsplit(str,' ',/extract)
    if (name eq strs[2]) then begin
      mtype=fix(strs[1])
      case mtype of
        2 : begin
          pos=strsplit(str,' ')
          value=strmid(str,pos[3])
        end
        4 : value=long(strs[3])
        5 : value=float(strs[3])
      endcase
      goto,label0
    endif
  endif
  endif

endwhile

label0:


close,unit
free_lun,unit

return
end
