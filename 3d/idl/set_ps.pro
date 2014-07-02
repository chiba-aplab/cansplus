pro set_ps,file,landscape=landscape

if( keyword_set(landscape) )then begin
    a = PSWINDOW(pagesize='a4',/landscape)
    set_plot,'ps'
    device,_Extra=a,/color,/times,filename=file,/enc,bits=8,/landscape
endif else begin
    a = PSWINDOW(pagesize='a4')
    set_plot,'ps'
    device,_Extra=a,/color,/times,filename=file,/enc,bits=8
endelse



return

end
