pro ldmyct,clr_index=clr,white=white,black=black,clr_bw=clr_bw,clr_wb=clr_wb

; Define private color map

n_clr=17
;clr=(!d.n_colors<256)/n_clr*lindgen(n_clr)
clr=lindgen(n_clr)


c_red  =bindgen(256)
c_green=bindgen(256)
c_blue =bindgen(256)

cbg=clr(0)  & ced=clr(1)-1  & c_red(cbg:ced)=0b   & c_green(cbg:ced)=0b   & c_blue(cbg:ced)=0b     ; black

cbg=clr(1)  & ced=clr(2)-1  & c_red(cbg:ced)=0b   & c_green(cbg:ced)=0b   & c_blue(cbg:ced)=84b    ; blue
cbg=clr(2)  & ced=clr(3)-1  & c_red(cbg:ced)=0b   & c_green(cbg:ced)=0b   & c_blue(cbg:ced)=168b
cbg=clr(3)  & ced=clr(4)-1  & c_red(cbg:ced)=0b   & c_green(cbg:ced)=0b   & c_blue(cbg:ced)=255b
cbg=clr(4)  & ced=clr(5)-1  & c_red(cbg:ced)=0b   & c_green(cbg:ced)=84b  & c_blue(cbg:ced)=255b
cbg=clr(5)  & ced=clr(6)-1  & c_red(cbg:ced)=0b   & c_green(cbg:ced)=168b & c_blue(cbg:ced)=255b
cbg=clr(6)  & ced=clr(7)-1  & c_red(cbg:ced)=0b   & c_green(cbg:ced)=200b & c_blue(cbg:ced)=255b

cbg=clr(7)  & ced=clr(8)-1  & c_red(cbg:ced)=0b   & c_green(cbg:ced)=225b & c_blue(cbg:ced)=170b   ; green
cbg=clr(8)  & ced=clr(9)-1  & c_red(cbg:ced)=0b   & c_green(cbg:ced)=255b & c_blue(cbg:ced)=0b
cbg=clr(9)  & ced=clr(10)-1 & c_red(cbg:ced)=120b & c_green(cbg:ced)=255b & c_blue(cbg:ced)=0b

cbg=clr(10) & ced=clr(11)-1 & c_red(cbg:ced)=220b & c_green(cbg:ced)=255b & c_blue(cbg:ced)=0b     ; red
cbg=clr(11) & ced=clr(12)-1 & c_red(cbg:ced)=255b & c_green(cbg:ced)=220b & c_blue(cbg:ced)=0b
cbg=clr(12) & ced=clr(13)-1 & c_red(cbg:ced)=255b & c_green(cbg:ced)=180b & c_blue(cbg:ced)=0b
cbg=clr(13) & ced=clr(14)-1 & c_red(cbg:ced)=255b & c_green(cbg:ced)=128b & c_blue(cbg:ced)=0b
cbg=clr(14) & ced=clr(15)-1 & c_red(cbg:ced)=255b & c_green(cbg:ced)=64b  & c_blue(cbg:ced)=0b
cbg=clr(15) & ced=clr(16)-1 & c_red(cbg:ced)=255b & c_green(cbg:ced)=0b   & c_blue(cbg:ced)=0b

;cbg=clr(16) & ced=255       & c_red(cbg:ced)=255b & c_green(cbg:ced)=255b & c_blue(cbg:ced)=255b   ; white

black=0b & white=!d.table_size-1
c_red(black)=0b   & c_green(black)=0b   & c_blue(black)=0b     ; black
c_red(white)=255b & c_green(white)=255b & c_blue(white)=255b   ; white

;red=byte(clr(15)) & green=byte(clr(8)) & blue=byte(clr(1))


tvlct,c_red,c_green,c_blue

;clr_bw=lindgen(17)*15
;clr_bw(16)=255

clr_bw=lindgen(17)*20<255
clr_wb=255-(lindgen(17)*14<255)
clr_wb(16)=0
;clr_wb=reverse([0,40+lindgen(15)*15,255])

return
end
