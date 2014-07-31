;; READING DATA
pushd,'data
dac_read,data,x,y,z,'0011*'
popd

;; PLOTS
pl1 = plot(x,data[*,5],layout=[2,2,1],margin=[0.2,0.1,0.05,0.1])
pl2 = plot(x,data[*,6],layout=[2,2,2],margin=[0.2,0.1,0.05,0.1],/current)
pl3 = plot(x,data[*,7],layout=[2,2,3],margin=[0.2,0.2,0.05,0.1],/current)
pl4 = plot(x,data[*,8],layout=[2,2,4],margin=[0.2,0.2,0.05,0.1],/current)

;; FOR ADJUSTING FORMATS
pl1.xtitle='$X$'
pl1.ytitle='$P$'
pl1.ystyle=2
pl1.thick=2
pl1.color='red'
pl1.font_name='Times'

pl2.xtitle='$X$'
pl2.ytitle='$\rho$'
pl2.ystyle=2
pl2.thick=2
pl2.color='blue'
pl2.font_name='Times'
pl2.symbol='circle'
pl2.sym_size=0.75

pl3.xtitle='$X$'
pl3.ytitle='$V_x$'
pl3.thick=2
pl3.color='purple'
pl3.font_name='Times'

pl4.xtitle='$X$'
pl4.ytitle='$V_y$'
pl4.thick=2
pl4.color='green'
pl4.font_name='Times'

pl1.save,'shktb.png',resolution=200 ;;DPI

end

