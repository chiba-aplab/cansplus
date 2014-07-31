;; READING DATA
pushd,'data'
dac_read,data,x,y,z,'0011_*'
popd
dx = x[1]-x[0]
dy = y[1]-y[0]

;; CALCULATION OF TEMPERATURE
temperature  = data[*,*,5]/data[*,*,6]
info = size(temperature)
nx = info[1]
ny = info[2]
pos = ny*0.64*!pi/(2*!pi)

;; WINDOW SIZE
dimensions = [nx*3.0,ny*1.5]*2

;; DRAW
img = image(temperature,x-dx/2.,y-dy/2.,rgb_table=3,axis_style=2,layout=[2,1,1], $
            dimensions=dimensions,xtickdir=1,ytickdir=1)
plt = plot(x,temperature[*,pos],layout=[2,1,2],/current)

;; ADJUSTING FORMATS
img.xstyle=1
img.ystyle=1
img.xtitle='X'
img.ytitle='Y'
img.font_size=16
img.font_name='Times'

plt.xstyle=1
plt.xtitle='X'
plt.ytitle='Temperature ($\rho/P$)'
plt.title='$Y=0.64\pi$'
plt.font_name='Times'
plt.symbol='circle'
plt.font_size=16
plt.thick=2
plt.aspect_ratio = (2*!pi)/0.7

img.save,'OTvortex.png',resolution=200

end
