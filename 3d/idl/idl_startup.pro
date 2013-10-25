!path = '~/idl/idl_cans:' + !path
;!path = '~/idl_ono/:' + !path
!path = '~/idl/:' + !path

DEVICE,DECOMPOSED=0,TRUE_COLOR=24,RETAIN=2
w=0
ii=0
;xmin=-1.0 & xmax=1.0
;ymin=0 & ymax=0.5
;xmin=0d0 & xmax=3.0d0
;ymin=-2.5d0 & ymax=2.5d0
ov='0'
vect=0
zauto=0
;wsx=1800 & wsy=600
;wsx=640 & wsy=480
;wsx=1120 & wsy=480
;wsx=1000 & wsy=1000
;wsx=1600 & wsy=800 & xmin=0.0 & xmax=1.5 & ymin=0.0 & ymax=0.8 ;横長ver
wsx=1000 & wsy=1000 & xmin=0.0 & xmax=6.0 & ymin=-3.0 & ymax=3.0
;wsx=1000 & wsy=1000 & xmin=-6.0 & xmax=6.0 & ymin=-6.0 & ymax=6.0
ans=''
png=0
eps=0
filename=''
FNAME=''
pxi=0.1 & pxe=0.9 & pyi=0.15 & pye=0.85
;pxi=0.07 & pxe=0.85 & pyi=0.2 & pye=0.95
dtout=0.628d0
yi(ii)=ii+20000
min_fill=1
ct=46
;ct=4
