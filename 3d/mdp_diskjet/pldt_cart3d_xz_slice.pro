;======================================================================;
; Name: plot_cart3d_xz_slice
;
; PURPOSE: Plot data of xz slice
;
;======================================================================;

;-- set default plot parameters --
!P.font=1
!P.thick=3 & !P.charsize=4 & !P.charthick=4
!x.style=1 & !y.style=1 & !z.style=1
!x.thick=2 & !y.thick=2

charsize=3.5

w=0
ii=0
ov='0'
;vect=0
vect=1
zauto = 0

wsx = 1000
wsy = 1000

png = 0
filename = ""

xmin =  0.0d0
xmax = 6.0d0
ymin = -3.0d0
ymax =  3.0d0

pxi = 0.1
pxe = 0.9
pyi = 0.15
pye = 0.85

dtout=0.628d0
yi=ii+20000
min_fill = 1
ct=46

;変数宣言
nxx=10 & ny=10 &length=50
xrange=[xmin,xmax] & yrange=[ymin,ymax]
max_value=1.0d23 & min_value=-1.0d23
nx=filenum
zp=fltarr(401)
xp=fltarr(2)
fp=fltarr(2,401)
xp(0)=0d0
xp(1)=1.0d0
x0=1.0d0

data2nlevels=100               ;data2の等高線のオーバープロット時の等高線の本数

vectnlevels=400             ;?
vectcolor=255               ;ベクトルの色
xgap1=15                    ;vect=1の時の、ベクトルのx方向の表示間隔（/セル数）
zgap1=15                    ;vect=1の時の、ベクトルのz方向の表示間隔（/セル数）
xgap2=0.2d0                 ;vect=2の時の、ベクトルのx方向の表示間隔（/長さ）
zgap2=0.2d0                 ;vect=2の時の、ベクトルのz方向の表示間隔（/長さ）

;ct=34				;カラーバーの色番号(34:虹色,23:虹色＋等高線)
ctfile="~/IDL/cans_idl/colors2.tbl" ;カラーバーのファイル指定



;************************************** 表示設定  *****************************************

;	プロットする物理量の選択、表示設定

;ここで選択した物理量はdata1に格納される
;data1の物理量はカラーマップで表示される

;この上に等高線としてオーバープロットしたい物理量がある場合は、
;予めopt.proを実行して、overplotで物理量を指定する必要がある(選択しない場合は0を入力)
;これはdata2に格納される

;******************************************************************************************

ans=''
;data1に格納される物理量の指定
read,' Variable for color-maps ? (lro,lpr,lte,v*,b*,ve,be,beta,ls  etc) : ',ans


;タイムステップを指定
;i:静止画を表示するステップ、もしくはアニメーションを開始するステップ
;nx:アニメーションを終了するステップ
read,' time,termination time = ',i,nx
iistart=i
iistop=nx
nx=nx+1


;表示の指定
;左から順に、静止画、アニメーション、最大値、最小値
read,'plot(0:picture,1:animation,2:max,3:min) : ',w


;ウインドウサイズの指定
;wsx,wsyはidl_startup.proで定義
;para.optを実行することでも変更可能
window,0,xsize=wsx,ysize=wsy





;***************************** data1に物理量を格納  **************************************


loadct,0
t=i/10d0

case ans of
   
   'ro' : begin
      data1=var.ro
      zmin=0d0 & zmax=1d0
      data1name='!9r!3'
   end

   'lro' : begin
      data1=alog10(var.ro)
      zmin=-3.0d0 & zmax=1d0
      data1name='log(!9r!3)'
   end

   'pr' : begin
      data1=var.pr
      zmin=-3d0 & zmax=1d0
      data1name='!16P!3'
   end

   'lpr' : begin
      data1=alog10(var.pr)
      zmin=-3.0d0 & zmax=0.5d0
      data1name='log(!16P!3)'
   end

   'te' : begin
      data1=var.te
      zmin=-3.0d0 & zmax=1.0d0
      data1name='!16T!3'
   end

   'lte' : begin
      data1=alog10(var.te)
      zmin=-4.0d0 & zmax=-1.0d0
      data1name='log(!16T!3)'
   end

   'cs' : begin
      data1=sqrt(gm*var.pr/var.ro)
      zmin=0d0 & zmax=1.5d0
      data1name='!16cs!3'
   end

   'lcs' : begin
      data1=sqrt(gm*var.pr/var.ro)
      data1=alog10(data1)
      zmin=-2d0 & zmax=1.0d0
      data1name='!16logcs!3'
   end

   'vv' : begin
      data1=var.vx^2+var.vy^2+var.vz^2
      data1=sqrt(data1)
      zmin=0d0 & zmax=5.0d0
      data1name='!16v!3'
   end

   'vx' : begin
      data1=var.vx
      zmin=-2.0d0 & zmax=2.0d0
      data1name='!16vx!3'
   end

   'vy' : begin
      data1=var.vy
      zmin=-1.0d0 & zmax=1.0d0
      data1name='!16vy!3'
   end

   'vz' : begin
      data1=var.vz
      zmin=-3.0d0 & zmax=3.0d0
      data1name='!16vz!3'
   end

   'bx' : begin
      data1=var.bx
      zmin=-1d-2 & zmax=1d-2
      data1name='!16Bx!3'
   end

   'by' : begin
      data1=var.by
      zmin=-1.0d-2 & zmax=1.0d-2
      data1name='!16By!3'
   end

   'bz' : begin
      data1=var.bz
      zmin=-1d-2 & zmax=1d-2
      data1name='!16Bz!3'
   end

   'lbz' : begin
      data1=alog10(var.bz)
      zmin=-1d-10 & zmax=1d-5
      data1name='!16logBz!3'
   end

   'ev' :begin
      data1=0.5*ro*(var.vx^2+var.vy^2+var.vz^2)
      zmin=0.0d0 & zmax=10.0d0
      data1name='!16E!Dknt!N!3'
   end 

   'eb' : begin
      data1=(var.bx^2+var.by^2+var.bz^2)/2.0d0
      zmin=0d0 & zmax=0.05d0
      data1name='!16E!Dmag!N!3'
   end 

   'bt' : begin
      data1=2.0d0*var.pr/(var.bx^2+var.by^2+var.bz^2+1.0d-40)
      zmin=0d0 & zmax=100.0d0
      data1name='!9b!3'
   end 

   'lbt' : begin
      data1=2.0d0*pr/(var.bx^2+var.by^2+var.bz^2+1.0d-40)
      data1=alog10(data1)
      zmin=0.0d0 & zmax=4.0d0
      data1name='log!9b!3'
   end 

   'pt' : begin
      data1=var.pr+(var.bx^2+var.by^2+var.bz^2)/2.0d0
      zmin=0d0 & zmax=1.0d0
      data1name='!16P+!16B!U2!N/2!3'
   end

   'ee' : begin
      data1=var.pr/(gm-1.0d0)+0.5d0*(var.ro*(var.vx^2+var.vy^2+var.vz^2)+var.bx^2+var.by^2+var.bz^2)
      zmin=-0.1d0 & zmax=10.0d0
      data1name='!16E!Dtotal!N!3'
   end

endcase



;****************** zauto ****************************
if zauto eq 1 then begin
zmin = min(data1) & zmax =max(data1)
endif
;*****************************************************


zrange1=[zmin,zmax]
df=(zmax-zmin+2.0d0)/400d0

for j=0,400 do begin
zp(j)=zmin+j*df-1.0d0
fp(*,j)=zmin+j*df-1.0d0
endfor

;********************************* data2に物理量を格納  *************************************


case ov of
    '0' : begin
    end
    'n' :begin
       data2=ro/1.6d-24
    end
  'lro' : begin
      data2=alog10(ro/1.6d-24)
      zrange2=[-2,1]
      nlevel=10
  end

  'pr' : begin
      data2=alog10(pr)
      zrange2=[-0.6,0.6]
  end

  'te' : begin
      data2=te
      zrange2=[-0.7,1.1]
  end

  'vx' : begin
      data2=vx
      zrange2=[-20.0,20.0]
  end

  'vy' : begin
      data2=vy
      zrange2=[-20.0,20.0]
  end

  'vz' : begin
      data2=vz
      zrange2=[-20.0,20.0]
  end

  'bx' : begin
      data2=bx
      zrange2=[-2.0,2.0]
  end

  'by' : begin
      data2=by
      zrange2=[-2.0,2.0]
  end

  'bz' : begin
      data2=bz
      zrange2=[-2.0,2.0]
  end
  
  've' :begin
      data2=0.5*ro*(vx^2+vy^2+vz^2)
      zrange2=[0,25]
  end

  'be' : begin
      data2=bx^2+by^2+bz^2
      zrange2=[0,5]
  end

  'beta' : begin
      data2=pr/(bx^2+by^2+bz^2)
      zrange2=[1,1]
;      ct=0
  end

  'ay' : begin
;     data2=ay
;     zrange2=[-5,5]
;     ct=0
      data2=ay
      level=[0,0.002,0.004,0.006,0.008,0.01]


  end 








endcase



if w le 1  then begin
if w eq 0 then nx=i+1









;************************************ data1の描画  ***************************************


for ii=iistart,nx-1 do begin


if(min_fill eq 1) then begin
 hoge=size(data1)
 sizex=hoge(1)
 sizez=hoge(2)
 for xx=0,sizex-1 do begin
  for zz=0,sizez-1 do begin
   if (data1(xx,zz,ii) lt zmin+abs(zmax-zmin)/50.0d0) then begin
 	data1(xx,zz,ii)=zmin+abs(zmax-zmin)/50.0d0
   endif
  endfor
 endfor
endif



loadct,0
;t=dt*ii/3d13
yi=ii+10000

;contour,data1(*,*,0),x*x0,z*x0,nlevels=2,xrange=xrange,yrange=yrange,zrange=zrange1,charsize=charsize,charthick=charthick,background=255,title='t='+strmid(ii*6.28d0/50.0d0,5,7),xtitle='r',ytitle='z',/isotropic,/fill,position=[pxi,pyi,pxe,pye],color=0

contour,data1(*,*,0),x*x0,z*x0,nlevels=2,xrange=xrange,yrange=yrange,zrange=zrange1,charsize=charsize,charthick=charthick,background=255,title='t='+strmid(ii*dt,5,7),xtitle='r',ytitle='z',/isotropic,/fill,position=[pxi,pyi,pxe,pye],color=0

loadct,ct,file=ctfile

contour,data1(*,*,ii),x*x0,z*x0,nlevels=100,xrange=xrange,yrange=yrange,zrange=zrange1,background=0,title='t='+strmid(ii,6,6)+'',max_value=max_value,/isotropic,/fill,/overplot


;カラーバーの表示
if (vect ne 1) then begin
if (vect ne 2) then begin
loadct,0
contour,fp,xp,zp,nlevels=400,xrange=[0,1],yrange=zrange1,zrange=zrange1,charsize=charsize,charthick=charthick,background=255,title=data1name,color=0,/fill,/noerase,position=[0.93,0.1,0.94,0.9],xstyle=1,xminor=1,xticks=1,xtickname=[' ',' ']
loadct,ct,file=ctfile
contour,fp,xp,zp,nlevels=400,xrange=[0,1],yrange=zrange1,zrange=zrange1,charsize=charsize,charthick=charthick,background=255,title=data1name,color=0,/fill,/noerase,position=[0.90,0.1,0.91,0.9],xstyle=1,xminor=1,xticks=1,xtickname=[' ',' '],/overplot
endif
endif






;**************** 速度ベクトルの描画 **********************

;opt.proを実行することで設定（初期設定はvect=0）

;vect=0:表示しない
;vect=1:セルの等間隔毎に表示（描画するセル間隔を一定にする）
;vect=2:x座標を等間隔に区切って表示（セルの幅は一定ではない）

;**********************************************************





;***********************************************

;   セルの等間隔毎にベクトルを表示するスタイル

;***********************************************


if vect eq 1 then begin

;xgap,zgap:描画対象にするセルの間隔
;ix,iz:元々の座標配列の要素数
xgap=xgap1
zgap=zgap1
ix=ix
iz=kx*npe-2*margin*(npe-1)

;xcount,zcount:間隔gap毎に間引いた後の配列の要素数
xcount=ix-(ix mod xgap)
xcount=xcount/xgap
zcount=iz-(iz mod zgap)
zcount=zcount/zgap

;vvx,vvz:間引いた配列の要素数を使って定義し直した速度ベクトルの成分
;xx,zz:同じく間引いた後の座標の配列
vvx=fltarr(xcount,zcount)
vvz=fltarr(xcount,zcount)
xx=fltarr(xcount)
zz=fltarr(zcount)

;xx,zzに対応する座標の値を代入
for l=0,xcount-1 do begin
xx(l)=x(l*xgap)
endfor

for m=0,zcount-1 do begin
zz(m)=z(m*zgap)
endfor

;vvx,vvzに対応する値を代入
for l=0,xcount-1 do begin
	for m=0,zcount-1 do begin
	vvx(l,m)=vx(l*xgap,m*zgap,ii)
	vvz(l,m)=vz(l*xgap,m*zgap,ii)
        endfor
endfor

;グラフの表示範囲に対応した配列番号を調べる
minxc=0
minzc=0
maxxc=0
maxzc=0

for l=0,xcount-1 do begin
	if(xx(l) le xmin) then begin
	minxc=l
	endif

	if(xx(l) le xmax) then begin
	maxxc=l
	endif
endfor

for m=0,zcount-1 do begin
	if(zz(m) le ymin) then begin
	minzc=m
	endif

	if(zz(m) le ymax) then begin	;縦軸はidl側ではzでなくyである
	maxzc=m
	endif
endfor
minxc=minxc+1
minzc=minzc+1

;velovectでベクトルをオーバープロト
loadct,0	;カラーセットを指定
velovect,vvx(minxc:maxxc,minzc:maxzc),vvz(minxc:maxxc,minzc:maxzc),xx(minxc:maxxc),zz(minzc:maxzc),/overplot,color=vectcolor

;カラーバーの表示
loadct,0
contour,fp,xp,zp,nlevels=vectnlevels,xrange=[0,1],yrange=zrange1,zrange=zrange1,charsize=charsize,charthick=charthick,background=255,title=data1name,color=0,/fill,/noerase,position=[0.90,0.1,0.91,0.9],xstyle=1,xminor=1,xticks=1,xtickname=[' ',' ']
loadct,ct,file=ctfile
contour,fp,xp,zp,nlevels=vectnlevels,xrange=[0,1],yrange=zrange1,zrange=zrange1,charsize=charsize,charthick=charthick,background=255,title=data1name,color=0,/fill,/noerase,position=[0.90,0.1,0.91,0.9],xstyle=1,xminor=1,xticks=1,xtickname=[' ',' '],/overplot



endif

;***********************************************

;   x座標の大体等間隔毎にベクトルを表示するスタイル

;***********************************************


if vect eq 2 then begin

;xgap,zgap:描画対象にするx座標の間隔
;ix,iz:元々の座標配列の要素数
;xxmin,xxmax:x座標の左端と右端の値
;zzmin,zzmax:z座標の下端と上端の値
xgap=xgap2
zgap=zgap2
ix=ix
iz=kx*npe-2*margin*(npe-1)
xxmin=x(0)
xxmax=x(ix-1)
zzmin=z(0)
zzmax=z(iz-1)

;xcount,zcount:間隔gap毎に間引いた後の配列の要素数

xarr=fltarr(ix)
zarr=fltarr(iz)
xcount=0
zcount=0

flag=0
bufx=xxmin+xgap
bufz=zzmin+zgap

xarr(xcount)=0
xcount=1
for l=1,ix-1 do begin
	if bufx le x(l) then begin
			xarr(xcount)=l
			bufx=x(l)+xgap
			xcount=xcount+1
	endif
endfor

zarr(zcount)=0
zcount=1
for m=1,iz-1 do begin
	if bufz le z(m) then begin
			zarr(zcount)=m
			bufz=z(m)+zgap
			zcount=zcount+1
	endif
endfor

;vvx,vvz:間引いた配列の要素数を使って定義し直した速度ベクトルの成分
;xx,zz:同じく間引いた後の座標の配列
vvx=fltarr(xcount,zcount)
vvz=fltarr(xcount,zcount)
xx=fltarr(xcount)
zz=fltarr(zcount)

;xx,zzに対応する座標の値を代入
for l=0,xcount-1 do begin
xx(l)=x(xarr(l))
endfor

for m=0,zcount-1 do begin
zz(m)=z(zarr(m))
endfor

;vvx,vvzに対応する値を代入
for l=0,xcount-1 do begin
	for m=0,zcount-1 do begin
	vvx(l,m)=vx(xarr(l),zarr(m),ii)
	vvz(l,m)=vz(xarr(l),zarr(m),ii)
        endfor
endfor

;グラフの表示範囲に対応した配列番号を調べる
minxc=0
minzc=0
maxxc=0
maxzc=0

for l=0,xcount-1 do begin
	if(xx(l) le xmin) then begin
	minxc=l
	endif

	if(xx(l) le xmax) then begin
	maxxc=l
	endif
endfor

for m=0,zcount-1 do begin
	if(zz(m) le ymin) then begin
	minzc=m
	endif

	if(zz(m) le ymax) then begin	;縦軸はidl側ではzでなくyである
	maxzc=m
	endif
endfor

;矢印が若干はみ出るので、表示する端のベクトルを若干個減らす
;minxc=minxc+1
;minzc=minzc+2
;maxxc=maxxc-2
;maxzc=maxzc-2

;velovectでベクトルをオーバープロト
loadct,0	;カラーセットを指定
velovect,vvx(minxc:maxxc,minzc:maxzc),vvz(minxc:maxxc,minzc:maxzc),xx(minxc:maxxc),zz(minzc:maxzc),/overplot,color=vectcolor

;カラーバーの表示
loadct,0
contour,fp,xp,zp,nlevels=vectnlevels,xrange=[0,1],yrange=zrange1,zrange=zrange1,charsize=charsize,charthick=charthick,background=255,title=data1name,color=0,/fill,/noerase,position=[0.90,0.1,0.91,0.9],xstyle=1,xminor=1,xticks=1,xtickname=[' ',' ']
loadct,ct,file=ctfile
contour,fp,xp,zp,nlevels=vectnlevels,xrange=[0,1],yrange=zrange1,zrange=zrange1,charsize=charsize,charthick=charthick,background=255,title=data1name,color=0,/fill,/noerase,position=[0.90,0.1,0.91,0.9],xstyle=1,xminor=1,xticks=1,xtickname=[' ',' '],/overplot

endif





;************************************ data2の描画  ***************************************


if ov ne '0' then begin
if vect eq 0 then begin
loadct,0
contour,data2(*,*,ii),x*x0,z*x0,nlevels=data2nlevels,xrange=xrange,yrange=yrange,zrange=zrange2,charsize=charsize,charthick=charthick,thick=3,background=0,title='',xtitle='',ytitle='',xtickname=replicate(' ',60),ytickname=replicate(' ',60),color=0,/isotropic,/noerase,position=[pxi,pyi,pxe,pye]
endif else begin
loadct,0
contour,data2(*,*,ii),x*x0,z*x0,nlevels=data2nlevels,xrange=xrange,yrange=yrange,zrange=zrange2,charsize=charsize,charthick=charthick,thick=2,background=0,title='',xtitle='',ytitle='',xtickname=replicate(' ',60),ytickname=replicate(' ',60),color=150,/isotropic,/noerase,position=[pxi,pyi,pxe,pye]
endelse
endif





;******************************************

;         画像をpng形式で保存

;******************************************
if png eq 1 then begin
;write_png,strmid(filename,0)+strmid(yi,5)+'.png',tvrd(/true)
str=string(yi,format='(i5)')
str=strmid(str,1,4)
write_png,filename+str+'.png',tvrd(/true)
endif

;wait,1
endfor
endif







;*************************** 最大値、最小値の出力等(未整理)  *********************************:

if w eq 2 then begin
maxdata=fltarr(nx)
for ii=0,nx-1 do begin
maxdata(ii)=max(alog10(data1(*,*,ii)))
endfor
plot,maxdata,background=255,color=0,charsize=2.0,charthick=2.0,xtitle="time",ytitle="log n",yrange=[0,4],thick=2

for ii=0,nx-1 do begin
maxdata(ii)=max(alog10(data1(30:ix-1,*,ii)))
endfor
plot,maxdata,background=255,color=0,charsize=2.0,charthick=2.0,xtitle="time",ytitle="log n",/noerase,yrange=[0,4],thick=2,linestyle=2


write_png,'maxn.png',tvrd(/true)
endif


if w eq 3 then begin
maxdata=fltarr(nx)
for ii=0,nx-1 do begin
maxdata(ii)=min(data1(*,*,ii))
endfor
plot,maxdata,background=255,color=0,charsize=2.0,charthick=2.0,xtitle="time",ytitle="min T",thick=3.0
write_png,'minte.png',tvrd(/true)
endif

if w eq 4 then begin
ii=i
cnx=101
cdata1=fltarr(cnx,cnx)
cx=fltarr(cnx)
cy=fltarr(cnx)
np=fltarr(cnx,cnx)
inix=-2d0 & enx=4d0 & iniy=2d0 & eny=5d0

stepx=(enx-inix)/(cnx-1) & stepy=(eny-iniy)/(cnx-1)
for ii=i,nx-1 do begin
loadct,0
contour,cdata1,alog10(cx),alog10(cy),yrange=[iniy,eny],xrange=[inix,enx],nlevel=100,background=255,color=0,thick=2,charsize=2.5,charthick=2.5,xtitle='log n',ytitle='log P/k!DB',/fill,zrange=[1d0,1d6],/isotropic,min_value=10

gmh=2.0d-26
lam1=7.3d-21
lam2=7.9d-27


for i=0,cnx-1 do begin
cx(i)=10.0d0^(inix+stepx*i)
endfor
for j=0,cnx-1 do begin
cy(j)=10.0d0^(iniy+stepy*j)
endfor
for i=1,cnx-2 do begin
for j=1,cnx-2 do begin
np(i,j)=0
tem=cy(j)/cx(i)
if tem le 1d4 then begin
cdata1(i,j)=cx(i)*(gmh-cx(i)*(lam1*exp(-118400d0/(tem+1500d0))+lam2*exp(-200d0/tem)))
endif else begin
cdata1(i,j)=exp(-(tem/1d4-1d0)^2)*cx(i)*(gmh-cx(i)*(lam1*exp(-118400d0/(tem+1500d0))+lam2*exp(-200d0/tem)))
endelse
;for jj=0,ix-1 do begin
;for kk=0,jkx-1 do begin
;if
oo=0
b=where(ro(oo:ix-1,*,ii) ge 0.5d0*(cx(i)+cx(i-1))*1.6d-24 and ro(oo:ix-1,*,ii) le 0.5d0*(cx(i+1)+cx(i))*1.6d-24 and pr(oo:ix-1,*,ii) ge 0.5d0*1.38d-16*(cy(j)+cy(j-1)) and pr(oo:ix-1,*,ii) le 0.5d0*1.38d-16*(cy(j+1)+cy(j)),bn)
;b=where(ro(*,*,ii) ge 0.5d0*(cx(i)+cx(i-1))*1.6d-24,bn)
np(i,j)=bn
;endif
;endfor
;endfor
print,ii,i,j
endfor
endfor

loadct,ct,file=ctfile

contour,np,alog10(cx),alog10(cy),yrange=[iniy,eny],xrange=[inix,enx],nlevel=500,thick=2,charsize=2.5,charthick=2.5,xtitle='log n',ytitle='log P/k!DB',/fill,/noerase,zrange=[1d0,1d2],/isotropic

loadct,0
contour,cdata1,alog10(cx),alog10(cy),yrange=[iniy,eny],xrange=[inix,enx],nlevel=1,background=255,color=0,thick=2,charsize=2.5,charthick=2.5,xtitle='log n',ytitle='log P/k!DB',level=[0],/noerase,/isotropic
;

if png eq 1 then begin
yi=ii+10000
write_png,strmid(filename,0)+strmid(yi,5)+'.png',tvrd(/true)
endif
endfor

endif
png=0
end
