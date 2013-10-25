function dev2nml,dev,xcord=xcord,ycord=ycord

nml=float(dev)/!d.x_size
if keyword_set(ycord) then nml=float(dev)/!d.y_size
  
return,nml
end
