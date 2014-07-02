pro set_window, xs, xe, ys, ye
common wdw, px,py,sx,sy

;; set the size of the window

!x.window(0) = xs
!y.window(0) = ys
!x.window(1) = xe
!y.window(1) = ye
px = !x.window * !d.x_vsize
py = !y.window * !d.y_vsize
sx = px(1)-px(0)
sy = py(1)-py(0)
sx = sx+0.5
sy = sy+0.5

end
