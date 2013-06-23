;======================================================================;
; Name: rddt_cart3d_xz_slice
;
; PURPOSE: Read data of xz slice
;
;======================================================================;

narg=-1

;----------------------------------------------------------------------;
;-- get parameters --;
dacgetparam, 'params.txt.0000', 'dtout', dt
dacgetparam, 'params.txt.0000', 'ix', ix
dacgetparam, 'params.txt.0000', 'jx', jx
dacgetparam, 'params.txt.0000', 'kx', kx
dacgetparam, 'params.txt.0000', 'margin', margin
dacgetparam, 'params.txt.0000', 'mpisize', npe
dacgetparam, 'params.txt.0000', 'te_factor', te_factor
dacgetparam, 'params.txt.0000', 'gm', gm

;-- get/set MPI coordinates in 3D Cartesian topology --;
dacgetparam, 'params.txt.0000', 'mpix', mpix
dacgetparam, 'params.txt.0000', 'mpiy', mpiy
dacgetparam, 'params.txt.0000', 'mpiz', mpiz
mpi_topo = lonarr(mpix,mpiy,mpiz)
files = 'params.txt.'+string(indgen(npe),form='(i4.4)')
for i=0, npe-1 do begin
   dacgetparam,files[i],'mpirank',mpirank
   dacgetparam,files[i],'mpitopox',mpitopox
   dacgetparam,files[i],'mpitopoy',mpitopoy
   dacgetparam,files[i],'mpitopoz',mpitopoz
   mpi_topo[mpitopox,mpitopoy,mpitopoz] = mpirank
endfor   

;-- read data (coord) and set global grid number --;
files = 'x.dac.'+string(mpi_topo[*,0,0],form='(i4.4)')
dacget1d,files,x,ix=igx,margin=margin

files = 'y.dac.'+string(mpi_topo[0,*,0],form='(i4.4)')
dacget1d,files,y,ix=jgx,margin=margin

files = 'z.dac.'+string(mpi_topo[0,0,*],form='(i4.4)')
dacget1d,files,z,ix=kgx,margin=margin

;-- read data (time step) --;
files = 't.dac.0000'
dacget0s, files, t, narg=narg



;----------------------------------------------------------------------;
;-- set parameters for reading --

;-- set default values of input parameters --;
nd_start_default = 1L
nd_end_default   = n_elements(t)
nd_inc_default   = 1L
y_slice_default  = long(jgx/2)

;-- set read time step --;
print,'Time step number:', n_elements(t)
nd_input=''
print, 'Input start/end/increment time step (separated by ,) [', $
       nd_start_default, ',', nd_end_default, ',',nd_inc_default,' ]', $
       format='(a55,i3,a1,i3,a1,i3,a2,$)'

read, nd_input
if nd_input eq '' then begin
   nd_start = nd_start_default
   nd_end   = nd_end_default
   nd_inc   = nd_inc_default
endif else begin
   nd_split = strsplit(nd_input, ',', /extract)
   nd_start = long(nd_split[0])
   nd_end   = long(nd_split[1])
   nd_inc   = long(nd_split[2])
endelse

;-- check input --;
while ( (nd_end gt n_elements(t)) or $
        (nd_start lt 1) or $
        (nd_start gt nd_end) $
      ) do begin
   print, 'Invalid time step. Reinput', format='(a26,$)'
   read, ans
   nd_split = strsplit(ans, ',', /extract)
   nd_start = long(nd_split[0])
   nd_end   = long(nd_split[1])
   nd_inc   = long(nd_split[2])
endwhile
print,'Start time step:', nd_start, format='(a20,i3)'
print,'End time step:', nd_end, format='(a20,i3)'
print,'Increment time step:', nd_inc, format='(a20,i3)'

;-- set y-grid for xz slice --;
y_slice_input = ''
print, 'y-grid:',jgx
print, 'Input y-grid for slice [', y_slice_default, ' ]', $
       format='(a25,i3,a2,$)'
read, y_slice_input
if y_slice_input eq '' then begin
   y_slice = y_slice_default
endif else begin
   while (long(y_slice_input) gt jgx) do begin
      print, 'Invalid y-grid. Reinput', format='(a24,$)'
      read, y_slice_input
   endwhile
   y_slice = long(y_slice_input)
endelse
print,'Slice y-grid:', y_slice, format='(a20,i3)'

;-- get variables name --;
var_name = file_search('*.dac.0001.0000')
for ivar = 0, n_elements(var_name)-1 do begin
   var_name_split = strsplit(var_name[ivar], '.', /extract)
   var_name[ivar] = var_name_split[0]
endfor

;-- set read variables --;
print,'Variables: ',var_name
print, 'Input variables (separated by ,) [all]', format='(a39,$)'
var_name_input=''
read, var_name_input

;-- check input --;
if (var_name_input ne 'all') and (var_name_input ne '') then begin

   var_name_input_split = strsplit(var_name_input, ',', /extract)
   flag_var_name = -1

   while flag_var_name eq -1 do begin
      flag_var_name = 0

      for ivar = 0, n_elements(var_name_input_split) -1 do begin
         if where(var_name eq var_name_input_split[ivar]) eq -1 then begin
            flag_var_name = -1
         endif
      endfor

      if flag_var_name eq -1 then begin
         print,'Invalid variable. Reinput', format='(a26,$)'
         read, var_name_input
         var_name_input_split = strsplit(var_name_input, ',', /extract)
      endif
   endwhile
   
   var_name = var_name_input_split
endif

print,'Read variablres:', var_name



;----------------------------------------------------------------------;
;-- read pysical variables --

;>>> ivar loop >>>;
for ivar = 0, n_elements(var_name) - 1 do begin

;>>> ind loop >>>;   
   for ind = nd_start, nd_end, nd_inc do begin

;-- set read file name --;      
      files = var_name[ivar] + '.dac.' + $
              string(ind,form='(i4.4)') + '.' + string(indgen(npe), $
                                                       form='(i4.4)')
      print,'Reading ',var_name[ivar] + '.dac.' + $
              string(ind,form='(i4.4)') + '.' + '****'

;-- read file --;
      dacget3s_cart3d, files, var_tmp_tmp, narg = narg $
                    , ix = ix, jx = jx, kx = kx, margin = margin $
                    , mpi_topo = mpi_topo
      if ind eq nd_start then begin
         var_tmp = reform(var_tmp_tmp[*, y_slice, *])
      endif else begin
         var_tmp = [ [ [var_tmp] ], [ [reform(var_tmp_tmp[*, y_slice, *])] ] ]
      endelse
   endfor
;<<< ind loop <<<;

;-- store data into structure for data --;
   if ivar eq 0 then begin
      var = create_struct(var_name[ivar],var_tmp)
   endif else begin
      var = create_struct(var, var_name[ivar], var_tmp)
   endelse

endfor
;<<< ivar loop <<<;

;-- add other physical variable --
var_tmp = te_factor * var.pr / var.ro
var = create_struct(var, 'te', var_tmp)

print,'Store read data ([x-grid, ygrid, time step]) into structure var'
help,/struct,var

filenum = 1

end
