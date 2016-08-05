;; Initilization when idl starts up.

;; set path environment
!path = expand_path('./') +':'+!path
!path = getenv('CANSPLUS_DIR')+'/idl/'+':'+!path

