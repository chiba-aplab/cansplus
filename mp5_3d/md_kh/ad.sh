#PBS -N ad5
#PBS -q XT4B
#PBS -j oe
#PBS -l mppwidth=128
#PBS -meb
#PBS -M ogawa@astro.s.chiba-u.ac.jp
cd /work/ogawatk/mp5_3d_2/
time aprun -n 128 ./a.out >& out.log
