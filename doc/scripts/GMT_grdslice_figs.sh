#!/usr/bin/env bash
# Make 5 illustrations used by grdslice documentation
# Make a grid of synthetic seamounts for testing grdslice
gmt begin
	gmt grdseamount -R-0:30/0:30/-0:30/0:32 -I1m -Gsynth.grd -Cg -E -fg <<- EOF
	0		0		30	80	40	131
	-0.3	-0.3	70	25	16	60
	0.3		0.3		50	25	25	60
	-0.35	0.35	60	25	20	45
	0.35	-0.35	-30	25	15	47
	EOF
	gmt makecpt -Cjet -T5/130/5
	gmt figure GMT_grdslice_view ps
	gmt grdview synth.grd -JM12c -JZ3c -Qs -C -Wc0.5p -p155/25 -B -Bz -BWSneZ
	gmt figure GMT_grdslice_contours ps
	gmt grdcontour synth.grd -JM12c -A20+f10p+c5% -C10 -S8 -B -BWSrt -Gl-0:30/0:05/0:0/0.05,-0:22/0/-0:22/0:30,0:22/0/0:22/-0:30 -L10/130 -T+d8p/2p+l.. --MAP_FRAME_TYPE=plain
	gmt figure GMT_grdslice_products ps
	gmt grdslice synth.grd -C10 -Isynth_index.txt -Esynth_slice.txt -Fsynth_foundation.txt -T10 -L10/130 -A5 -nl -Q8 -S8 > synth_center.txt
	gmt set MAP_FRAME_TYPE plain MAP_TITLE_OFFSET 0
	gmt subplot begin 1x3 -R-0:30/0:30/-0:30/0:32 -JM5c -Fs5c/0 -SCb+tc -SRl -BWrtS -M16p
    	gmt plot synth_foundation.txt -W1p -Glightgray -c -B+tFoundations
    	gmt convert synth_index.txt -o0,1,4 | gmt text -F+f12p
    	gmt plot synth_slice.txt -W0.25p -c -B+tSlices
    	gmt convert synth_center.txt -Ef | gmt plot -Sc5p -Gred
    	gmt convert synth_center.txt -Ef | gmt text -F+f12p+l+jTR -Gwhite -Dj5p
		gmt convert synth_center.txt | gmt plot -SE -W0.25 -i0,1,7,5+s2,6+s2 -B+tEllipses -c
		gmt convert synth_center.txt -Ef | gmt plot -Sc5p -Gred
	gmt subplot end
	gmt figure GMT_grdslice_connect ps
	gmt plot3d synth_center.txt -R-0:30/0:30/-0:30/0:32/0/130 -JM12c -JZ3c -W2p -B -Bz -BWSneZ -p165/25
	gmt grdview synth.grd -Qs -C -Wc0.5p -p -t50
gmt end show
