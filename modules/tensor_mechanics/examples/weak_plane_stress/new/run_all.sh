#!/bin/sh
./run_elas_ref_study.sh
./run_plas_ref_study.sh
./run_interp_options.sh

find . -name '*.png'|xargs -L 1 mogrify -trim -border 5 -bordercolor white
