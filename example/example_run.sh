#!/usr/bin/bash

DATASET_DIR=./dataset

# Example01
EX_DIR01=${DATASET_DIR}/example01
cd $EX_DIR01
MGCplotter -r Mgallisepticum.gbff -o ./example_result01 --assign_cog_color
cd ../../

# Example02
EX_DIR02=${DATASET_DIR}/example02
cd $EX_DIR02
MGCplotter -r Mgallisepticum.gbff -o ./example_result02 --assign_cog_color \
           --query_files ./example02/*.gbff
cd ../../

# Gallery01
GALLERY_DIR01=${DATASET_DIR}/gallery01
cd $GALLERY_DIR01
MGCplotter -r ./ecoli.gbk -o ./gallery_result01 --rrna_color blue --trna_color red \
           --gc_content_p_color orange --gc_content_n_color blue \
           --gc_skew_p_color pink --gc_skew_n_color green
cd ../../

# Gallery02
GALLERY_DIR02=${DATASET_DIR}/gallery02
cd $GALLERY_DIR02
MGCplotter -r ./ecoli.gbk -o ./gallery_result02 --assign_cog_color \
           --query_files ./gallery02/NC_011751.gbk ./gallery02/NC_017634.gbk ./gallery02/NC_018658.gbk \
           --ticks_labelsize 50
cd ../../

# Gallery03
GALLERY_DIR03=${DATASET_DIR}/gallery03
cd $GALLERY_DIR03
MGCplotter -r ./Mgallisepticum.gbff -o ./gallery_result03 --assign_cog_color \
           --query_files ./gallery03/*.gbff --conserved_cds_color '#dc143c' \
           --rrna_r 0 --trna_r 0 --conserved_cds_r 0.01
cd ../../

# Gallery04
GALLERY_DIR04=${DATASET_DIR}/gallery04
cd $GALLERY_DIR04
MGCplotter -r ./Malvi.gbk -o ./gallery_result04 --assign_cog_color \
           --query_files ./gallery04/*.faa --conserved_cds_r 0.05 \
           --gc_content_r 0 --gc_skew_r 0
cd ../../

# Gallery05
GALLERY_DIR05=${DATASET_DIR}/gallery05
cd $GALLERY_DIR05
MGCplotter -r ./Mgallisepticum.gbff -o ./gallery_result05 --assign_cog_color \
           --cog_color_json ./cog_color.json
cd ../../
