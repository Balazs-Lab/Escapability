#!/bin/bash

echo "==============================="
echo " Setting Up Environment"
echo "==============================="

conda init
conda activate viral_escape

echo "==============================="
echo " Escape Summary Analysis"
echo "==============================="

python3 analysis/1-average_mutations.py

echo "==============================="
echo " Haplotype Analysis"
echo "==============================="

python3 analysis/2-haplotype_analysis.py

echo "==============================="
echo " Sample Specific Analysis"
echo "==============================="

python3 analysis/3-sample_summary.py
