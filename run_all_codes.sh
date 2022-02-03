# This is the master script for all codes.

Rscript ./code/step-1a__random-forest.R ./data/input-1__scrna_class.txt
wait
python3 ./code/step-1b__random-forest.py
wait

python3 ./code/step-2__gate_coverage_calc.py \
	./data/input-2a__scrna_annotation.txt \
	./data/input-2b__scrna_gc-matrix.txt
wait

python3 ./code/step-3a__cnn_gradcam.py \
	./data/input-2a__scrna_annotation.txt \
	./data/input-2b__scrna_gc-matrix.txt
wait
python3 ./code/step-3b__cnn_gradcam.py
wait
python3 ./code/step-3c__cnn_gradcam.py
wait

echo "All Jobs Done"

