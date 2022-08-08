#!/bin/bash

echo "=============================================================================" >> results.txt
echo "Started at: $(date)" >> results.txt
echo "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =" >> results.txt
echo "" >> results.txt

cd arclength
./benchmark.out --benchmark_filter=ErrorEstimateArcLenClad >> ../results.txt
./benchmark.out --benchmark_filter=ErrorEstimateArcLenAdapt >> ../results.txt
cd ..
echo "" >> results.txt

cd blackscholes
./benchmark.out --benchmark_filter=ErrorEstimateBlkSolClad >> ../results.txt
./benchmark.out --benchmark_filter=ErrorEstimateBlkSolAdapt >> ../results.txt
cd ..
echo "" >> results.txt

cd HPCCG
./benchmark.out --benchmark_filter=ErrorEstimateHPCCGClad >> ../results.txt
./benchmark.out --benchmark_filter=ErrorEstimateHPCCGAdapt >> ../results.txt
cd ..
echo "" >> results.txt

cd kmeans
./benchmark.out --benchmark_filter=ErrorEstimateKMeansClad >> ../results.txt
./benchmark.out --benchmark_filter=ErrorEstimateKMeansAdapt >> ../results.txt
cd ..

echo "" >> results.txt
echo "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =" >> results.txt
echo "Finished at: $(date)" >> results.txt
echo "=============================================================================" >> results.txt
echo "" >> results.txt
echo "" >> results.txt
