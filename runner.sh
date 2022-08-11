#!/bin/bash

rm /tmp/b.txt /tmp/m.txt /tmp/o.txt

bench_gen() {
    benchall -O3 ${1}.cpp -o ${1}.out -L$ws_dir/benchmark/build/src
}

gen_adc () {
    bench_gen benchmark
}

gen_mix () {
    bench_gen benchmark-mixed
}

gen_ovr () {
    bench_gen benchmark-overhead
}

bench_run() {
    $1 --benchmark_filter=$2 | tail -n +4 >> $3
}

run_adc() {
    bench_run ./benchmark.out $1 /tmp/b.txt
}

run_mix() {
    bench_run ./benchmark-mixed.out $1 /tmp/m.txt
}

run_ovr() {
    bench_run ./benchmark-overhead.out $1 /tmp/o.txt
}

echo "=============================================================================" >> results.txt
echo "Started at: $(date)" >> results.txt
echo "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =" >> results.txt
echo "" >> results.txt


cd arclength
gen_adc
gen_mix
run_adc ErrorEstimateArcLen
run_adc ErrorEstimateArcLenClad
run_adc ErrorEstimateArcLenAdapt
run_mix ArcLenLowerPrec
run_mix ArcLenHighPrec
cd ..

cd blackscholes
gen_adc
run_adc ErrorEstimateBlkSol
run_adc ErrorEstimateBlkSolClad
run_adc ErrorEstimateBlkSolAdapt
cd ..

cd HPCCG
gen_adc
gen_mix
run_adc ErrorEstimateHPCCG
run_adc ErrorEstimateHPCCGClad
run_adc ErrorEstimateHPCCGAdapt
run_mix HPCCGLowerPrec
run_mix HPCCGHighPrec
cd ..

cd kmeans
gen_adc
run_adc ErrorEstimateKMeans
run_adc ErrorEstimateKMeansClad
run_adc ErrorEstimateKMeansAdapt
cd ..

cd simpsons
gen_adc
gen_mix
run_adc ErrorEstimateSimp
run_adc ErrorEstimateSimpClad
run_adc ErrorEstimateSimpAdapt
run_mix SimpLowerPrec
run_mix SimpHighPrec
cd ..

echo "------------------------------- Adapt VS Clad -------------------------------" >> results.txt
cat /tmp/b.txt >> results.txt
echo "" >> results.txt
echo "------------------------------ Mixed Precision ------------------------------" >> results.txt
cat /tmp/m.txt >> results.txt
echo "" >> results.txt

echo "" >> results.txt
echo "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =" >> results.txt
echo "Finished at: $(date)" >> results.txt
echo "=============================================================================" >> results.txt
echo "" >> results.txt
echo "" >> results.txt
