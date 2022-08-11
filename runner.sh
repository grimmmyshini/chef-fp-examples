#!/bin/bash

rm /tmp/b.txt

bench_gen() {
    benchall -O3 ${1}.cpp -o ${1}.out -L$ws_dir/benchmark/build/src
}

alias gen_adc=bench_gen benchmark
alias gen_mix=bench_gen benchmark-mixed
alias gen_ovr=bench_gen benchmark-overhead

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
gen_ovr
run_adc ErrorEstimateArcLenClad
run_adc ErrorEstimateArcLenAdapt
run_mix ArcLenLowerPrec
run_mix ArcLenHighPrec
run_ovr ArcLen
rub_ovr ArcLenWithClad
cd ..

cd blackscholes
gen_adc
gen_ovr
run_adc ErrorEstimateBlkSolClad
run_adc ErrorEstimateBlkSolAdapt
run_ovr BlkSol
rub_ovr BlkSolWithClad
cd ..

cd HPCCG
gen_adc
gen_mix
gen_ovr
run_adc ErrorEstimateHPCCGClad
run_adc ErrorEstimateHPCCGAdapt
run_mix HPCCGLowerPrec
run_mix HPCCGHighPrec
run_ovr HPCCG
rub_ovr HPCCGWithClad
cd ..
echo "" >> results.txt

cd kmeans
gen_adc
gen_ovr
run_adc ErrorEstimateKMeansClad
run_adc ErrorEstimateKMeansAdapt
run_ovr KMeans
rub_ovr KMeansWithClad
cd ..

cd simpsons
gen_adc
gen_mix
gen_ovr
run_adc ErrorEstimateSimpClad
run_adc ErrorEstimateSimpAdapt
run_mix SimpLowerPrec
run_mix SimpHighPrec
run_ovr Simp
rub_ovr SimpWithClad
cd ..

cat /tmp/b.txt >> results.txt

echo "" >> results.txt
echo "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =" >> results.txt
echo "Finished at: $(date)" >> results.txt
echo "=============================================================================" >> results.txt
echo "" >> results.txt
echo "" >> results.txt
