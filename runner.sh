#!/bin/bash

rm /tmp/b.txt /tmp/m.txt /tmp/o.txt

bench_gen() {
    benchall -O3 ${1}.cpp -o ${1}.out -L$ws_dir/benchmark/build/src
}

gen_adc() {
    bench_gen benchmark
}

gen_mix() {
    bench_gen benchmark-mixed
}

gen_ovr() {
    bench_gen benchmark-overhead
}

bench_run() {
    $1 --benchmark_filter=${2}\$ | tail -n +4 >>$3
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

echo "=============================================================================" >>results.txt
echo "Started at: $(date)" >>results.txt
echo "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =" >>results.txt
echo "" >>results.txt

cd arclength
gen_adc
gen_mix
run_adc ErrorEstimateArcLen/10000
run_adc ErrorEstimateArcLen/100000
run_adc ErrorEstimateArcLen/1000000
run_adc ErrorEstimateArcLen/10000000
run_adc ErrorEstimateArcLen/100000000
run_adc ErrorEstimateArcLenClad/10000
run_adc ErrorEstimateArcLenClad/100000
run_adc ErrorEstimateArcLenClad/1000000
run_adc ErrorEstimateArcLenClad/10000000
run_adc ErrorEstimateArcLenClad/100000000
run_adc ErrorEstimateArcLenAdapt/10000
run_adc ErrorEstimateArcLenAdapt/100000
run_adc ErrorEstimateArcLenAdapt/1000000
run_adc ErrorEstimateArcLenAdapt/10000000
run_mix ArcLenLowerPrec
run_mix ArcLenHighPrec
cd ..

cd blackscholes
gen_adc
ErrorEstimateBlkSol/0
ErrorEstimateBlkSol/1
ErrorEstimateBlkSol/2
ErrorEstimateBlkSol/3
ErrorEstimateBlkSol/4
ErrorEstimateBlkSolClad/0
ErrorEstimateBlkSolClad/1
ErrorEstimateBlkSolClad/2
ErrorEstimateBlkSolClad/3
ErrorEstimateBlkSolClad/4
ErrorEstimateBlkSolAdapt/0
ErrorEstimateBlkSolAdapt/1
ErrorEstimateBlkSolAdapt/2
cd ..

cd HPCCG
gen_adc
gen_mix
run_adc ErrorEstimateHPCCG/20/30/10
run_adc ErrorEstimateHPCCG/20/30/20
run_adc ErrorEstimateHPCCG/20/30/40
run_adc ErrorEstimateHPCCG/20/30/80
run_adc ErrorEstimateHPCCG/20/30/160
run_adc ErrorEstimateHPCCG/20/30/320
run_adc ErrorEstimateHPCCGClad/20/30/10
run_adc ErrorEstimateHPCCGClad/20/30/20
run_adc ErrorEstimateHPCCGClad/20/30/40
run_adc ErrorEstimateHPCCGClad/20/30/80
run_adc ErrorEstimateHPCCGClad/20/30/160
run_adc ErrorEstimateHPCCGClad/20/30/320
run_adc ErrorEstimateHPCCGAdapt/20/30/10
run_adc ErrorEstimateHPCCGAdapt/20/30/20
run_adc ErrorEstimateHPCCGAdapt/20/30/40
run_adc ErrorEstimateHPCCGAdapt/20/30/80
run_adc ErrorEstimateHPCCGAdapt/20/30/160
run_adc ErrorEstimateHPCCGAdapt/20/30/320
run_mix HPCCGLowerPrec
run_mix HPCCGHighPrec
cd ..

cd kmeans
gen_adc
run_adc ErrorEstimateKMeans/0
run_adc ErrorEstimateKMeans/1
run_adc ErrorEstimateKMeans/2
run_adc ErrorEstimateKMeans/3
run_adc ErrorEstimateKMeans/4
run_adc ErrorEstimateKMeansClad/0
run_adc ErrorEstimateKMeansClad/1
run_adc ErrorEstimateKMeansClad/2
run_adc ErrorEstimateKMeansClad/3
run_adc ErrorEstimateKMeansClad/4
run_adc ErrorEstimateKMeansAdapt/0
run_adc ErrorEstimateKMeansAdapt/1
run_adc ErrorEstimateKMeansAdapt/2
run_adc ErrorEstimateKMeansAdapt/3
cd ..

cd simpsons
gen_adc
gen_mix
run_adc ErrorEstimateSimp/10000
run_adc ErrorEstimateSimp/100000
run_adc ErrorEstimateSimp/1000000
run_adc ErrorEstimateSimp/10000000
run_adc ErrorEstimateSimp/100000000
run_adc ErrorEstimateSimpClad/10000
run_adc ErrorEstimateSimpClad/100000
run_adc ErrorEstimateSimpClad/1000000
run_adc ErrorEstimateSimpClad/10000000
run_adc ErrorEstimateSimpClad/100000000
run_adc ErrorEstimateSimpAdapt/10000
run_adc ErrorEstimateSimpAdapt/100000
run_adc ErrorEstimateSimpAdapt/1000000
run_adc ErrorEstimateSimpAdapt/10000000
run_mix SimpLowerPrec/10000
run_mix SimpLowerPrec/100000
run_mix SimpLowerPrec/1000000
run_mix SimpLowerPrec/10000000
run_mix SimpLowerPrec/100000000
run_mix SimpHighPrec/10000
run_mix SimpHighPrec/100000
run_mix SimpHighPrec/1000000
run_mix SimpHighPrec/10000000
run_mix SimpHighPrec/100000000
cd ..

echo "------------------------------- Adapt VS Clad -------------------------------" >>results.txt
cat /tmp/b.txt >>results.txt
echo "" >>results.txt
echo "------------------------------ Mixed Precision ------------------------------" >>results.txt
cat /tmp/m.txt >>results.txt
echo "" >>results.txt

echo "" >>results.txt
echo "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =" >>results.txt
echo "Finished at: $(date)" >>results.txt
echo "=============================================================================" >>results.txt
echo "" >>results.txt
echo "" >>results.txt
