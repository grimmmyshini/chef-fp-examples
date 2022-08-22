#!/bin/bash

rm /tmp/b.txt /tmp/m.txt /tmp/o.txt

bench_gen() {
    benchall -O3 ${1}.cpp ${2} -o ${1}.out -L$ws_dir/benchmark/build/src
}

gen_adc() {
    bench_gen benchmark -DMEMORY_PROFILER
}

gen_mix() {
    bench_gen benchmark-mixed
}

bench_run() {
    $1 --benchmark_filter=${2}\$ > $3
}

run_adc() {
    bench_run ./benchmark.out $1 /tmp/b.json
    python3 process_json.py /tmp/b.json >> $ws_dir/results.txt
}

run_mix() {
    bench_run ./benchmark-mixed.out $1 /tmp/m.json
    python3 process_json.py /tmp/m.json >> $ws_dir/results.txt
}

echo "" > $ws_dir/results.txt

cd arclength
gen_adc
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
cd ..

cd blackscholes
gen_adc
run_adc ErrorEstimateBlkSol/0
run_adc ErrorEstimateBlkSol/1
run_adc ErrorEstimateBlkSol/2
run_adc ErrorEstimateBlkSol/3
run_adc ErrorEstimateBlkSol/4
run_adc ErrorEstimateBlkSolClad/0
run_adc ErrorEstimateBlkSolClad/1
run_adc ErrorEstimateBlkSolClad/2
run_adc ErrorEstimateBlkSolClad/3
run_adc ErrorEstimateBlkSolClad/4
run_adc ErrorEstimateBlkSolAdapt/0
run_adc ErrorEstimateBlkSolAdapt/1
run_adc ErrorEstimateBlkSolAdapt/2
cd ..

cd HPCCG
gen_adc
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
cd ..

cd arclength
gen_mix
run_mix ArcLenLowerPrec2/10000    
run_mix ArcLenLowerPrec2/100000   
run_mix ArcLenLowerPrec2/1000000  
run_mix ArcLenLowerPrec2/10000000 
run_mix ArcLenLowerPrec2/100000000
run_mix ArcLenLowerPrec/10000     
run_mix ArcLenLowerPrec/100000    
run_mix ArcLenLowerPrec/1000000   
run_mix ArcLenLowerPrec/10000000  
run_mix ArcLenLowerPrec/100000000 
run_mix ArcLenHighPrec/10000      
run_mix ArcLenHighPrec/100000     
run_mix ArcLenHighPrec/1000000    
run_mix ArcLenHighPrec/10000000   
run_mix ArcLenHighPrec/100000000  
cd ..

cd simpsons
gen_mix
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
