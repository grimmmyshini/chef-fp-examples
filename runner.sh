#!/bin/bash

OUTPUT_PROCESSING_SCRIPT=$ws_dir/clad-fp-error-est-examples/process_output.py

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
    # Takes 3 arguments in order: benchmark executable name; benchmark_filter input; benchmark_out file path
    /usr/bin/time -v $1 --benchmark_filter=${2}\$ --benchmark_out_format=json --benchmark_out=$3 2>&1 | grep "Maximum resident set size" | sed "s/[^0-9]//g"
}

run_bench() {
    # Takes 3 arguments in order: benchmark executable name; benchmark_filter input/benchmark name; json file path/benchmark_out file path
    BENCHMARK_EXEC=$1
    BENCHMARK_NAME=$2
    JSON_FILE=$3
    FILE=$ws_dir/results.txt

    echo "==================== Running: $BENCHMARK_NAME ADAPT Vs. CLAD ===================="
    MEM_USAGE=$(bench_run $BENCHMARK_EXEC $BENCHMARK_NAME $JSON_FILE)
    python3 $OUTPUT_PROCESSING_SCRIPT $JSON_FILE $MEM_USAGE >>$FILE
}

run_adc() {
    # Takes 1 argument: benchmark_filter input/benchmark name;
    run_bench ./benchmark.out $1 /tmp/b.json
}

run_mix() {
    # Takes 1 argument: benchmark_filter input/benchmark name;
    run_bench ./benchmark-mixed.out $1 /tmp/m.json
}

echo "" >$ws_dir/results.txt

# cd arclength
# gen_adc
# run_adc ErrorEstimateArcLen/10000
# run_adc ErrorEstimateArcLen/100000
# run_adc ErrorEstimateArcLen/1000000
# run_adc ErrorEstimateArcLen/10000000
# run_adc ErrorEstimateArcLen/100000000
# run_adc ErrorEstimateArcLenClad/10000
# run_adc ErrorEstimateArcLenClad/100000
# run_adc ErrorEstimateArcLenClad/1000000
# run_adc ErrorEstimateArcLenClad/10000000
# run_adc ErrorEstimateArcLenClad/100000000
# run_adc ErrorEstimateArcLenAdapt/10000
# run_adc ErrorEstimateArcLenAdapt/100000
# run_adc ErrorEstimateArcLenAdapt/1000000
# run_adc ErrorEstimateArcLenAdapt/10000000
# cd ..

# echo "" >$ws_dir/results.txt

cd blackscholes
gen_adc
run_adc BlackScholes/0
run_adc BlackScholes/1
run_adc BlackScholes/2
run_adc BlackScholes/3
run_adc BlackScholes/4
run_adc BlackScholes_Clad/0
run_adc BlackScholes_Clad/1
run_adc BlackScholes_Clad/2
run_adc BlackScholes_Clad/3
run_adc BlackScholes_Clad/4
run_adc BlackScholes_Adapt/0
run_adc BlackScholes_Adapt/1
run_adc BlackScholes_Adapt/2
cd ..

# echo "" >$ws_dir/results.txt

# cd HPCCG
# gen_adc
# run_adc ErrorEstimateHPCCG/20/30/10
# run_adc ErrorEstimateHPCCG/20/30/20
# run_adc ErrorEstimateHPCCG/20/30/40
# run_adc ErrorEstimateHPCCG/20/30/80
# run_adc ErrorEstimateHPCCG/20/30/160
# run_adc ErrorEstimateHPCCG/20/30/320
# run_adc ErrorEstimateHPCCGClad/20/30/10
# run_adc ErrorEstimateHPCCGClad/20/30/20
# run_adc ErrorEstimateHPCCGClad/20/30/40
# run_adc ErrorEstimateHPCCGClad/20/30/80
# run_adc ErrorEstimateHPCCGClad/20/30/160
# run_adc ErrorEstimateHPCCGClad/20/30/320
# run_adc ErrorEstimateHPCCGAdapt/20/30/10
# run_adc ErrorEstimateHPCCGAdapt/20/30/20
# run_adc ErrorEstimateHPCCGAdapt/20/30/40
# run_adc ErrorEstimateHPCCGAdapt/20/30/80
# run_adc ErrorEstimateHPCCGAdapt/20/30/160
# run_adc ErrorEstimateHPCCGAdapt/20/30/320
# cd ..

# echo "" >$ws_dir/results.txt

# cd kmeans
# gen_adc
# run_adc ErrorEstimateKMeans/0
# run_adc ErrorEstimateKMeans/1
# run_adc ErrorEstimateKMeans/2
# run_adc ErrorEstimateKMeans/3
# run_adc ErrorEstimateKMeans/4
# run_adc ErrorEstimateKMeansClad/0
# run_adc ErrorEstimateKMeansClad/1
# run_adc ErrorEstimateKMeansClad/2
# run_adc ErrorEstimateKMeansClad/3
# run_adc ErrorEstimateKMeansClad/4
# run_adc ErrorEstimateKMeansAdapt/0
# run_adc ErrorEstimateKMeansAdapt/1
# run_adc ErrorEstimateKMeansAdapt/2
# run_adc ErrorEstimateKMeansAdapt/3
# cd ..

# echo "" >$ws_dir/results.txt

# cd simpsons
# gen_adc
# run_adc ErrorEstimateSimp/10000
# run_adc ErrorEstimateSimp/100000
# run_adc ErrorEstimateSimp/1000000
# run_adc ErrorEstimateSimp/10000000
# run_adc ErrorEstimateSimp/100000000
# run_adc ErrorEstimateSimpClad/10000
# run_adc ErrorEstimateSimpClad/100000
# run_adc ErrorEstimateSimpClad/1000000
# run_adc ErrorEstimateSimpClad/10000000
# run_adc ErrorEstimateSimpClad/100000000
# run_adc ErrorEstimateSimpAdapt/10000
# run_adc ErrorEstimateSimpAdapt/100000
# run_adc ErrorEstimateSimpAdapt/1000000
# run_adc ErrorEstimateSimpAdapt/10000000
# cd ..

# echo "" >$ws_dir/results.txt

# cd arclength
# gen_mix
# run_mix ArcLenLowerPrec2/10000
# run_mix ArcLenLowerPrec2/100000
# run_mix ArcLenLowerPrec2/1000000
# run_mix ArcLenLowerPrec2/10000000
# run_mix ArcLenLowerPrec2/100000000
# run_mix ArcLenLowerPrec/10000
# run_mix ArcLenLowerPrec/100000
# run_mix ArcLenLowerPrec/1000000
# run_mix ArcLenLowerPrec/10000000
# run_mix ArcLenLowerPrec/100000000
# run_mix ArcLenHighPrec/10000
# run_mix ArcLenHighPrec/100000
# run_mix ArcLenHighPrec/1000000
# run_mix ArcLenHighPrec/10000000
# run_mix ArcLenHighPrec/100000000
# cd ..

# echo "" >$ws_dir/results.txt

# cd simpsons
# gen_mix
# run_mix SimpLowerPrec/10000
# run_mix SimpLowerPrec/100000
# run_mix SimpLowerPrec/1000000
# run_mix SimpLowerPrec/10000000
# run_mix SimpLowerPrec/100000000
# run_mix SimpHighPrec/10000
# run_mix SimpHighPrec/100000
# run_mix SimpHighPrec/1000000
# run_mix SimpHighPrec/10000000
# run_mix SimpHighPrec/100000000
# cd ..