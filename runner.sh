#!/bin/bash
if [ -z "${ws_dir}" ]; then
    export ws_dir=/code
    echo "Setting ws_dir"
fi
OUTPUT_PROCESSING_SCRIPT=$ws_dir/chef-fp-examples/process_output.py

benchall() {
    clang++-13 \
        -Xclang -add-plugin -Xclang clad -Xclang -load -Xclang clad.so \
        -Xclang -plugin-arg-clad -Xclang -fcustom-estimation-model \
        -Xclang -plugin-arg-clad -Xclang $ws_dir/chef-fp-examples/PrintModel/libPrintModel.so \
        -lstdc++ -lm -std=c++11 -O3 \
        "$@" \
        -L$ws_dir/benchmark/build/src \
        -lpthread -lbenchmark -DCODI_ZeroAdjointReverse=0
}

bench_gen() {
    benchall -O3 ${1}.cpp ${@:2} -o ${1}.out -L$ws_dir/benchmark/build/src
}

gen_adc() {
    bench_gen benchmark -DMEMORY_PROFILER $@
}

gen_mix() {
    bench_gen benchmark-mixed -DMEMORY_PROFILER $@
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
    rm -f JSON_FILE

    MEM_USAGE=$(bench_run $BENCHMARK_EXEC $BENCHMARK_NAME $JSON_FILE)
    python3 $OUTPUT_PROCESSING_SCRIPT $JSON_FILE $MEM_USAGE >>$FILE
}

run_adc() {
    # Takes 1 argument: benchmark_filter input/benchmark name;
    echo "==================== Running: $1 ADAPT Vs. CHEF-FP ===================="
    run_bench ./benchmark.out $1 /tmp/b.json
}

run_mix() {
    # Takes 1 argument: benchmark_filter input/benchmark name;
    echo "==================== Running: $1 Mixed Precision ===================="
    run_bench ./benchmark-mixed.out $1 /tmp/m.json
}

python3 $OUTPUT_PROCESSING_SCRIPT >$ws_dir/results.txt

# --------------------------------- ADAPT vs. CHEF-FP -----------------------------------

cd arclength
gen_adc
run_adc ArcLength/10000
run_adc ArcLength/100000
run_adc ArcLength/1000000
run_adc ArcLength/10000000
run_adc ArcLength/100000000
run_adc ArcLength_Clad/10000
run_adc ArcLength_Clad/100000
run_adc ArcLength_Clad/1000000
run_adc ArcLength_Clad/10000000
run_adc ArcLength_Clad/100000000
run_adc ArcLength_Adapt/10000
run_adc ArcLength_Adapt/100000
run_adc ArcLength_Adapt/1000000
run_adc ArcLength_Adapt/10000000
run_adc ArcLength_Adapt/100000000
cd ..

echo "" >>$ws_dir/results.txt

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
run_adc BlackScholes_Adapt/3
run_adc BlackScholes_Adapt/4
cd ..

echo "" >>$ws_dir/results.txt

cd HPCCG
gen_adc
run_adc HPCCG/20/30/10
run_adc HPCCG/20/30/20
run_adc HPCCG/20/30/40
run_adc HPCCG/20/30/80
run_adc HPCCG/20/30/160
run_adc HPCCG/20/30/320
run_adc HPCCG_Clad/20/30/10
run_adc HPCCG_Clad/20/30/20
run_adc HPCCG_Clad/20/30/40
run_adc HPCCG_Clad/20/30/80
run_adc HPCCG_Clad/20/30/160
run_adc HPCCG_Clad/20/30/320
run_adc HPCCG_Adapt/20/30/10
run_adc HPCCG_Adapt/20/30/20
run_adc HPCCG_Adapt/20/30/40
run_adc HPCCG_Adapt/20/30/80
run_adc HPCCG_Adapt/20/30/160
run_adc HPCCG_Adapt/20/30/320
cd ..

echo "" >>$ws_dir/results.txt

cd kmeans
gen_adc
run_adc KMeans/0
run_adc KMeans/1
run_adc KMeans/2
run_adc KMeans/3
run_adc KMeans/4
run_adc KMeans_Clad/0
run_adc KMeans_Clad/1
run_adc KMeans_Clad/2
run_adc KMeans_Clad/3
run_adc KMeans_Clad/4
run_adc KMeans_Adapt/0
run_adc KMeans_Adapt/1
run_adc KMeans_Adapt/2
run_adc KMeans_Adapt/3
run_adc KMeans_Adapt/4
cd ..

echo "" >>$ws_dir/results.txt

cd simpsons
gen_adc
run_adc Simpsons/10000
run_adc Simpsons/100000
run_adc Simpsons/1000000
run_adc Simpsons/10000000
run_adc Simpsons/100000000
run_adc Simpsons_Clad/10000
run_adc Simpsons_Clad/100000
run_adc Simpsons_Clad/1000000
run_adc Simpsons_Clad/10000000
run_adc Simpsons_Clad/100000000
run_adc Simpsons_Adapt/10000
run_adc Simpsons_Adapt/100000
run_adc Simpsons_Adapt/1000000
run_adc Simpsons_Adapt/10000000
run_adc Simpsons_Adapt/100000000
cd ..

# # -------------------------------- Mixed Precision -------------------------------

echo "" >>$ws_dir/results.txt

cd arclength
gen_mix
run_mix ArcLength_MixedPrecision/10000
run_mix ArcLength_MixedPrecision/100000
run_mix ArcLength_MixedPrecision/1000000
run_mix ArcLength_MixedPrecision/10000000
run_mix ArcLength_MixedPrecision/100000000
run_mix ArcLength_HighPrecision/10000
run_mix ArcLength_HighPrecision/100000
run_mix ArcLength_HighPrecision/1000000
run_mix ArcLength_HighPrecision/10000000
run_mix ArcLength_HighPrecision/100000000
cd ..

echo "" >>$ws_dir/results.txt

cd blackscholes
gen_mix
run_mix Blackscholes_MP/0/0
run_mix Blackscholes_MP/0/1
run_mix Blackscholes_MP/0/2
run_mix Blackscholes_MP/0/3
run_mix Blackscholes_MP/0/4
run_mix Blackscholes_MP/1/0
run_mix Blackscholes_MP/1/1
run_mix Blackscholes_MP/1/2
run_mix Blackscholes_MP/1/3
run_mix Blackscholes_MP/1/4
gen_mix -DFASTEXP=On
echo "With fast exp:" >>$ws_dir/results.txt
run_mix Blackscholes_MP/1/0
run_mix Blackscholes_MP/1/1
run_mix Blackscholes_MP/1/2
run_mix Blackscholes_MP/1/3
run_mix Blackscholes_MP/1/4
cd ..

echo "" >>$ws_dir/results.txt

cd HPCCG
gen_mix
run_mix HPCCG_MixedPrecision/20/30/10
run_mix HPCCG_MixedPrecision/20/30/20
run_mix HPCCG_MixedPrecision/20/30/40
run_mix HPCCG_MixedPrecision/20/30/80
run_mix HPCCG_MixedPrecision/20/30/160
run_mix HPCCG_MixedPrecision/20/30/320
run_mix HPCCG_HighPrecision/20/30/10
run_mix HPCCG_HighPrecision/20/30/20
run_mix HPCCG_HighPrecision/20/30/40
run_mix HPCCG_HighPrecision/20/30/80
run_mix HPCCG_HighPrecision/20/30/160
run_mix HPCCG_HighPrecision/20/30/320
cd ..

echo "" >>$ws_dir/results.txt

cd kmeans
gen_mix
run_mix EuclidDist_MixedPrecision/10/3
run_mix EuclidDist_MixedPrecision/10/10
run_mix EuclidDist_MixedPrecision/100/10
run_mix EuclidDist_MixedPrecision/100/100
run_mix EuclidDist_MixedPrecision/1000/100
run_mix EuclidDist_HighPrecision/10/3
run_mix EuclidDist_HighPrecision/10/10
run_mix EuclidDist_HighPrecision/100/10
run_mix EuclidDist_HighPrecision/100/100
run_mix EuclidDist_HighPrecision/1000/100
run_mix KMeans_MixedPrecision/0
run_mix KMeans_MixedPrecision/1
run_mix KMeans_MixedPrecision/2
run_mix KMeans_MixedPrecision/3
run_mix KMeans_MixedPrecision/4
run_mix KMeans_HighPrecision/0
run_mix KMeans_HighPrecision/1
run_mix KMeans_HighPrecision/2
run_mix KMeans_HighPrecision/3
run_mix KMeans_HighPrecision/4
cd ..

echo "" >>$ws_dir/results.txt

cd simpsons
gen_mix
run_mix Simpsons_MixedPrecision/10000
run_mix Simpsons_MixedPrecision/100000
run_mix Simpsons_MixedPrecision/1000000
run_mix Simpsons_MixedPrecision/10000000
run_mix Simpsons_MixedPrecision/100000000
run_mix Simpsons_HighPrecision/10000
run_mix Simpsons_HighPrecision/100000
run_mix Simpsons_HighPrecision/1000000
run_mix Simpsons_HighPrecision/10000000
run_mix Simpsons_HighPrecision/100000000
cd ..