#!/bin/bash
if [ -z "${ws_dir}" ]; then
    export ws_dir=/code
    echo "Setting ws_dir"
fi
OUTPUT_PROCESSING_SCRIPT=$ws_dir/chef-fp-examples/process_output.py

compile() {
    clang++-13 \
        -Xclang -add-plugin -Xclang clad -Xclang -load -Xclang clad.so \
        -Xclang -plugin-arg-clad -Xclang -fcustom-estimation-model \
        -Xclang -plugin-arg-clad -Xclang $ws_dir/chef-fp-examples/PrintModel/libPrintModel.so \
        -lstdc++ -lm -std=c++11 -O3 \
        "$@"
}

runmain() {
    cd $1
    compile clad-main.cpp -o cheffp-bin
     echo "[[[[[[[[[[[[[[[[[[[[[[[[[[[[ Running: $1 ]]]]]]]]]]]]]]]]]]]]]]]]]]]]"
    ./cheffp-bin ${@:2}
    cd - > /dev/null
}

runmain arclength
runmain kmeans -i data/100.txt
runmain simpsons
runmain HPCCG 20 30 10
runmain blackscholes data/100.txt
