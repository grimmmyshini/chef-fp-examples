FROM ubuntu:jammy

RUN apt -y update
RUN apt -y upgrade
RUN apt install -y build-essential
RUN apt install -y wget cmake gcc git python3 python3-pip curl time
RUN apt install -y clang-13 libclang-13-dev llvm-13-tools llvm-13-dev

WORKDIR /code
RUN git clone https://github.com/google/benchmark -b v1.6.1
RUN cd benchmark/ && \
    cmake -E make_directory "build" && \
    cmake -E chdir "build" cmake -DBENCHMARK_DOWNLOAD_DEPENDENCIES=on -DCMAKE_BUILD_TYPE=Release ../ && \
    cmake --build "build" --config Release
RUN git clone https://github.com/LLNL/adapt-fp
RUN cd adapt-fp && \
    git checkout be4cda8
RUN wget https://github.com/SciCompKL/CoDiPack/archive/refs/tags/v1.9.1.tar.gz
RUN tar -xzvf v1.9.1.tar.gz
RUN rm v1.9.1.tar.gz
RUN pip3 install lit
RUN git clone https://github.com/vgvassilev/clad --depth=1 --single-branch --branch=v1.1
RUN cd clad && \
    mkdir build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release -DLLVM_DIR=/usr/lib/llvm-13 -DClang_DIR=/usr/lib/llvm-13 -DCMAKE_INSTALL_PREFIX=../inst -DLLVM_EXTERNAL_LIT="$(which lit)" && \
    make -j6 && make install -j6

COPY . clad-fp-error-est-examples/

WORKDIR /code/clad-fp-error-est-examples

RUN cd PrintModel && \
    clang++-13 -I$PWD/../../clad/inst/include \
        -I/usr/include/llvm-13 -I/usr/lib/llvm-13/include \
        -fPIC -shared -fno-rtti -Wl,-undefined -Wl,suppress PrintModel.cpp \
        -o libPrintModel.so
RUN cd blackscholes/data && \
    clang++-13 inputGen.c -o inputGen.out && \
    ./inputGen.out 100 100.txt && \
    ./inputGen.out 1000 1000.txt && \
    ./inputGen.out 10000 10000.txt && \
    ./inputGen.out 100000 100000.txt && \
    ./inputGen.out 1000000 1000000.txt
RUN cd kmeans/data/inpuGen/ && \
    clang++-13 datagen.cpp -o datagen.out && \
    ./datagen.out 100 -f && mv 100_34f.txt ../100.txt && \
    ./datagen.out 1000 -f && mv 1000_34f.txt ../1000.txt && \
    ./datagen.out 10000 -f && mv 10000_34f.txt ../10000.txt && \
    ./datagen.out 100000 -f && mv 100000_34f.txt ../100000.txt && \
    ./datagen.out 1000000 -f && mv 1000000_34f.txt ../1000000.txt

WORKDIR /code

RUN echo "export CPLUS_INCLUDE_PATH=\$CPLUS_INCLUDE_PATH:/code/CoDiPack-1.9.1/include:/code/benchmark/include:/code/adapt-fp:/code/clad/include:/code/clad-fp-error-est-examples/PrintModel" >> $HOME/.bashrc
RUN echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/code/clad/build/lib" >> $HOME/.bashrc
RUN echo "export LIBRARY_PATH=\$LIBRARY_PATH:/code/benchmark/build/src" >> $HOME/.bashrc
RUN echo "export CODIPACK_HOME=/code/CoDiPack-1.9.1"