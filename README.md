# CHEF-FP examples repository

This repository contains benchmarks for CHEF-FP comparing it against the state of the art ADAPT-FP tool.

## Building a docker image for setting up CHEF-FP

To build the docker image run:

```bash
docker build -t "cheffp:Dockerfile" .
```

This installs CHEF-FP, ADAPT-FP along with their dependencies.

## Using the docker image

Run the following command to open a bash terminal in a container:

```bash
docker run -it --rm cheffp:Dockerfile /bin/bash
```

## Run the CHEF-FP vs ADAPT benchmarks in the docker image

To run the benchmarks run the docker image and execute:

```bash
cd /code/clad-fp-error-est-examples
./runner.sh
```

The results of the benchmark will be in the `/code/results.txt` file.

## Using CHEF-FP for floating point error estimation

Lets take the following code as an example:

```c++
double f(double x) {
    double y = x * x;
    return y;
}

int main() {
    f(2.3);
}
```

Include the `ErrorFunc.h` header in `PrintModel/` directory.
Then, Mark the function you want to run error estimation on using
`clad::error_estimate` and run the generated function using `execute`.
In this example we use the error estimation model used in the paper so we
can print the mixed precision analysis using `clad::printErrorReport()`.
So finally it will look like something like this:

```c++
#include "clad/Differentiator/Differentiator.h"
#include "ErrorFunc.h"

double f(double x) {
    double y = x * x;
    return y;
}

int main() {
    // Generate the floating point error estimation code for 'f'
    auto df = clad::estimate_error(f);

    // Declare the necessary variables and call the generated code
    double dx = 0, final_error = 0;
    df.execute(2.3, &dx, final_error);

    // Print the mixed precision configuration
    clad::printErrorReport();
}
```

To compile use the following command:
```bash
clang++-13 \
        -Xclang -add-plugin -Xclang clad -Xclang -load -Xclang clad.so \
        -Xclang -plugin-arg-clad -Xclang -fcustom-estimation-model \
        -Xclang -plugin-arg-clad -Xclang /code/clad-fp-error-est-examples/PrintModel/libPrintModel.so \
        -lstdc++ -lm -std=c++11 -O3 \
        source.cpp
```

And run it:
```bash
./a.out
```

## Create your own custom models

To create your own custom models go through `/code/clad/demos/ErrorEstimation/CustomModel/README.md` or it's [online version](https://github.com/vgvassilev/clad/blob/v1.1/demos/ErrorEstimation/CustomModel/README.md)