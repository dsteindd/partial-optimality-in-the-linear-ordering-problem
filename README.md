# Code for the Paper "Partial Optimality in the Linear Ordering Problem"

This is the code for "Partial Optimality in the Linear Ordering Problem", accepted at ICML'24.

## Requirements
You need to have [Gurobi](https://gurobi.com) installed. 

## Data Preparation
Run ```prepare_lolib.sh```. It just unpacks lolib_2010.zip in the ```./data``` directory.

## How to build?

```shell
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=RELEASE ..
make
```

## How to run?
Depending on the instances and experiments you want to run, execute one of the following commands.

### Joint Application (without solving ILP of smaller problems)

#### Synthetic Instances

```shell
./synthetic-problem n --maxElements 200 --alphas 0.4 0.65 0.7 1.0
./synthetic-problem a --ns 20 40 60 80 100 140 200
```
The results can be viewed at ```./results/syntethic/alphas``` and ```./results/syntethic/ns```.

#### LOLIB Instances

```shell
./lolib-problem p --problem IO
./lolib-problem p --problem SGB
./lolib-problem p --problem Spec
```
The results can be viewed at ```./results/lolib/persistency```

### Joint Application (with solving ILP of smaller problems)
#### Synthetic Instances

```shell
./synthetic-problem presolve-and-ilp-n --maxElements 200 --alphas 0.4 0.65 0.7 1.0
./synthetic-problem presolve-and-ilp-a --ns 20 40 60 80 100 140 200
```
The results can be viewed at ```./results/syntethic/ilp-solver-alphas``` and ```./results/syntethic/ilp-solver-ns```.

#### LOLIB Instances

```shell
./lolib-problem pi --problem IO
./lolib-problem pi --problem Spec
```
The results can be viewed at ```./results/lolib/persistency-and-ilp```


### Individual Application of each Condition
#### Synthetic Instances

```shell
./synthetic-problem-individual-conditions n --maxElements 200 --alphas 0.4 0.65 0.7 1.0
./synthetic-problem-individual-conditions a --ns 20 40 60 80 100 140 200
```
The results can be viewed at ```./results/synthetic-individual/alphas``` and ```./results/synthetic-individual/ns```.

```shell
./lolib-problem-individual-conditions --problem IO
./lolib-problem-individual-conditions --problem Spec
```
The results can be viewed at ```./results/lolib/individual```.