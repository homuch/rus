# RUS -- Repeat-Until-Success quantum circuit synthesis

## BUILD

You will need the following libraries installed on your system: 
1. Boost 1.48
-- program_options 
-- chrono
-- timer
-- system
2. The GNU Multiple Precision Arithmetic Library (gmp and gmpxx)
3. The GNU MPFR Library (mpfr)

Also C++ compiler supporting C++11 is necessary.
Information about program use available through --help option.

## ABOUT 
The program code based on results of https://arxiv.org/abs/1409.3552. 
In addition to Boost, The GNU Multiple Precision Arithmetic Library, The GNU MPFR Library the library 
mpfr::real by Christian Schneider <software(at)chschneider(dot)eu> is used for high precision
In addition, a significant portion of the code leverages from [SQCT](https://github.com/vadym-kl/sqct).

## DIRECTORY STRUCTURE 
* sk -- implementation of the Solovay-Kitaev algorithm
* es -- exact synthesis algorithm
* theory -- numerical proof of result from arXiv:1206.5236, tests of exact synthesis algorithm 
* appr -- optimal round off of unitaries
* rus -- the rus implementation

## USAGE

### Build Single Gate

(todo)

### Build Gate Database
Add `--database`

For example, 

```
./rusSyn -O output_folder -F 5 -E 1e-5 -C g-count -P pi/32 --qubit-name q[0] --ancil-name a[0] --database
```

will generate Rz gate from $+\pi$ to $-\pi$ per $\pi/32$ (Centered at 0). The output files name will be `${output_folder}/out${i}_pi|32_1e-5`, where $i$ is integer (which means $\theta=i*\pi/32$).

(The program may die for high effort and precision currently)