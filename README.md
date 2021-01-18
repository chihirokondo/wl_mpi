# MPI parallelization of the Wang-Landau method
## How to use
The single file `inlcude/wl_mpi.hpp` is required to use this package.
When you use this package you need to create directory named `log` on your running
directory.

~~~shell-session
$ mkdir log
~~~

### How to build `sample.cpp`
You can build `sample.cpp` with

~~~shell-session
$ mkdir build
$ cd build
$ cmake ..
$ make
~~~

## Requirements
### Modules
This package needs OpenMPI and Json for Modern C++.
The latter is embedded as the git submodule.
Code for sample model depends on Eigen/Dense.
### Member functions you must prepare
If you would like to apply this package to your model, you need to prepare the following member functions in your model class (`YourModel`).
```c++
// Just propose not update.
double YourModel::Propose(std::mt19937 &engine);

// Accept propose and update the state.
void YourModel::Update();

// Exchange informations with given partner in given communicator.
// e.g. energy and configuration (case of estimating density of state).
void YourModel::Exchange(int partner, MPI_Comm local_comm);

// Store intermediate state in the log file, which is necessary to restart the execution.
void YourModel::StoreLog(std::ofstream *ofs_ptr);

// Read intermediate state from the log file.
void YourModel::SetFromLog(std::ifstream &ifs);

// Return current value which belongs to searching space by the Wang-Landau method.
// e.g. energy (case of estimating density of state).
double YourModel::val();

// Return definition of 1MCS on your model.
size_t YourModel::sweeps();
```
Implementation examples are in `model_sample/ferro_ising.hpp`.

## License
The code is licensed under the [MIT License](https://opensource.org/licenses/MIT):

Copyright &copy; 2020-2021 Chihiro Kondo

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Thanks
- [Replica-Exchange Wang-Landau sampling tutorial lectures at IX Brazilian Meeting on Simulational Physics, 2017](https://github.com/yingwaili/bmsp2017)
- [JSON for Modern C++](https://github.com/nlohmann/json#)
