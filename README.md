# Parallelization of the Wang-Landau method by MPI
## How to build
You can build sample.cpp with

~~~shell-session
$ mkdir build
$ cd build
$ cmake ..
$ make
~~~

When you use this package you need to create directory named log on your running
directory.

~~~shell-session
$ mkdir log
~~~

## Requirements
### Modules
This package needs OpenMPI and Json for Modern C++.
The latter is included as the git submodule.
Code for sample model depends on Eigen/Dense.
### API of your own model
```c++
double YourModel::Propose(std::mt19937 &engine);
void YourModel::Update();
void YourModel::ExchangeConfig(int partner, MPI_Comm local_comm);
void YourModel::StoreLog(std::ofstream *ofs_ptr);
void YourModel::SetFromLog(std::ifstream &ifs);
double YourModel::val();
size_t YourModel::sweeps();
void YourModel::set_val(double energy_new);
```
Implementation examples are in `model_sample/ferro_ising.hpp`.

## License
The code is licensed under the [MIT License](https://opensource.org/licenses/MIT):

Copyright &copy; 2020 Chihiro Kondo

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Thanks
- [Replica-Exchange Wang-Landau sampling tutorial lectures at IX Brazilian Meeting on Simulational Physics, 2017](https://github.com/yingwaili/bmsp2017)
- [JSON for Modern C++](https://github.com/nlohmann/json#)
