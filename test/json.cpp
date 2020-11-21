#include <iostream>
#include <vector>
#include "../include/nlohmann_json_cmake_fetchcontent/include/nlohmann/json.hpp"
using json = nlohmann::json;


int main() {
  json j;
  j["list"] = {1, 2, 3};
  std::cout << j.dump(4) << std::endl;

  std::cout << j["list"].get<std::vector<int>>().size() << std::endl;
  
  std::vector<int> hoge = j["list"].get<std::vector<int>>();
  std::cout << hoge[0] << std::endl;
  return 0;  
}
