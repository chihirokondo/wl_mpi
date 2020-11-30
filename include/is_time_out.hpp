#ifndef WANGLANDAU_IS_TIME_OUT_H_
#define WANGLANDAU_IS_TIME_OUT_H_


#include <iostream>
#include <chrono>


class IsTimeOut {
 public:
  IsTimeOut(double timelimit_secs)
      : start_(std::chrono::system_clock::now()),
        timelimit_secs_(timelimit_secs) {}
  bool operator()() const;
 private:
  const std::chrono::system_clock::time_point start_;
  const double timelimit_secs_;
};


inline bool IsTimeOut::operator()() const {
  // true = (time is out), false = (time is "in").
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  bool is_time_out =
      std::chrono::duration_cast<std::chrono::seconds>(end-start_).count() >=
      timelimit_secs_;
  return is_time_out;
}


#endif // WANGLANDAU_IS_TIME_OUT_H_
