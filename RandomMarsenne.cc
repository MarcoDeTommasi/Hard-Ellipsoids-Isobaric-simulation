#ifndef RandomMarsenne_cc
#define RandomMarsenne_cc

#include <iostream>
#include <random>
#include <time.h>
#include <chrono>

using namespace std;
using unidist = uniform_real_distribution<double>;
class RandomMarsenne{
 public:
  RandomMarsenne(){
  }
  double RandM01(){
    unidist xsi(0,1.0);
    random_device rd;
    seed_seq rand_seq{rd()^time(NULL),rd()^time(NULL)};
    mt19937 mt;
    mt.seed(rand_seq);
    return xsi(mt);
  }
};

#endif
