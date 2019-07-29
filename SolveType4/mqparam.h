#ifndef PROJECT_MQPARAM_H
#define PROJECT_MQPARAM_H

#include <iostream>

// structure to store one term
struct term {
  uint32_t data[3];
};

// structure to store a MQ polynomial
struct poly {
  // term array
  term *p;
  int length;
};

#define M 67
#define N 67
#define LEN 10001
#define VARIABLE_NUM 10
#define EQUATION_NUM 22

#endif //PROJECT_MQPARAM_H
