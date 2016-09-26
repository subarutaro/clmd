#include <iostream>
#include <sys/time.h>

class Timer{
private:
  double elapsed_time;
  double beg,end;
  void gettimeofday();
public:
  Timer(){
    flush();
  }
  ~Timer(){};
  void flush(){
    elapsed_time = 0.0;
  }
  void beg(){
    
  }
};
