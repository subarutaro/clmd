#ifndef HXX_CLMD
#define HXX_CLMD

#include "md.hxx"
#include "clmanager.hxx"

class clMD : public MD{
private:
  clManager *clmanager;
public:
  cl_mem r_dev,v_dev,f_dev;
  clMD();
  ~clMD();
  void Initialize(const int,const DorF,const DorF);
  void CalcForce();
};


#endif
