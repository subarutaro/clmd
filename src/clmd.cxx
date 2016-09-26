#include "clmd.hxx"

clMD::clMD(){
  clParams params("input.txt");
  clmanager = new clManager(params);
};
clMD::~clMD(){};

void clMD::Initialize(const int _nmol,const DorF _dens,const DorF _temp){
  MD::Initialize(_nmol,_dens,_temp);

  r_dev = clCreateBuffer(clmanager->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(DorF)*3*nmol, r, &clmanager->err);
  f_dev = clCreateBuffer(clmanager->context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(DorF)*3*nmol, f, &clmanager->err);

  clSetKernelArg(clmanager->kernel,0,sizeof(cl_mem),&r_dev);
  clSetKernelArg(clmanager->kernel,1,sizeof(cl_mem),&f_dev);
  clSetKernelArg(clmanager->kernel,2,sizeof(DorF),&rcut);
  clSetKernelArg(clmanager->kernel,3,sizeof(DorF),&cs);
  clSetKernelArg(clmanager->kernel,4,sizeof(DorF),&nmol);

  clmanager->CreateCommandQueue();
}

void clMD::CalcForce(){
  clmanager->err = clEnqueueWriteBuffer(clmanager->queue,r_dev,CL_TRUE, 0,sizeof(float)*3*nmol,r, 0,NULL,NULL);
  clmanager->EnqueueNDRangeKernel(nmol,nmol/4);
  clmanager->err = clEnqueueReadBuffer(clmanager->queue,f_dev,CL_TRUE, 0,sizeof(float)*3*nmol,f, 0,NULL,NULL);
  //for(int i=0;i<nmol;i++) cout << i << ' ' << f[i] << ' ' << f[i+nmol] << ' ' << f[i+2*nmol] << endl;
}
