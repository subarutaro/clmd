#include "md.hxx"

MD::MD(){};

MD::~MD(){
  if(r) free(r);
  if(v) free(v);
  if(f) free(f);
}

void MD::Initialize(const int _nmol,const DorF _dens,const DorF _temp){
  cout << "# ---construct MD---" << endl;

  nmol = _nmol;

  dt   = 0.000930f;
  dth  = 0.5f * dt;
  dens = _dens;
  temp = _temp;

  cs = powf((DorF)nmol/dens,1./3.);
  csh = 0.5f * cs;

  if(csh > 5.f) rcut = 5.f;
  else          rcut = csh;

  r = (DorF*)calloc(3*nmol,sizeof(DorF));
  assert(r != NULL);
  rx = &r[0];
  ry = &r[nmol];
  rz = &r[2*nmol];

  v = (DorF*)calloc(3*nmol,sizeof(DorF));
  assert(v != NULL);
  vx = &v[0];
  vy = &v[nmol];
  vz = &v[2*nmol];

  f = (DorF*)calloc(3*nmol,sizeof(DorF));
  assert(f != NULL);
  fx = &f[0];
  fy = &f[nmol];
  fz = &f[2*nmol];
}

void MD::SetCoorFCC(){
  int ndim=1;
  while(nmol > 4*ndim*ndim*ndim) ndim++;
  DorF unit = cs / (DorF)ndim;
  {
    int i = 0;
    for(int z=0;z<ndim;z++){
      for(int y=0;y<ndim;y++){
	for(int x=0;x<ndim;x++){
	  rx[i] = unit * x;
	  ry[i] = unit * y;
	  rz[i] = unit * z;
	  i++;if(i==nmol) break;
	  rx[i] = unit * (x+0.5);
	  ry[i] = unit * (y+0.5);
	  rz[i] = unit * z;
	  i++;if(i==nmol) break;
	  rx[i] = unit * x;
	  ry[i] = unit * (y+0.5);
	  rz[i] = unit * (z+0.5);
	  i++;if(i==nmol) break;
	  rx[i] = unit * (x+0.5);
	  ry[i] = unit * y;
	  rz[i] = unit * (z+0.5);
	  i++;if(i==nmol) break;
	}
	if(i==nmol) break;
      }
      if(i==nmol) break;
    }
  }
  for(int i=0;i<nmol;i++){
    rx[i] -= 0.5*cs;
    ry[i] -= 0.5*cs;
    rz[i] -= 0.5*cs;
  }
}

void MD::SetVelRandom(){
  for(int i=0;i<nmol;i++){
    vx[i] = (DorF)rand()/(DorF)RAND_MAX;
    vy[i] = (DorF)rand()/(DorF)RAND_MAX;
    vz[i] = (DorF)rand()/(DorF)RAND_MAX;
  }
  KillMomentum();
  VelocityScaling();
}

void MD::VelocityScaling(){
  CalcKineticEnergy();
  const DorF scale = sqrt(1.5*nmol*temp/kin);
  for(int i=0;i<nmol;i++){
    vx[i] *= scale;
    vy[i] *= scale;
    vz[i] *= scale;
  }
  kin *= scale*scale;
  KillMomentum();
}

void MD::KillMomentum(){
  DorF mom[3] = {0.f,0.f,0.f};
  for(int i=0;i<nmol;i++){
    mom[0] += vx[i];
    mom[1] += vy[i];
    mom[2] += vz[i];
  }
  mom[0] /= (DorF)nmol;
  mom[1] /= (DorF)nmol;
  mom[2] /= (DorF)nmol;
  for(int i=0;i<nmol;i++){
    vx[i] -= mom[0];
    vy[i] -= mom[1];
    vz[i] -= mom[2];
  }
}

void MD::CalcKineticEnergy(){
  kin = 0.0;
  for(int i=0;i<nmol;i++){
    kin += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
  }
  kin *= 0.5f;
}

void MD::IntegrateCoor(){
  for(int i=0;i<nmol;i++){
    rx[i] += vx[i] * dt;
    ry[i] += vy[i] * dt;
    rz[i] += vz[i] * dt;
    if(rx[i] <= -csh) rx[i] += cs; if(rx[i] > csh) rx[i] -= cs;
    if(ry[i] <= -csh) ry[i] += cs; if(ry[i] > csh) ry[i] -= cs;
    if(rz[i] <= -csh) rz[i] += cs; if(rz[i] > csh) rz[i] -= cs;
  }
}

void MD::IntegrateVel(){
  for(int i=0;i<nmol;i++){
    vx[i] += fx[i] * 0.5*dt;
    vy[i] += fy[i] * 0.5*dt;
    vz[i] += fz[i] * 0.5*dt;
  }
}

void MD::CalcForce(){
  for(int i=0;i<3*nmol;i++) f[i] = 0.0;
  const DorF csh = 0.5*cs;
  for(int i=0;i<nmol;i++){
    const DorF rxi = rx[i];
    const DorF ryi = ry[i];
    const DorF rzi = rz[i];
    DorF fxi = 0.0;
    DorF fyi = 0.0;
    DorF fzi = 0.0;
    for(int j=0;j<nmol;j++){
      if(i==j) continue;
      DorF dx = rxi - rx[j];
      DorF dy = ryi - ry[j];
      DorF dz = rzi - rz[j];
      if(dx <= -csh) dx += cs; if(dx > csh) dx -= cs;
      if(dy <= -csh) dy += cs; if(dy > csh) dy -= cs;
      if(dz <= -csh) dz += cs; if(dz > csh) dz -= cs;

      const DorF r02 = dx*dx + dy*dy + dz*dz;
      const DorF r02i = 1.0 / r02;
      const DorF r06i = r02i * r02i *r02i;

      const DorF ftmp = 48.0 * r06i * (r06i - 0.5) * r02i;
      fxi += ftmp * dx;
      fyi += ftmp * dy;
      fzi += ftmp * dz;
    }
    fx[i] = fxi;
    fy[i] = fyi;
    fz[i] = fzi;
  }
}

void MD::CalcPotentialEnergy(){
  pot = 0.0;
  const DorF csh = 0.5*cs;
  for(int i=0;i<nmol;i++){
    const DorF rxi = rx[i];
    const DorF ryi = ry[i];
    const DorF rzi = rz[i];
    DorF fxi = 0.0;
    DorF fyi = 0.0;
    DorF fzi = 0.0;
    for(int j=0;j<nmol;j++){
      if(i==j) continue;
      DorF dx = rxi - rx[j];
      DorF dy = ryi - ry[j];
      DorF dz = rzi - rz[j];
      if(dx <= -csh) dx += cs; if(dx > csh) dx -= cs;
      if(dy <= -csh) dy += cs; if(dy > csh) dy -= cs;
      if(dz <= -csh) dz += cs; if(dz > csh) dz -= cs;

      const DorF r02 = dx*dx + dy*dy + dz*dz;
      const DorF r02i = 1.0 / r02;
      const DorF r06i = r02i * r02i *r02i;

      pot += 2.0 * r06i * (r06i - 1.0);
    }
  }
}

void MD::OutputCDV(const string filename){
  ofstream ofs(filename.c_str());
  ofs << "'box_sx=" << -csh << ",box_ex="<< csh;
  ofs << ",box_sy=" << -csh << ",box_ey="<< csh;
  ofs << ",box_sz=" << -csh << ",box_ez="<< csh;
  ofs << endl;
  for(int i=0;i<nmol;i++){
    ofs << ' ' << i     << " " << "0";
    ofs << ' ' << rx[i] << ' ' << ry[i] << ' ' << rz[i];
    ofs << ' ' << vx[i] << ' ' << vy[i] << ' ' << vz[i];
    ofs << ' ' << fx[i] << ' ' << fy[i] << ' ' << fz[i];
    ofs << endl;
  }
}

void MD::DisplayEnergies(ostream &os){
  static int s = 0;
  os << s++ << ' ' << pot << ' ' << kin << ' ' << pot+kin << ' ' << (pot+kin - tot)/tot << endl;
}

void MD::DisplayConditions(ostream &os){
  os << "# Number of molecules:\t" << nmol   << endl;
  os << "# Density:\t"        << dens        << endl;
  os << "# Temperature:\t"    << temp        << endl;
  os << "# Cell size:\t"      << cs          << endl;
  os << "# Cutoff radii:\t"   << rcut        << endl;
  os << "# Delta time:\t"     << dt << ' ' << dth << endl;
  os << endl;
}
