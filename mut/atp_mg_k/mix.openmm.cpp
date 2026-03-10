#include <string>
#include <vector>
#include <set>
#include <type_traits>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <cstring>

// Zlib Required
#include "zlib.h"
#include "OpenMM.h"

// Compile option: g++ -O3 mix.openmm.cpp -lz -lOpenMM

// Path to openmm plugin librarys (libOpenMMCPU.so, libOpenMMCUDA.so etc.)
#ifndef XCG_OPENMM_PLUGPATH
#define XCG_OPENMM_PLUGPATH "/opt/openmm-8.1.2/exe/lib/plugins"
#endif

// MD simulation variables:
//   Cutoff distance in nm
const double cutoff_lj = 1.00000;
//   Total timesteps
const int ntotalSteps = 1000000;
//   Output frequency
const int output_freq = 1000;
//   Temp (K)
const double temp = 300.0;
//   Pres (bar)
const double pres = 1.00000;
//   Timestep (ps)
const double timestep = 0.00200;

// MD simulation macros:
// Define XCG_OPENMM_PME for PME charge calc
//#define XCG_OPENMM_PME

// Define XCG_OPENMM_NPT for Pressure controling brodstat
//#define XCG_OPENMM_NPT

// Define XCG_OPENMM_NOSHAKE for disabling Shake
//#define XCG_OPENMM_NOSHAKE

// Define XCG_OPENMM_PLUMED for turning on plumed forces from plumed.dat (Need openmm-plumed plugin and compile with extra -lOpenMMPlumed option)
//#define XCG_OPENMM_PLUMED

// Define XCG_OPENMM_SOFTCORE for turning on FEP the soft core potential (-DXCG_OPENMM_SOFTCORE=SOFTCORE or FEP)
//#define XCG_OPENMM_SOFTCORE FEP
//#define XCG_OPENMM_SOFTCORE SOFTCORE

#ifdef XCG_OPENMM_PLUMED
#include "PlumedForce.h"
#endif


inline void fwrite_int32(FILE *fp, uint32_t i) {
  fwrite(&i, sizeof(uint32_t), 1, fp);
}

std::string dump_file(std::string file) {
  FILE *fp = fopen(file.c_str(), "r");
  std::string ret;
  char buff[1024];
  while(fgets(buff, 1024, fp) != nullptr) {
    ret += buff;
  }
  fclose(fp);
  return ret;
}

void write_header(FILE *fp, int natoms, std::string remarks, int nevery = 10000, float dt = 0.002, int nstep = 0) {
  char title_string[200];
  fwrite_int32(fp, 84);
  strcpy(title_string, "CORD");
  fwrite(title_string, 4, 1, fp);
  fwrite_int32(fp, 0);
  fwrite_int32(fp, nstep);
  fwrite_int32(fp, nevery);
  fwrite_int32(fp, nstep);
  fwrite_int32(fp, 0);
  fwrite_int32(fp, 0);
  fwrite_int32(fp, 0);
  fwrite_int32(fp, 0);
  fwrite_int32(fp, 0);
  fwrite(&dt, sizeof(float), 1, fp);
  fwrite_int32(fp, 1);
  fwrite_int32(fp, 0);
  fwrite_int32(fp, 0);
  fwrite_int32(fp, 0);
  fwrite_int32(fp, 0);
  fwrite_int32(fp, 0);
  fwrite_int32(fp, 0);
  fwrite_int32(fp, 0);
  fwrite_int32(fp, 0);
  fwrite_int32(fp, 24);
  fwrite_int32(fp, 84);
  fwrite_int32(fp, 164);
  fwrite_int32(fp, 2);
  strncpy(title_string, remarks.c_str(), 80);
  title_string[79] = '\0';
  fwrite(title_string, 80, 1, fp);
  fwrite(title_string, 80, 1, fp);
  fwrite_int32(fp, 164);
  fwrite_int32(fp, 4);
  fwrite_int32(fp, natoms);
  fwrite_int32(fp, 4);
  fflush(fp);
}

void write_box(FILE *fp, double dimx, double dimy, double dimz) {
  double dim[6];
  dim[0] = dimx * 10.0;
  dim[2] = dimy * 10.0;
  dim[5] = dimz * 10.0;
  dim[1] = dim[3] = dim[4] = 0.0;
  uint32_t out_integer = 48;
  fwrite_int32(fp, out_integer);
  fwrite(dim, out_integer, 1, fp);
  fwrite_int32(fp, out_integer);
  fflush(fp);
}

#define NFILE_POS 8L
#define NSTEP_POS 20L

void write_frame(FILE *fp, int natoms, const std::vector<OpenMM::Vec3> &x, int nframes, int nevery) {
  uint32_t out_integer = natoms * sizeof(float);
  fwrite_int32(fp, out_integer);
  for (int i = 0; i < natoms; i++) {
    float v = x[i][0] * 10.0;
    fwrite((void*)&v, sizeof(float), 1, fp);
  }
  fwrite_int32(fp, out_integer);
  fwrite_int32(fp, out_integer);
  for (int i = 0; i < natoms; i++) {
    float v = x[i][1] * 10.0;
    fwrite((void*)&v, sizeof(float), 1, fp);
  }
  fwrite_int32(fp, out_integer);
  fwrite_int32(fp, out_integer);
  for (int i = 0; i < natoms; i++) {
    float v = x[i][2] * 10.0;
    fwrite((void*)&v, sizeof(float), 1, fp);
  }
  fwrite_int32(fp, out_integer);

  out_integer = nframes;
  fseek(fp, NFILE_POS, SEEK_SET);
  fwrite_int32(fp, out_integer);
  out_integer = nframes * nevery;
  fseek(fp, NSTEP_POS, SEEK_SET);
  fwrite_int32(fp, out_integer);
  fseek(fp, 0, SEEK_END);
}

std::set<int> get_linspace(int first, int last) {
  std::set<int> ret;
  for(int i = first; i <= last; i++) {
    ret.insert(i);
  }
  return ret;
}

std::set<int> vec2set_int(const std::vector<int> &v) {
  std::set<int> ret;
  for(auto &p : v) {
    ret.insert(p);
  }
  return ret;
}

void writeXYZ(std::string file, const std::vector<OpenMM::Vec3> &coord) {
  FILE *fp = fopen(file.c_str(), "a");
  fprintf(fp, "%d\n\n", (int)coord.size());
  for (int i = 0; i < (int)coord.size(); i++) {
    fprintf(fp, "He %.5f %.5f %.5f\n", coord[i][0] * 10.0, coord[i][1] * 10.0, coord[i][2] * 10.0);
  }
  fclose(fp);
}
std::vector<std::vector<int>> cat_vec_int(const std::vector<std::vector<int>> &v1, const std::vector<std::vector<int>> &v2) {
  std::set<std::set<int>> unique_set;
  std::vector<std::vector<int>> unique_vec;
  for (auto &p : v1) {
    unique_set.insert({p[0], p[1]});
  }
  for (auto &p : v2) {
    unique_set.insert({p[0], p[1]});
  }
  for (auto &p : unique_set) {
    std::vector<int> v;
    for (auto &q : p) {
      v.push_back(q);
    }
    unique_vec.push_back(v);
  }
  return unique_vec;
}

std::vector<double> interpol_exc(std::vector<double> param, double min, double max, int nint) {
  double exc = param[0];
  double sig = param[1];
  double gridSize = (max - min) / (nint - 1);
  double A = 4.0 * exc * pow(sig, 12);
  std::vector<double> pot(nint);
  for (int i = 0; i < nint; i++) {
    double r = min + gridSize * i;
    double rsq = r * r;
    double r6 = rsq * rsq * rsq;
    pot[i] = A / r6 / r6;
  }
  double offset = pot.back();
  for (int i = 0; i < nint; i++) {
    pot[i] -= offset;
  }
  return pot;
}

std::vector<double> interpol_lj(std::vector<double> param, double min, double max, int nint) {
  double eps = param[0];
  double sig = param[1];
  double gridSize = (max - min) / (nint - 1);
  double A = 5.0 * eps * pow(sig, 12);
  double B = 6.0 * eps * pow(sig, 10);
  std::vector<double> pot(nint);
  for (int i = 0; i < nint; i++) {
    double r = min + gridSize * i;
    double rsq = r * r;
    double r4 = rsq * rsq;
    double r6 = r4 * rsq;
    pot[i] = A / r6 / r6 - B / r4 / r6;
  }
  double offset = pot.back();
  for (int i = 0; i < nint; i++) {
    pot[i] -= offset;
  }
  return pot;
}

std::vector<double> interpol_lj_coul(std::vector<double> param, double min, double max, int nint) {
  double C12 = param[0];
  double C6 = param[1];
  double q1 = param[2];
  double q2 = param[3];
  double gridSize = (max - min) / (nint - 1);
  double A = 4.184 * 332.05221729 * 0.1;
  std::vector<double> pot(nint);
  for (int i = 0; i < nint; i++) {
    double r = min + gridSize * i;
    double rsq = r * r;
    double r6 = rsq * rsq * rsq;
    pot[i] = C12 / r6 / r6 - C6 / r6 + A * q1 * q2 / r;
  }
  double offset = pot.back();
  for (int i = 0; i < nint; i++) {
    pot[i] -= offset;
  }
  return pot;
}

std::vector<double> interpol_zero(int nint) {
  return std::vector<double>(nint, 0.0);
}
std::vector<double> interpol_minus(const std::vector<double> &v1, const std::vector<double> &v2) {
  std::vector<double> ret;
  for (int i = 0; i < v1.size(); i++) {
    ret.push_back(v1[i] - v2[i]);
  }
  return ret;
}

std::vector<double> interpol_dw(std::vector<double> param, std::vector<int> knm, double min, double max, int nint) {
  int nparam = param.size() / 2;
  int nts = (nparam - 1) / 2;
  int K = knm[0], N = knm[1], M = knm[2];
  std::vector<double> eps(nparam), sig(nparam);
  for (int i = 0; i < nparam; i++) {
    eps[i] = param[2 * i];
    sig[i] = param[2 * i + 1];
  }
  double gridSize = (max - min) / (nint - 1);
  double A_1 = eps[0] * pow(sig[0], 2 * K);
  double B_1 = 2 * eps[0] * pow(sig[0], K);
  double A_n = 5 * eps.back() * pow(sig.back(), 12);
  double B_n = 6 * eps.back() * pow(sig.back(), 10);
  std::vector<std::vector<double>> tmp1(nts, std::vector<double>(5)), tmp2(nts, std::vector<double>(5));
  for (int i = 0; i < nts; i++) {
    int i1 = 2 * i;
    int its = 2 * i + 1;
    int i2 = 2 * i + 2;
    double e1 = eps[i1];
    double ets = eps[its];
    double e2 = eps[i2];
    double s1 = sig[i1];
    double sts = sig[its];
    double s2 = sig[i2];
    tmp1[i][0] = s1;
    tmp1[i][1] = sts;
    tmp1[i][2] = e1 - ets;
    tmp1[i][3] = s1 - sts;
    tmp1[i][4] = ets;
    tmp2[i][0] = sts;
    tmp2[i][1] = s2;
    tmp2[i][2] = ets - e2;
    tmp2[i][3] = s2 - sts;
    tmp2[i][4] = ets;
  }
  std::vector<double> pot(nint);
  for (int i = 0; i < nint; i++) {
    double r = min + gridSize * i;
    if (r <= sig[0]) {
      double rk = pow(r, K);
      double rk2 = rk * rk;
      pot[i] = A_1 / rk2 - B_1 / rk;
    } else if (r >= sig.back()) {
      double rsq = r * r;
      double r4 = rsq * rsq;
      double r6 = r4 * rsq;
      pot[i] = A_n / r6 / r6 - B_n / r4 / r6;
    } else {
      for (int j = 0; j < nts; j++) {
        if (r >= tmp1[j][0] && r <= tmp1[j][1]) {
          double R1 = (r - tmp1[j][1]) / tmp1[j][3];
          double R1N = pow(R1, N);
          double R12N = R1N * R1N;
          pot[i] = tmp1[j][2] * (R12N - 2 * R1N) - tmp1[j][4];
          continue;
        }
        if (r >= tmp2[j][0] && r <= tmp2[j][1]) {
          double R2 = (r - tmp2[j][0]) / tmp2[j][3];
          double R22 = R2 * R2;
          double R22M = pow(R2, M);
          pot[i] = tmp2[j][2] * M * R2 / (R22M + M - 1) - tmp2[j][4];
          continue;
        }
      }
    }
  }
  double offset = pot.back();
  for (int i = 0; i < nint; i++) {
    pot[i] -= offset;
  }
  return pot;
}

/*
std::vector<double> interpol_dwbond(std::vector<double> param, double delta, double min, double max, int nint) {
  int nstates = param.size() / 2;
  double gridSize = (max - min) / (nint - 1);
  std::vector<double> eng_states(nstates);
  std::vector<double> pot(nint);
  for (int i = 0; i < nint; i++) {
    double r = min + gridSize * i;
    for (int j = 0; j < nstates; j++) {
      double k = param[j * 2];
      double r0 = param[j * 2 + 1];
      double dr = r - r0;
      eng_states[j] = k * dr * dr;
    }
    sort(eng_states.begin(), eng_states.end());
    double e0 = *(eng_states.begin());
    double e1 = *(eng_states.begin() + 1);
    pot[i] = 0.5 * (e1 + e0) - sqrt(pow(e0 - e1, 2) / 4.0 + delta * delta);
  }
  double offset = *(min_element(pot.begin(), pot.end()));
  for (int i = 0; i < nint; i++) {
    pot[i] -= offset;
  }
  return pot;
}

std::vector<double> interpol_dwangle(std::vector<double> param, double delta, double min, double max, int nint) {
  int nstates = param.size() / 3;
  double gridSize = (max - min) / (nint - 1);
  std::vector<double> eng_states(nstates);
  std::vector<double> pot(nint);
  for (int i = 0; i < nint; i++) {
    double r = min + gridSize * i;
    for (int j = 0; j < nstates; j++) {
      double k = param[j * 3];
      double d0 = param[j * 3 + 1];
      double t0 = param[j * 3 + 2];
      double dr = r - t0;
      eng_states[j] = k * (1.0 - exp(-0.5 * dr * dr / d0 / d0));
    }
    sort(eng_states.begin(), eng_states.end());
    double e0 = *(eng_states.begin());
    double e1 = *(eng_states.begin() + 1);
    pot[i] = 0.5 * (e1 + e0) - sqrt(pow(e0 - e1, 2) / 4.0 + delta * delta);
  }
  double offset = *(min_element(pot.begin(), pot.end()));
  for (int i = 0; i < nint; i++) {
    pot[i] -= offset;
  }
  return pot;
}

std::vector<double> interpol_dwdihedral(std::vector<double> param, double delta, double min, double max, int nint) {
  int nstates = param.size() / 3;
  double gridSize = (max - min) / (nint - 1);
  std::vector<double> eng_states(nstates);
  std::vector<double> pot(nint);
  for (int i = 0; i < nint; i++) {
    double r = min + gridSize * i;
    for (int j = 0; j < nstates; j++) {
      double k = param[j * 2];
      double t0 = param[j * 2 + 1];
      double dr = r - t0;
      eng_states[j] = k * (cos(dr) - 1.0);
    }
    sort(eng_states.begin(), eng_states.end());
    double e0 = *(eng_states.begin());
    double e1 = *(eng_states.begin() + 1);
    pot[i] = 0.5 * (e1 + e0) - sqrt(pow(e0 - e1, 2) / 4.0 + delta * delta);
  }
  double offset = *(min_element(pot.begin(), pot.end()));
  for (int i = 0; i < nint; i++) {
    pot[i] -= offset;
  }
  return pot;
}
*/

template <class ValType,
          typename std::enable_if<std::is_same<double, ValType>::value ||
                                  std::is_same<int, ValType>::value>::type * = nullptr>
void readGZ(gzFile &fp, std::vector<ValType> &vec) {
  int size1;
  gzread(fp, (void *)&size1, sizeof(int));
  vec.resize(size1);
  for (auto &v : vec) {
    gzread(fp, (void *)&v, sizeof(ValType));
  }
}

template <class ValType,
          typename std::enable_if<std::is_same<double, ValType>::value ||
                                  std::is_same<int, ValType>::value>::type * = nullptr>
void readGZ(gzFile &fp, std::vector<std::vector<ValType>> &vec) {
  int size1;
  gzread(fp, (void *)&size1, sizeof(int));
  vec.resize(size1);
  for (auto &line : vec) {
    int size2;
    gzread(fp, (void *)&size2, sizeof(int));
    line.resize(size2);
    for (auto &v : line) {
      gzread(fp, (void *)&v, sizeof(ValType));
    }
  }
}

int main(int argv, char **argc) {

  // MD data file
  gzFile gp = gzopen("mix.openmm.gz", "r");
  // OpenMM system:
  OpenMM::System system;
  int natoms = 72156;
  int ntatoms = 528;
  double total_mass = 0.0;
  std::vector<double> v_coul_ff_mass;
  std::vector<std::vector<double>> coord_vec;
  std::vector<OpenMM::Vec3> coord_init(natoms);
  system.setDefaultPeriodicBoxVectors(OpenMM::Vec3(18.970,0,0), OpenMM::Vec3(0,17.962,0), OpenMM::Vec3(0,0,11.629));
  readGZ(gp, v_coul_ff_mass);
  readGZ(gp, coord_vec);
  for (int i = 0; i < natoms; i++) {
     system.addParticle(v_coul_ff_mass[i]);
     total_mass += v_coul_ff_mass[i];
     coord_init[i] = OpenMM::Vec3(coord_vec[i][0], coord_vec[i][1], coord_vec[i][2]);
  }
#if (!defined XCG_OPENMM_SOFTCORE) || (XCG_OPENMM_SOFTCORE == SOFTCORE)
  // Add Coulumb
  OpenMM::NonbondedForce *coul_ff = new OpenMM::NonbondedForce();
  system.addForce(coul_ff);
#ifdef XCG_OPENMM_PME
  coul_ff->setNonbondedMethod(OpenMM::NonbondedForce::PME);
#else
  coul_ff->setNonbondedMethod(OpenMM::NonbondedForce::CutoffPeriodic);
#endif
  coul_ff->setCutoffDistance(cutoff_lj);
  std::vector<int> atIdx;
  std::vector<double> neigh14sceeList, neigh14scljList;
  std::vector<std::vector<int>> neigh14List;
  std::vector<std::vector<int>> neigh1213List;
  std::vector<std::vector<double>> v_coul_ff_atomic;
  std::vector<std::vector<double>> pairCoeff_14_eps, pairCoeff_14_sig;
  readGZ(gp, v_coul_ff_atomic);
  for (int i = 0; i < natoms; i++) {
#if (defined XCG_OPENMM_SOFTCORE) && (XCG_OPENMM_SOFTCORE == SOFTCORE)
     coul_ff->addParticle(0.0, v_coul_ff_atomic[i][1], v_coul_ff_atomic[i][2]);
#else
     coul_ff->addParticle(v_coul_ff_atomic[i][0], v_coul_ff_atomic[i][1], v_coul_ff_atomic[i][2]);
#endif
  }
#endif
#if (defined XCG_OPENMM_SOFTCORE) && (XCG_OPENMM_SOFTCORE == SOFTCORE)
  // Softcore Coulomb potential
  OpenMM::CustomNonbondedForce *coul_soft_ff = new OpenMM::CustomNonbondedForce("138.9306*q1*q2/sqrt(r*r+0.10)-offset;"
                                                                                "offset=138.9306*q1*q2/sqrt(rcut*rcut+0.10)");
  coul_soft_ff->addGlobalParameter("rcut", cutoff_lj);
  coul_soft_ff->addPerParticleParameter("q");
  coul_soft_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
  coul_soft_ff->setCutoffDistance(cutoff_lj);
  system.addForce(coul_soft_ff);
  for (int i = 0; i < natoms; i++) {
    coul_soft_ff->addParticle({v_coul_ff_atomic[i][0]});
  }
#endif

#if (!defined XCG_OPENMM_SOFTCORE) || (XCG_OPENMM_SOFTCORE == SOFTCORE)
  // Add 1-4 interactions
  readGZ(gp, pairCoeff_14_eps);
  readGZ(gp, pairCoeff_14_sig);
  readGZ(gp, atIdx);
  readGZ(gp, neigh14List);
  readGZ(gp, neigh14sceeList);
  readGZ(gp, neigh14scljList);
  readGZ(gp, neigh1213List);
  std::vector<std::vector<int>> exc_list = cat_vec_int(neigh14List, neigh1213List);
  for (int i = 0; i < exc_list.size(); i++) {
    int s1 = exc_list[i][0], s2 = exc_list[i][1];
    coul_ff->addException(s1, s2, 0.0, 0.0, 0.0);
#if (defined XCG_OPENMM_SOFTCORE) && (XCG_OPENMM_SOFTCORE == SOFTCORE)
    coul_soft_ff->addExclusion(s1, s2);
#endif
  }
  for(int i = 0; i < neigh14List.size(); i++) {
    int s1 = neigh14List[i][0], s2 = neigh14List[i][1];
    int s1t = atIdx[s1], s2t = atIdx[s2];
    double eps = pairCoeff_14_eps[s1t][s2t];
    double sig = pairCoeff_14_sig[s1t][s2t];
    coul_ff->addException(s1, s2, neigh14sceeList[i] * v_coul_ff_atomic[s1][0] * v_coul_ff_atomic[s2][0], sig, neigh14scljList[i] * eps, true);
  }
#endif

  // Add LJ potential

#if (!defined XCG_OPENMM_SOFTCORE) || (XCG_OPENMM_SOFTCORE == SOFTCORE)
  // LJ potential
#if (!defined XCG_OPENMM_SOFTCORE)
  OpenMM::CustomNonbondedForce *lj_ff = new OpenMM::CustomNonbondedForce("C12(type1, type2)/(r^12)-C6(type1, type2)/(r^6) - offset;"
                                                                         "offset=C12(type1, type2)/(rcut^12)-C6(type1, type2)/(rcut^6)");
#elif (XCG_OPENMM_SOFTCORE == SOFTCORE)
  OpenMM::CustomNonbondedForce *lj_ff = new OpenMM::CustomNonbondedForce("C12(type1, type2)/((ratio6)^2)-C12(type1, type2)/(ratio6) - offset;"
                                                                         "offset=C12(type1, type2)/((ratiocut6)^2)-C12(type1, type2)/(ratiocut6);"
                                                                         "ratio6=(r^6)/C6(type1, type2) + 0.10;"
                                                                         "ratiocut6=(rcut^6)/C6(type1, type2) + 0.10;");
#endif
  lj_ff->addGlobalParameter("rcut", cutoff_lj);
  lj_ff->addPerParticleParameter("type");
  lj_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
  lj_ff->setCutoffDistance(cutoff_lj);
  std::vector<double> C12values, C6values;
  readGZ(gp, C12values);
  readGZ(gp, C6values);
#if (defined XCG_OPENMM_SOFTCORE) && (XCG_OPENMM_SOFTCORE == SOFTCORE)
  for(int i = 0; i < C12values.size(); i++) {
    if (C12values[i] > 0.0 && C6values[i] == 0.0) {
       C6values[i] = sqrt(C12values[i] / 0.01);
       C12values[i] = 0.01;
    } else if (C12values[i] > 0.0 && C6values[i] > 0.0) {
       double sig6 = C12values[i] / C6values[i];
       C12values[i] = C12values[i] / sig6 / sig6;
       C6values[i] = sig6;
    } else if (C12values[i] == 0.0 && C6values[i] == 0.0) {
       C12values[i] = 0.0;
       C6values[i] = 1.0;
    }
  }
#endif
  OpenMM::Discrete2DFunction *tf_C12 = new OpenMM::Discrete2DFunction(528, 528, C12values);
  OpenMM::Discrete2DFunction *tf_C6 = new OpenMM::Discrete2DFunction(528, 528, C6values);
  lj_ff->addTabulatedFunction("C12", tf_C12);
  lj_ff->addTabulatedFunction("C6", tf_C6);
  system.addForce(lj_ff);
  for (int i = 0; i < natoms; i++) {
    lj_ff->addParticle({(double)atIdx[i]});
  }
  for (int i = 0; i < exc_list.size(); i++) {
    lj_ff->addExclusion(exc_list[i][0], exc_list[i][1]);
  }

#endif
  // Dipole-dipole & Dipole-charge potential
  double pref = 4.184*pow(1.82223, 4)/(8.314*temp/1000.0/4.184);
  std::vector<double> dipole_list;
  OpenMM::CustomNonbondedForce *dip_ff = new OpenMM::CustomNonbondedForce("dtmpi/(r^6)+dtmpj/(r^4) - offset;"
                                                                          "offset=dtmpi/(rcut^6)+dtmpj/(rcut^4);"
                                                                          "dtmpi=-2.0/3.0*pref*dsqi*dsqj; dtmpj=-0.33333*pref*dsqi*csqj-0.33333*pref*csqi*dsqj;"
                                                                          "csqi=C1*C1; csqj=C2*C2;"
                                                                          "dsqi=D1*D1; dsqj=D2*D2");
  dip_ff->addGlobalParameter("rcut", cutoff_lj);
  dip_ff->addGlobalParameter("pref", pref);
  dip_ff->addPerParticleParameter("C");
  dip_ff->addPerParticleParameter("D");
  dip_ff->setNonbondedMethod(OpenMM::CustomNonbondedForce::CutoffPeriodic);
  dip_ff->setCutoffDistance(cutoff_lj);
  //system.addForce(dip_ff);
  readGZ(gp, dipole_list);
  for (int i = 0; i < natoms; i++) {
    dip_ff->addParticle({v_coul_ff_atomic[i][0], dipole_list[i]});
  }
  for (int i = 0; i < exc_list.size(); i++) {
    dip_ff->addExclusion(exc_list[i][0], exc_list[i][1]);
  }

  // Bond potential (Harmonic) / Elastic Network Model
  OpenMM::HarmonicBondForce *bo_ff = new OpenMM::HarmonicBondForce();
  system.addForce(bo_ff);
  std::vector<std::vector<double>> bondList;
  readGZ(gp, bondList);
  for (auto &p : bondList) {
    if ((int)p[2] == 0) {
      if (v_coul_ff_mass[p[0] - 1] > 2.0 && v_coul_ff_mass[p[1] - 1] > 2.0) {
        bo_ff->addBond((int)p[0] - 1, (int)p[1] - 1, p[3], p[4]);
      } else {
#ifndef XCG_OPENMM_NOSHAKE
        system.addConstraint((int)p[0] - 1, (int)p[1] - 1, p[3]);
#else
        bo_ff->addBond((int)p[0] - 1, (int)p[1] - 1, p[3], p[4]);
#endif
      }
    }
  }

  // Angle potential (Harmonic)
  OpenMM::HarmonicAngleForce *an_harm_ff = new OpenMM::HarmonicAngleForce();
  system.addForce(an_harm_ff);
  std::vector<std::vector<double>> angleList;
  readGZ(gp, angleList);
  for (auto &p : angleList) {
    if ((int)p[3] == 0) {
      an_harm_ff->addAngle((int)p[0] - 1, (int)p[1] - 1, (int)p[2] - 1, p[4], p[5]);
    }
  }

  // Dihedral Potential (Fourier)
  OpenMM::PeriodicTorsionForce *di_prop = new OpenMM::PeriodicTorsionForce();
  system.addForce(di_prop);
  std::vector<std::vector<double>> dihedralList;
  readGZ(gp, dihedralList);
  for (auto &p : dihedralList) {
    if ((int)p[4] == 0 && (int)p[5] > 0) {
      di_prop->addTorsion((int)p[0] - 1, (int)p[1] - 1, (int)p[2] - 1, (int)p[3] - 1, p[5], p[6], p[7]);
    }
  }

  // Improper Dihedral Potential (Fourier)
  OpenMM::PeriodicTorsionForce *di_fo_impr = new OpenMM::PeriodicTorsionForce();
  system.addForce(di_fo_impr);
  std::vector<std::vector<double>> improperList;
  readGZ(gp, improperList);
  for (auto &p : improperList) {
    if ((int)p[4] == 0 && (int)p[5] > 0) {
      di_fo_impr->addTorsion((int)p[0] - 1, (int)p[1] - 1, (int)p[2] - 1, (int)p[3] - 1, p[5], p[6], p[7]);
    }
  }

#ifdef XCG_OPENMM_PLUMED
  auto script = dump_file("plumed.dat");
  PlumedPlugin::PlumedForce *rst_plumed = new PlumedPlugin::PlumedForce(script);
  system.addForce(rst_plumed);
#endif

#ifdef XCG_OPENMM_NPT
  OpenMM::MonteCarloBarostat *p_stat = new OpenMM::MonteCarloBarostat(pres, temp);
  system.addForce(p_stat);
#endif

  OpenMM::LangevinIntegrator *langevin = new OpenMM::LangevinIntegrator(temp, 2.0, timestep);

  // Path to your plugin directory
  OpenMM::Platform::loadPluginsFromDirectory(XCG_OPENMM_PLUGPATH);
  OpenMM::Context *context = new OpenMM::Context(system, *langevin);
  context->setPositions(coord_init);

  // L-BFGS optimizer, occasionally crashes due to the improper initial confomation
  OpenMM::LocalEnergyMinimizer *min = new OpenMM::LocalEnergyMinimizer();
  min->minimize(*context, 10, 1000);

  // Overall # of traj snapshots
  int noutputs = ntotalSteps / output_freq;

  // Output dcd
  FILE *fp_dcd = fopen("mix.openmm.dcd", "wb");

  // Output box vectors
  OpenMM::Vec3 a, b, c;

  write_header(fp_dcd, natoms, "XCG generated dcd", output_freq);
  auto t_start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < noutputs; i++) {
    OpenMM::State state = context->getState((OpenMM::State::Positions | OpenMM::State::Energy), true);
    const std::vector<OpenMM::Vec3>& positionsInNm = state.getPositions();
    state.getPeriodicBoxVectors(a, b, c);
    if (output_freq <= 1000 && i % (10000/output_freq) == 0) {
      const double pot_eng = state.getPotentialEnergy();
      double density = total_mass * (1.0e21) / (a[0] * b[1] * c[2] * 6.022e23);
      auto t_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t_start);
      if (std::isnan(pot_eng) || std::isinf(pot_eng)) {
        fclose(fp_dcd);
        return 1;
      }
      printf("progress=%.3f(%%) step=%d pot_eng=%.6e(kJ/mol) density=%.3e(g/cm^3) time_elapsed=%.3f(s)\n", ((double)i)/noutputs*100, i * output_freq, pot_eng, density, double(t_elapsed.count())/1000.0);
    }
    // writeXYZ("mix.openmm.traj.xyz", positionsInNm);
    write_box(fp_dcd, a[0], b[1], c[2]);
    write_frame(fp_dcd, natoms, positionsInNm, i, output_freq);
    langevin->step(output_freq);
  }
  auto t_end = std::chrono::high_resolution_clock::now();
  auto t_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);
  printf("Total time elapsed: %.3f(s)\n", double(t_elapsed.count())/1000.0);
  std::system("touch Stamp_OpenmmNormalEnd_mix");

  fclose(fp_dcd);
}
