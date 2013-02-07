#ifndef __PEDVCF_H__
#define __PEDVCF_H__

#include "NucFamGenotypeLikelihood.h"
#include "StringMap.h"
#include <stdlib.h>
#include <iostream>
#include "CmdLinePar.h"
#include "MutationModel.h"
#include <map>
#include "FamilyLikelihoodES.h"
#include "FamilyLikelihoodSeq.h"

#include <cassert>
#include <string>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include "Argument.h"
#include "IO.h"
#include "tabix.h"
#include "PeopleSet.h"
#include "RangeList.h"
#include "VCFUtil.h"

class PedVCF
{
 public:
  std::string vcfInFile;
  CmdLinePar par;
  Pedigree *ped;
  PedigreeGLF *pedGLF;
  double tstv_ratio;
  std::string meta_data;
 public:
  PedVCF();
  ~PedVCF();
  double GetPrior_ts(double tstv_ratio);
  double GetPrior_tv(double tstv_ratio);
  bool isTs(int a1, int a2);
  bool isTv(int a1, int a2);
  void SetPedigreePtr(Pedigree *p){ ped = p ;}
  void SetPedigreeGLF(PedigreeGLF *p) { pedGLF = p; }
  void SetPar(CmdLinePar&);
  void SetVCFInput(String& vcfInput);
  void VarCallFromVCF();
};

#endif
