#ifndef __FAMILYLIKELIHOODSEQ_VCF_H__
#define __FAMILYLIKELIHOODSEQ_VCF_H__

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

class FamilyLikelihoodSeq_VCF : public FamilyLikelihoodSeq
{
 public:
  Pedigree *ped;
  std::vector<std::string> pids;
  std::vector<std::string> includedPIDs;
  std::map<std::string, std::pair<int, int> > pid2traverse;
  std::map<std::string, int> included;
  std::string vcfInFile;
  int GL_idx;
  int PL_idx;
  int ref;
  int alt;
  int position;
  std::string refStr;
  std::string altStr;
  VCFInputFile *vin;
  VCFRecord *r;
  VCFPeople *people;
  int nFields;
  double QUAL;
  int GT_index, GQ_index, DP_index;
  std::vector<std::vector<int> > sexes; 
  int maleFounders, femaleFounders;
  double PL2LK_table[1024];
  bool isBiallelic;
  bool pid_match;
  int withdata_cnt;
  int nSamplesWithData;

 public:
  FamilyLikelihoodSeq_VCF();
  ~FamilyLikelihoodSeq_VCF();
  double f(double);
  void SetPedigree(Pedigree *p);
  void SetVinPtr(VCFInputFile *v) { vin = v; }
  void MapPID2Traverse();
  double PL2LK(int);
  void FillPenetrance();
  int Allele2Int(std::string&);
  void FillZeroPenetrance(int fam_idx, int person_idx, int geno_idx);
  void GetSexes();
  void PrintPenetranceMatrix();
  void PrintLogLKMatrix();
  double MonomorphismLogLikelihood(int);
  double PolymorphismLogLikelihood(int  a1, int a2);
  double CalcAllFamLogLikelihood(double);
  double CalcSingleFamLogLikelihood_Founders(int i, double freq);
  void CalcPostProb(double freq);
  void CalcPostProb_SingleExtendedPed_BA_VCF(int i, double freq);
  void CalcPostProb_SingleFam_BA_Founders_VCF(int, double);
  void CalcPostProb_SinglePerson(int, int, double);
  double lkSinglePerson(int i, int j, double freq);
  void SetQUAL(double qual) { QUAL = qual; }
  void OutputVCF(FILE*);
  void OutputVCF2(FILE*);
  double lkSingleFam(int i, double freq);
  double lkSingleFam_denovo(int i, double freq);
  void CalcParentMarginal(int i, double freq);
  void CalcParentMarginal_denovo(int i, double freq);
  void CalcPostProb_SingleNucFam(int i, double freq);
  JointGenoLk KidJointGenoLikelihood(int famIdx, int kidIdx);
  JointGenoLk_denovo KidJointGenoLikelihood_denovo(int famIdx, int kidIdx);
  double likelihoodKids(int fa1, int fa2, int ma1, int ma2, int famIdx);
  double likelihoodONEKid(int famIdx, int personIdx, int fa1, int fa2, int ma1, int ma2);
  JointGenoLk likelihoodKidGenotype(int fa1, int fa2, int ma1, int ma2, int famIdx, int kidIndex);
};

#endif
