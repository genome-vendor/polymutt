#ifndef __CMDLINEPAR_H__
#define __CMDLINEPAR_H__

#include "StringArray.h"
#include  <string>

class CmdLinePar
{
public:
std::string cmd;
int minTotalDepth;
int maxTotalDepth;
double minAvgDepth;
double maxAvgDepth;
double minPS;
double minMapQuality;
double posterior;
double theta;
double theta_indel;
double precision;
bool  denovo;
double denovo_mut_rate;
double denovo_tstv_ratio;
double denovoLR;
bool output_denovo_only;
bool gl_off;
String chrX_label;
String chrY_label;
String MT_label;
String vcfInFile;
String vcfOutFile;
bool force_call;
bool out_all_sites;

public:
CmdLinePar() {
  theta = .0;
  minTotalDepth = 0;
  maxTotalDepth = 0;
  minPS = 0;
  minAvgDepth = 0.;
  maxAvgDepth = 0.;
  minMapQuality = 0.;
  posterior = 0.5;
  theta = 0.001;
  theta_indel = 0.001;
  precision = 0.0001;
  denovo = false;
  denovo_mut_rate = 0.0;
  denovo_tstv_ratio = 0.5;
  denovoLR = 0;
  output_denovo_only = false;
  gl_off = false;
  String chrX_label = "";
  String chrY_label = "";
  String MT_label = "";
  force_call = false;
  out_all_sites = false;
 }
};
        
#endif
