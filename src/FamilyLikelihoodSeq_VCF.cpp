#include "FamilyLikelihoodSeq_VCF.h"

#define GLLIM 255

FamilyLikelihoodSeq_VCF::FamilyLikelihoodSeq_VCF()
{
 ped = NULL;
 GL_idx = -1;
 PL_idx = -1;
 GT_index = -1;
 GQ_index = -1;
 DP_index = -1;
 QUAL = 0;
 isIndel = false;
 pid_match = false;
 isChrX = false;
 isChrY = false;
 isMT = false;
 withdata_cnt = 0;
 nSamplesWithData = 0;
 for(int i=0; i<256; i++)
  PL2LK_table[i] = pow(10, -double(i)/10.0);
}

FamilyLikelihoodSeq_VCF::~FamilyLikelihoodSeq_VCF(){}


void FamilyLikelihoodSeq_VCF::SetPedigree(Pedigree *p)
{
  ped = p;
}

double FamilyLikelihoodSeq_VCF::f(double freq)
{
  return(-(CalcAllFamLogLikelihood(freq)));
}

void FamilyLikelihoodSeq_VCF::MapPID2Traverse()
{
    pid2traverse.clear();
    pids.clear();

    for(int i=0; i<ped->familyCount; i++)
        for(int j=0; j<ped->families[i]->count; j++)
        {
            std::pair<int, int> traverse_pair;
            traverse_pair.first = i;
            traverse_pair.second = j;

            Person *p;
            p = ped->persons[ped->families[i]->path[j]];
            pid2traverse[std::string(p->pid.c_str())] = traverse_pair;
            pids.push_back(std::string(p->pid.c_str()) );
        }
}

double FamilyLikelihoodSeq_VCF::PL2LK(int pl)
{
  if(pl<0) error("Phred-scaled likelihood %d can not be negative\n", pl);
  if(pl>GLLIM) return(PL2LK_table[GLLIM]);
  //if(pl<50) return(pow(10, -0.1 * pl)); else{ printf("PL255: %g\n", PL2LK_table[255]); return(PL2LK_table[255]); }
  return(PL2LK_table[pl]);
}

int FamilyLikelihoodSeq_VCF::Allele2Int(std::string& allele)
{
 if(!allele.compare("A") || !allele.compare("a") ) return (1);
 if(!allele.compare("C") || !allele.compare("c") ) return (2);
 if(!allele.compare("G") || !allele.compare("g") ) return (3);
 if(!allele.compare("T") || !allele.compare("t") ) return (4);
 return(0);
}

double FamilyLikelihoodSeq_VCF::MonomorphismLogLikelihood(int refBase)
{
  double llkRef = 0.0;
  
  int homoRefIdx = glfHandler::GenotypeIndex(refBase, refBase);
  for(int i=0; i<nFam; i++)
    for(int j=0; j<ped->families[i]->count; j++)
        llkRef += fam[i].loglk[j][homoRefIdx]; 
  return(llkRef);
}

double FamilyLikelihoodSeq_VCF::PolymorphismLogLikelihood(int  a1, int a2)
{
  SetAlleles(a1, a2);
  OptimizeFrequency();
  return(GetMaxLogLikelihood());
}

double FamilyLikelihoodSeq_VCF::CalcAllFamLogLikelihood(double freq)
{
  double loglk = 0.0;

#pragma omp parallel for reduction(+:loglk)
  for(int i=0; i<nFam; i++)
  {
   if(ped->families[i]->count==ped->families[i]->founders)
	loglk += CalcSingleFamLogLikelihood_Founders(i, freq);
   else if(ped->families[i]->isNuclear() && nFam>1 && !isChrX && !isChrY && !isMT)
	loglk += log10(lkSingleFam(i, freq) );
   else loglk +=  CalcSingleFamLogLikelihood_BA(i, freq);

   //printf("freq=%f a1=%d a2=%d famidx=%d loglk=%f %f\n", freq, allele1, allele2, i, log10(lkSingleFam(i, freq) ), CalcSingleFamLogLikelihood_BA(i, freq) );
  }

  return(loglk);
}

double FamilyLikelihoodSeq_VCF::CalcSingleFamLogLikelihood_Founders(int i, double freq)
{
 double llk = 0.0;

 for(int j=0; j<ped->families[i]->count; j++)
    llk += log10(lkSinglePerson(i, j, freq) );

 return(llk);
}

double FamilyLikelihoodSeq_VCF::lkSinglePerson(int i, int j, double freq)
{
  double sum=0.0;
  double lk11, lk12, lk22;
  lk11 = fam[i].penetrances[j][geno11];
  lk12 = fam[i].penetrances[j][geno12];
  lk22 = fam[i].penetrances[j][geno22];

  double priors[3];
  priors[0] = freq*freq;
  priors[1] = freq*(1-freq)*2;
  priors[2] = (1-freq)*(1-freq);

  if(isChrX){if(pedGLF->sexes[i][j]==MALE) { lk12=0; priors[0]=freq; priors[1]=0; priors[2]=1-freq;} }
  if(isChrY){if(pedGLF->sexes[i][j]==MALE) { lk12=0; priors[0]=freq; priors[1]=0; priors[2]=1-freq;} else return(1.0); } 
  if(isMT){ lk12=0; priors[0] = freq; priors[1]=0; priors[2]=1-freq;}

  sum = sum + lk11*priors[0] + lk12*priors[1] + lk22*priors[2];
  return(sum);
}

void FamilyLikelihoodSeq_VCF::CalcPostProb(double freq)
{
  for(int i=0; i<nFam; i++)
  {
    if(ped->families[i]->count==ped->families[i]->founders)
	CalcPostProb_SingleFam_BA_Founders_VCF(i, freq);
    else if(ped->families[i]->isNuclear() && nFam>1 && !isChrX && !isChrY && !isMT)
	CalcPostProb_SingleNucFam(i, freq);
    else CalcPostProb_SingleExtendedPed_BA_VCF(i, freq);
 }
  CalcGQ();
}

void FamilyLikelihoodSeq_VCF::CalcPostProb_SingleFam_BA_Founders_VCF(int i, double freq)
{
   for(int j=0; j<pedGLF->ped->families[i]->founders; j++)
   {
    sex = pedGLF->sexes[i][j];
    CalcPostProb_SinglePerson(i, j, freq);
   }
}

void FamilyLikelihoodSeq_VCF::CalcPostProb_SinglePerson(int i, int j, double freq)
{
  double lk11, lk12, lk22;
  double priors[3];
  priors[0] = freq*freq;
  priors[1] = freq*(1-freq)*2;
  priors[2] = (1-freq)*(1-freq);

  lk11 = fam[i].penetrances[j][geno11];
  lk12 = fam[i].penetrances[j][geno12];
  lk22 = fam[i].penetrances[j][geno22];

  if(isChrX)
  { if(pedGLF->sexes[i][j]==MALE){ priors[0] = freq; priors[1]=0.; priors[2]=1-freq; }
    else { priors[0] = freq*freq; priors[1] = 2*freq*(1-freq); priors[2] = (1-freq)*(1-freq);}
  }
  if(isChrY)
  { if(pedGLF->sexes[i][j]==MALE) { priors[0] = freq; priors[1] = 0.; priors[2] = 1-freq; }
    else { priors[0]=priors[1]=priors[2]=1.0; }
  }
  if(isMT)
  { priors[0] = freq; priors[1] = 0; priors[2] = 1-freq; 
  }

  double mlk11 = lk11*priors[0];
  double mlk12 = lk12*priors[1];
  double mlk22 = lk22*priors[2];
  double sum = mlk11+mlk12+mlk22;

  if(sum==0) postProb[i][j][0]=postProb[i][j][1]=postProb[i][j][2]=1/3;
  else {
    postProb[i][j][0] =  mlk11/sum;
    postProb[i][j][1] =  mlk12/sum;
    postProb[i][j][2] =  mlk22/sum;
  }

  if(isChrY && pedGLF->sexes[i][j]==FEMALE) { postProb[i][j][0] = postProb[i][j][1] = postProb[i][j][2] = 0.0; }  

  int best = GetBestGenoIdx(mlk11, mlk12, mlk22);
  bestGenoIdx[i][j] =  best;
  bestGenoLabel[i][j] = isChrY && pedGLF->sexes[i][j]==FEMALE ? "." : GetBestGenoLabel_vcfv4(best);
  //bestGenoLabel[i][j] = GetBestGenoLabel(best);
  dosage[i][j] = CalcDosage(i, j);
}


void FamilyLikelihoodSeq_VCF::CalcPostProb_SingleExtendedPed_BA_VCF(int i, double freq)
{
  double lk11, lk12, lk22;
  double sum;
  int best;

  Matrix tmp;
  for(int j=0; j<ped->families[i]->count; j++)
    {
      sex = pedGLF->sexes[i][j];

      if(isChrY && pedGLF->sexes[i][j]==FEMALE)
      {
        bestGenoIdx[i][j] = 0;
        bestGenoLabel[i][j] = ".";
        postProb[i][j][0]=postProb[i][j][1]=postProb[i][j][2]=0.0;
        dosage[i][j] = 0;
        continue;
      }



      tmp = fam[i].penetrances; //copy is not nessary and consider optimization

      FillZeroPenetrance(i, j, geno11);
      lk11 = CalcSingleFamLikelihood_BA(i, freq);

      fam[i].penetrances = tmp;

      FillZeroPenetrance(i, j, geno12);
      lk12 = CalcSingleFamLikelihood_BA(i, freq);

      fam[i].penetrances = tmp;

      FillZeroPenetrance(i, j, geno22);
      lk22 = CalcSingleFamLikelihood_BA(i, freq);

      fam[i].penetrances = tmp;

      sum = lk11+lk12+lk22;
      if(sum==0)
        postProb[i][j][0] = postProb[i][j][1] = postProb[i][j][2] = 0.0;
      else {
        postProb[i][j][0] = lk11/sum;
        postProb[i][j][1] = lk12/sum;
        postProb[i][j][2] = lk22/sum;
      }

      best = GetBestGenoIdx(lk11, lk12, lk22);
      bestGenoIdx[i][j] = best; 
      //bestGenoLabel[i][j] = GetBestGenoLabel(best);
      //bestGenoLabel[i][j] = isChrY && pedGLF->sexes[i][j]==FEMALE ? "./." : GetBestGenoLabel_vcfv4(best);
      bestGenoLabel[i][j] = GetBestGenoLabel_vcfv4(best); //printf("bg: %s %f\n", bestGenoLabel[i][j].c_str(), QUAL);
      dosage[i][j] = CalcDosage(i,j);
    }
}

void FamilyLikelihoodSeq_VCF::FillPenetrance()
{
 r = &vin->getVCFRecord();
 people = &r->getPeople();

 refStr = r->getRef();
 altStr = r->getAlt();
 position = r->getPos();

 for(int i=0; i<nFam; i++)
 {
  fam[i].loglk.Set(0.0);
  fam[i].penetrances.Set(1.0);
 }

  std::map<std::string, std::pair<int, int> >::iterator it;
  StringArray gl;
  VCFIndividual* indv;

  if(nSamplesWithData==0)
  for(int i=0; i<people->size(); i++)
  {
      indv = (*people)[i];
//      it = pid2traverse.find( (*indv).getName() ); //Note: this should be done only once on the header line. Will optimaize later
//      if(it==pid2traverse.end() ) continue; //Note: after the above, this should be checked a vector of bools
	if(pid2traverse.count((*indv).getName())==0) {printf("Sample ID \"%s\" not included in the analysis!\n", (*indv).getName().c_str());continue; }
	nSamplesWithData++;
  }

 if(refStr.compare(altStr)==0) { isBiallelic = false; return; } // Ignore monomorphic sites.

 isBiallelic = true;
 for(int s=0; s<altStr.size();s++)
  if(altStr[s]==',')
    { isBiallelic = false; return; }

 isIndel =  (refStr.size()>1 || altStr.size()>1) ? true : false;

 ref = isIndel ? 1 : Allele2Int(refStr);
 alt = isIndel ? 2 : Allele2Int(altStr);
 allele1 = ref;
 allele2 = alt;
 pedGLF->refBase = ref;

  if(DP_index<0)
  {
	DP_index = r->getFormatIndex("DP");
  }

  if(GL_idx<0 && PL_idx<0)
  {
      GL_idx = r->getFormatIndex("GL");
      PL_idx = r->getFormatIndex("PL");

      if(GL_idx<0 && PL_idx<0) {
         fprintf(stderr, "NO GL or PL field was found. Please check the vcf file at chr:%s and position:%d", r->getChrom(), r->getPos());
          exit(1);
       }

     for(int i=0; i<people->size(); i++)
      {
	indv = (*people)[i];
        it = pid2traverse.find( (*indv).getName() );
	if(it!=pid2traverse.end()) { pid_match = true; break;}
      }

    if(pid_match==false) error("NO individual IDs match in the ped and vcf file!\n");
   }

   if(isBiallelic==false) return;

   withdata_cnt = 0;

    for(int i=0; i<people->size(); i++)
    {
        indv = (*people)[i];
        it = pid2traverse.find( (*indv).getName() ); //Note: this should be done only once on the header line. Will optimaize later
        if(it==pid2traverse.end() ) continue; //Note: after the above, this should be checked a vector of bools

	bool missing = false;

	int gt_idx0 = glfHandler::GenotypeIndex(ref, ref);
        int gt_idx1 = glfHandler::GenotypeIndex(ref, alt);
        int gt_idx2 = glfHandler::GenotypeIndex(alt, alt);

	gl.ReplaceTokens(String( GL_idx>0 ? (*indv).get(GL_idx, &missing).toStr() : (*indv).get(PL_idx, &missing).toStr() ), ",");

        if(missing) {
          fam[it->second.first].loglk[it->second.second][gt_idx0] = 0.0;
          fam[it->second.first].loglk[it->second.second][gt_idx1] = 0.0;
          fam[it->second.first].loglk[it->second.second][gt_idx2] = 0.0;

          fam[it->second.first].penetrances[it->second.second][gt_idx0] = 1.0;
          fam[it->second.first].penetrances[it->second.second][gt_idx1] = 1.0;
          fam[it->second.first].penetrances[it->second.second][gt_idx2] = 1.0; 

	  return;
         }

	if(gl.Length()!=3) error("GL or PL filed does not have 3 values separated by commas at: %s %d!\n", r->getChrom(), r->getPos());

        double GL0 = gl[0].AsDouble();
        double GL1 = gl[1].AsDouble();
        double GL2 = gl[2].AsDouble();

	if(GL0!=0.0 || GL1!=0.0 || GL2!=0.0) withdata_cnt++;

        fam[it->second.first].loglk[it->second.second][gt_idx0] = PL_idx>0 ? (GL0 > GLLIM ? -GLLIM/10.0 : -GL0/10.0) : (-10*GL0>GLLIM ? -GLLIM/10.0 : GL0);
        fam[it->second.first].loglk[it->second.second][gt_idx1] = PL_idx>0 ? (GL1 > GLLIM ? -GLLIM/10.0 : -GL1/10.0) : (-10*GL1>GLLIM ? -GLLIM/10.0 : GL1);
        fam[it->second.first].loglk[it->second.second][gt_idx2] = PL_idx>0 ? (GL2 > GLLIM ? -GLLIM/10.0 : -GL2/10.0) : (-10*GL2>GLLIM ? -GLLIM/10.0 : GL2);

        fam[it->second.first].penetrances[it->second.second][gt_idx0] = PL2LK(int(PL_idx>0 ? GL0 : -10*GL0));
        fam[it->second.first].penetrances[it->second.second][gt_idx1] = PL2LK(int(PL_idx>0 ? GL1 : -10*GL1));
        fam[it->second.first].penetrances[it->second.second][gt_idx2] = PL2LK(int(PL_idx>0 ? GL2 : -10*GL2));

    } //end of iterating all individuals in a VCF record
}

void FamilyLikelihoodSeq_VCF::FillZeroPenetrance(int fam_idx, int person_idx, int geno_idx)
{
	if(geno_idx==geno11) { fam[fam_idx].penetrances[person_idx][geno12] = fam[fam_idx].penetrances[person_idx][geno22] =  0.0; }
	else if(geno_idx==geno12) { fam[fam_idx].penetrances[person_idx][geno11] = fam[fam_idx].penetrances[person_idx][geno22] =  0.0; }
	else if(geno_idx==geno22) { fam[fam_idx].penetrances[person_idx][geno11] = fam[fam_idx].penetrances[person_idx][geno12] =  0.0; }
  	else error("The genotpe index %d is not valide\n");
}

void FamilyLikelihoodSeq_VCF::GetSexes()
{
 sexes.resize(ped->familyCount);
 for(int i=0; i<ped->familyCount; i++)
 {
  sexes[i].resize(ped->families[i]->count);
  for(int j=0; j<ped->families[i]->count; j++)
  {
   Person *p = ped->persons[ped->families[i]->path[j]];
   sexes[i][j] = p->sex;
   if(p->sex==MALE && p->isFounder()) maleFounders++;
   if(p->sex==FEMALE && p->isFounder()) femaleFounders++;
  }
 }
}

void FamilyLikelihoodSeq_VCF::PrintPenetranceMatrix()
{
 for(int i=0; i<ped->familyCount; i++)
  for(int j=0; j<ped->families[i]->count; j++)
  {
    printf("%d\t%d:", i, j);
    for(int k=0; k<10; k++)
      printf("\t%f", (fam[i].penetrances[j][k]));
    printf("\n");
  }
}


void FamilyLikelihoodSeq_VCF::PrintLogLKMatrix()
{
 for(int i=0; i<ped->familyCount; i++)
  for(int j=0; j<ped->families[i]->count; j++)
  {
    printf("%d\t%d:", i, j);
    for(int k=0; k<10; k++)
      printf("\t%f", (fam[i].loglk[j][k]));
    printf("\n");
  }
}

// This will update the individuals in the VCF file that are also in the .ped file
// Those in VCF but not in .ped file are kept untouched
// The order of individuals is the same as in the original VCF file
void FamilyLikelihoodSeq_VCF::OutputVCF(FILE *ofh)
{

int bestIdx;
int GTQual;

double polyPrior = GetPolyPrior();
std::map<std::string, std::pair<int, int> >::iterator it;

if(isBiallelic==false) return;

bool missing = false;
//double dp = r->getInfoTag("DP", &missing).toDouble(); printf("DP=%f\n", dp);
//std::string dpStr = r->getInfoTag("DP", &missing).toStr();
//std::string dpInfo("");
//if(!missing) { dpInfo = std::string(";DP=") + dpStr;}
//else dp = 0.0;

VCFIndividual* indv;
AC = 0;
int totalSamplesWithDepth = 0;
int totalDepth = 0;
std::string pid;

for(int i=0; i<people->size(); i++) 
{
        indv = (*people)[i];
	pid = (*indv).getName();

	if(included[pid]>0==false) continue; 
	std:pair<int, int> fam_pid = pid2traverse[pid];
	AC += bestGenoIdx[fam_pid.first][fam_pid.second]; 

	int dp = DP_index>0 ? (*indv).get(DP_index, &missing).toInt() : 0;
	if(missing) continue;

	totalDepth+=dp;
	totalSamplesWithDepth++;
}

char buff[100];
sprintf(buff,"%d", totalDepth);
std::string dpInfo(";DP=");
dpInfo += std::string(buff);

fprintf(ofh, "%s\t%d\t%s\t%s\t%s\t%.2f\t%s\tAF=%.2f;AC=%d%s\t%s",
			r->getChrom(),
                         r->getPos(),
                         r->getID(),
                         //r->getQual(),
                         r->getRef(),
                         r->getAlt(),
			QUAL,
                         //r->getQual(),
                        // totalDepth<=1 || double(totalDepth)/double(totalSamplesWithDepth)<=1 ? "LOWDP" : ".", 
			r->getFilt(),
                         1-GetMinimizer(),
			 AC,
			 dpInfo.c_str(),
                         PL_idx>0 ? "GT:GQ:DP:PL" : "GT:GQ:DP:GL"//r->getFormat()
	);

for(int i=0; i<people->size(); i++)
   {
        indv = (*people)[i];
	pid = (*indv).getName();
        if(included[pid]>0==false) continue;

	std::pair<int, int> fam_pid;
	fam_pid = pid2traverse[pid];
	//fprintf(ofh, "%d:", bestGenoIdx[it->second.first][it->second.second] );
	fprintf(ofh, "\t");
	fprintf(ofh, "%s:", bestGenoQual[fam_pid.first][fam_pid.second]>0 || bestGenoLabel[fam_pid.first][fam_pid.second]=="." ? bestGenoLabel[fam_pid.first][fam_pid.second].c_str() : "./." );
	fprintf(ofh, "%d:", bestGenoQual[fam_pid.first][fam_pid.second] );

	const char * indvDP = DP_index>0 ? (*indv).get(DP_index, &missing).toStr() : ".";
	fprintf(ofh, "%s:", missing ? "." : indvDP);

	const char *PL_value = (*indv).get(PL_idx>0 ? PL_idx : GL_idx, &missing).toStr();
	fprintf(ofh, "%s",  missing ? "." : PL_value );
    }

   fprintf(ofh, "\n");
   fflush(ofh);
}

//################################################################
//#### The following functions are for nuclear families
//#################################################################
double FamilyLikelihoodSeq_VCF::lkSingleFam(int i, double freq)
{
  double sum = 0.0;

  CalcParentMarginal(i, freq);
  for(int idx=0; idx<9; idx++)
    sum += parentMarginal[i][idx]; //Note that Marginals here are joint likelihood of parents marginalizing offspring

  return(sum);
}

double FamilyLikelihoodSeq_VCF::lkSingleFam_denovo(int i, double freq)
{
  if(pedGLF->ped->families[i]->count==pedGLF->ped->families[i]->founders)
    {
   double lk = 1.0;
        for(int j=0; j<pedGLF->ped->families[i]->founders; j++)
     lk *= lkSinglePerson(i, j, freq);
   return lk;
  }
  double sum = 0.0;

  CalcParentMarginal_denovo(i, freq);
  for(int idx=0; idx<9; idx++)
    sum += parentMarginal[i][idx]; //Note that Marginals here are joint likelihood of parents marginalizing offspring

  return(sum);
}


void FamilyLikelihoodSeq_VCF::CalcParentMarginal(int i, double freq)
{
  double lkF11, lkF12, lkF22, lkM11, lkM12, lkM22; // lkC11, lkC12, lkC22;
  lkF11=lkF12=lkF22=lkM11=lkM12=lkM22=1.0;

  lkF11 = fam[i].penetrances[0][geno11];
  lkF12 = fam[i].penetrances[0][geno12];
  lkF22 = fam[i].penetrances[0][geno22];
  lkM11 = fam[i].penetrances[1][geno11];
  lkM12 = fam[i].penetrances[1][geno12];
  lkM22 = fam[i].penetrances[1][geno22];

 if(isChrX) lkF12 = 0.0;
 if(isChrY) { lkM11=lkM12=lkM22=1.0; lkF12=0.0; }
 if(isMT) lkF12=lkM12=0.0;

  parentGLF[i][0] = lkF11*lkM11;
  parentGLF[i][1] = lkF11*lkM12;
  parentGLF[i][2] = lkF11*lkM22;
  parentGLF[i][3] = lkF12*lkM11;
  parentGLF[i][4] = lkF12*lkM12;
  parentGLF[i][5] = lkF12*lkM22;
  parentGLF[i][6] = lkF22*lkM11;
  parentGLF[i][7] = lkF22*lkM12;
  parentGLF[i][8] = lkF22*lkM22;


  if(nFam>1)
    SetParentPrior(freq);
  else 
    SetParentPriorSingleTrio();

  parentConditional[i][0] = likelihoodKids(allele1,allele1,allele1,allele1, i)*parentGLF[i][0];
  parentConditional[i][1] = likelihoodKids(allele1,allele1,allele1,allele2, i)*parentGLF[i][1];
  parentConditional[i][2] = likelihoodKids(allele1,allele1,allele2,allele2, i)*parentGLF[i][2];
  parentConditional[i][3] = likelihoodKids(allele1,allele2,allele1,allele1, i)*parentGLF[i][3];
  parentConditional[i][4] = likelihoodKids(allele1,allele2,allele1,allele2, i)*parentGLF[i][4];
  parentConditional[i][5] = likelihoodKids(allele1,allele2,allele2,allele2, i)*parentGLF[i][5];
  parentConditional[i][6] = likelihoodKids(allele2,allele2,allele1,allele1, i)*parentGLF[i][6];
  parentConditional[i][7] = likelihoodKids(allele2,allele2,allele1,allele2, i)*parentGLF[i][7];
  parentConditional[i][8] = likelihoodKids(allele2,allele2,allele2,allele2, i)*parentGLF[i][8];

  //Mariginal is joint likelihood of parents marginalizing all kids
  //This is to avoid redundant calculation
  for(int j=0; j<9; j++)
   parentMarginal[i][j] = parentConditional[i][j]*parentPrior[j];

}

void FamilyLikelihoodSeq_VCF::CalcParentMarginal_denovo(int i, double freq)
{
  double lkF11, lkF12, lkF22, lkM11, lkM12, lkM22; // lkC11, lkC12, lkC22;
  lkF11=lkF12=lkF22=lkM11=lkM12=lkM22=1.0;

  lkF11 = fam[i].penetrances[0][geno11];
  lkF12 = fam[i].penetrances[0][geno12];
  lkF22 = fam[i].penetrances[0][geno22];
  lkM11 = fam[i].penetrances[1][geno11];
  lkM12 = fam[i].penetrances[1][geno12];
  lkM22 = fam[i].penetrances[1][geno22];

  parentGLF[i][0] = lkF11*lkM11;
  parentGLF[i][1] = lkF11*lkM12;
  parentGLF[i][2] = lkF11*lkM22;
  parentGLF[i][3] = lkF12*lkM11;
  parentGLF[i][4] = lkF12*lkM12;
  parentGLF[i][5] = lkF12*lkM22;
  parentGLF[i][6] = lkF22*lkM11;
  parentGLF[i][7] = lkF22*lkM12;
  parentGLF[i][8] = lkF22*lkM22;


  if(nFam>1)
    SetParentPrior_denovo(freq);
  else 
    SetParentPriorSingleTrio_denovo(freq);

  parentConditional[i][0] = likelihoodKids_denovo(allele1,allele1,allele1,allele1, i)*parentGLF[i][0];
  parentConditional[i][1] = likelihoodKids_denovo(allele1,allele1,allele1,allele2, i)*parentGLF[i][1];
  parentConditional[i][2] = likelihoodKids_denovo(allele1,allele1,allele2,allele2, i)*parentGLF[i][2];
  parentConditional[i][3] = likelihoodKids_denovo(allele1,allele2,allele1,allele1, i)*parentGLF[i][3];
  parentConditional[i][4] = likelihoodKids_denovo(allele1,allele2,allele1,allele2, i)*parentGLF[i][4];
  parentConditional[i][5] = likelihoodKids_denovo(allele1,allele2,allele2,allele2, i)*parentGLF[i][5];
  parentConditional[i][6] = likelihoodKids_denovo(allele2,allele2,allele1,allele1, i)*parentGLF[i][6];
  parentConditional[i][7] = likelihoodKids_denovo(allele2,allele2,allele1,allele2, i)*parentGLF[i][7];
  parentConditional[i][8] = likelihoodKids_denovo(allele2,allele2,allele2,allele2, i)*parentGLF[i][8];

  //Mariginal is joint likelihood of parents marginalizing all kids
  //This is to avoid redundant calculation
  parentMarginal[i][0] = parentConditional[i][0]*parentPrior[0];
  parentMarginal[i][1] = parentConditional[i][1]*parentPrior[1];
  parentMarginal[i][2] = parentConditional[i][2]*parentPrior[2];
  parentMarginal[i][3] = parentConditional[i][3]*parentPrior[3];
  parentMarginal[i][4] = parentConditional[i][4]*parentPrior[4];
  parentMarginal[i][5] = parentConditional[i][5]*parentPrior[5];
  parentMarginal[i][6] = parentConditional[i][6]*parentPrior[6];
  parentMarginal[i][7] = parentConditional[i][7]*parentPrior[7];
  parentMarginal[i][8] = parentConditional[i][8]*parentPrior[8];

}


void FamilyLikelihoodSeq_VCF::CalcPostProb_SingleNucFam(int i, double freq)
{
  double p11, p12, p22;
  double sum;
  int best;

  if(pedGLF->ped->families[i]->count<=pedGLF->ped->families[i]->founders)
  {
   for(int j=0; j<pedGLF->ped->families[i]->founders; j++)
   {
    sex = pedGLF->sexes[i][j];
    CalcPostProb_SinglePerson(i, j, freq);
   }
    return;
  }

      CalcParentMarginal(i, freq);

      for(int j=0; j<pedGLF->ped->families[i]->count; j++)
           {
            sex = pedGLF->sexes[i][j];
            if(j==0){
            p11 = parentMarginal[i][0] + parentMarginal[i][1] + parentMarginal[i][2];
            p12 = parentMarginal[i][3] + parentMarginal[i][4] + parentMarginal[i][5];
            p22 = parentMarginal[i][6] + parentMarginal[i][7] + parentMarginal[i][8];
            sum = p11 + p12 + p22;

            if(sum==0) 
                   postProb[i][j][0] = postProb[i][j][1] = postProb[i][j][2] = 1/3;
            else {
              postProb[i][j][0] = p11/sum;
              postProb[i][j][1] = p12/sum;
              postProb[i][j][2] = p22/sum;
            }
            best = GetBestGenoIdx(p11, p12, p22);
            bestGenoIdx[i][j] = best; 
            bestGenoLabel[i][j] = isChrY && pedGLF->sexes[i][j]==FEMALE ? "." : GetBestGenoLabel_vcfv4(best);
            //bestGenoLabel[i][j] = GetBestGenoLabel(best); 
            dosage[i][j] = CalcDosage(i,j); 
          } //first parent
          else if(j==1){
            p11 = parentMarginal[i][0] + parentMarginal[i][3] + parentMarginal[i][6];
            p12 = parentMarginal[i][1] + parentMarginal[i][4] + parentMarginal[i][7];
            p22 = parentMarginal[i][2] + parentMarginal[i][5] + parentMarginal[i][8];
            sum = p11 + p12 + p22;
            if(sum==0) 
                postProb[i][j][0] = postProb[i][j][1] = postProb[i][j][2] = 1/3;
           else {
              postProb[i][j][0] = p11/sum;
              postProb[i][j][1] = p12/sum;
              postProb[i][j][2] = p22/sum;
            }
            best = GetBestGenoIdx(p11, p12, p22);
            bestGenoIdx[i][j] = best; 
            bestGenoLabel[i][j] = isChrY && pedGLF->sexes[i][j]==FEMALE ? "." : GetBestGenoLabel_vcfv4(best); 
            //bestGenoLabel[i][j] = GetBestGenoLabel(best);
            dosage[i][j] = CalcDosage(i,j);

          } //second parent
          else {
            JointGenoLk JGLK = KidJointGenoLikelihood(i,j); //Joint likelihood of the kid j marginalizing other kids

            JGLK.CalcPost();

            postProb[i][j][0] = JGLK.post11;
            postProb[i][j][1] = JGLK.post12;
            postProb[i][j][2] = JGLK.post22;

            best = GetBestGenoIdx(JGLK.post11, JGLK.post12, JGLK.post22);

            bestGenoIdx[i][j] = best;

            //sex = pedGLF->sexes[i][j];

            bestGenoLabel[i][j] = (isChrY && pedGLF->sexes[i][j]==FEMALE) ? "." : GetBestGenoLabel_vcfv4(best); 
            //bestGenoLabel[i][j] = GetBestGenoLabel(best);
            dosage[i][j] = CalcDosage(i,j);
          }  //kids
      } //END of for loop
}

// Likelihood of genotypes of a kid marginalizing parents and his/her siblings
JointGenoLk FamilyLikelihoodSeq_VCF::KidJointGenoLikelihood(int famIdx, int kidIdx)
{ 
  double lkF11, lkF12, lkF22, lkM11, lkM12, lkM22;

  lkF11 = fam[famIdx].penetrances[0][geno11];
  lkF12 = fam[famIdx].penetrances[0][geno12];
  lkF22 = fam[famIdx].penetrances[0][geno22];
  lkM11 = fam[famIdx].penetrances[1][geno11];
  lkM12 = fam[famIdx].penetrances[1][geno12];
  lkM22 = fam[famIdx].penetrances[1][geno22];

  JointGenoLk JGLK1111 = likelihoodKidGenotype(allele1,allele1,allele1,allele1, famIdx, kidIdx);
  JointGenoLk JGLK1112 = likelihoodKidGenotype(allele1,allele1,allele1,allele2, famIdx, kidIdx);
  JointGenoLk JGLK1122 = likelihoodKidGenotype(allele1,allele1,allele2,allele2, famIdx, kidIdx);
  JointGenoLk JGLK1211 = likelihoodKidGenotype(allele1,allele2,allele1,allele1, famIdx, kidIdx);
  JointGenoLk JGLK1212 = likelihoodKidGenotype(allele1,allele2,allele1,allele2, famIdx, kidIdx);
  JointGenoLk JGLK1222 = likelihoodKidGenotype(allele1,allele2,allele2,allele2, famIdx, kidIdx);
  JointGenoLk JGLK2211 = likelihoodKidGenotype(allele2,allele2,allele1,allele1, famIdx, kidIdx);
  JointGenoLk JGLK2212 = likelihoodKidGenotype(allele2,allele2,allele1,allele2, famIdx, kidIdx);
  JointGenoLk JGLK2222 = likelihoodKidGenotype(allele2,allele2,allele2,allele2, famIdx, kidIdx);

  JGLK1111.multiplyParentLikelihood(parentGLF[famIdx][0]*parentPrior[0]);
  JGLK1112.multiplyParentLikelihood(parentGLF[famIdx][1]*parentPrior[1]);
  JGLK1122.multiplyParentLikelihood(parentGLF[famIdx][2]*parentPrior[2]);
  JGLK1211.multiplyParentLikelihood(parentGLF[famIdx][3]*parentPrior[3]);
  JGLK1212.multiplyParentLikelihood(parentGLF[famIdx][4]*parentPrior[4]);
  JGLK1222.multiplyParentLikelihood(parentGLF[famIdx][5]*parentPrior[5]);
  JGLK2211.multiplyParentLikelihood(parentGLF[famIdx][6]*parentPrior[6]);
  JGLK2212.multiplyParentLikelihood(parentGLF[famIdx][7]*parentPrior[7]);
  JGLK2222.multiplyParentLikelihood(parentGLF[famIdx][8]*parentPrior[8]);
  
  double JGLK11 = JGLK1111.g11 + JGLK1112.g11 + JGLK1122.g11 + JGLK1211.g11 + JGLK1212.g11 + JGLK1222.g11 + JGLK2211.g11 + JGLK2212.g11 + JGLK2222.g11;
  double JGLK12 = JGLK1111.g12 + JGLK1112.g12 + JGLK1122.g12 + JGLK1211.g12 + JGLK1212.g12 + JGLK1222.g12 + JGLK2211.g12 + JGLK2212.g12 + JGLK2222.g12;
  double JGLK22 = JGLK1111.g22 + JGLK1112.g22 + JGLK1122.g22 + JGLK1211.g22 + JGLK1212.g22 + JGLK1222.g22 + JGLK2211.g22 + JGLK2212.g22 + JGLK2222.g22;
  
  JointGenoLk jglk;
  jglk.g11 = JGLK11;
  jglk.g12 = JGLK12;
  jglk.g22 = JGLK22;

  return(jglk);
}


// Likelihood of genotypes of a kid marginalizing parents and (s)his siblings allowing for de novo mutations
JointGenoLk_denovo FamilyLikelihoodSeq_VCF::KidJointGenoLikelihood_denovo(int famIdx, int kidIdx)
{ 
  double lkF11, lkF12, lkF22, lkM11, lkM12, lkM22;

  lkF11 = fam[famIdx].penetrances[0][geno11];
  lkF12 = fam[famIdx].penetrances[0][geno12];
  lkF22 = fam[famIdx].penetrances[0][geno22];
  lkM11 = fam[famIdx].penetrances[1][geno11];
  lkM12 = fam[famIdx].penetrances[1][geno12];
  lkM22 = fam[famIdx].penetrances[1][geno22];

  JointGenoLk_denovo JGLK[9];
  JGLK[0] = likelihoodKidGenotype_denovo(allele1,allele1,allele1,allele1, famIdx, kidIdx);
  JGLK[1] = likelihoodKidGenotype_denovo(allele1,allele1,allele1,allele2, famIdx, kidIdx);
  JGLK[2] = likelihoodKidGenotype_denovo(allele1,allele1,allele2,allele2, famIdx, kidIdx);
  JGLK[3] = likelihoodKidGenotype_denovo(allele1,allele2,allele1,allele1, famIdx, kidIdx);
  JGLK[4] = likelihoodKidGenotype_denovo(allele1,allele2,allele1,allele2, famIdx, kidIdx);
  JGLK[5] = likelihoodKidGenotype_denovo(allele1,allele2,allele2,allele2, famIdx, kidIdx);
  JGLK[6] = likelihoodKidGenotype_denovo(allele2,allele2,allele1,allele1, famIdx, kidIdx);
  JGLK[7] = likelihoodKidGenotype_denovo(allele2,allele2,allele1,allele2, famIdx, kidIdx);
  JGLK[8] = likelihoodKidGenotype_denovo(allele2,allele2,allele2,allele2, famIdx, kidIdx);
  
  for(int i=0; i<9; i++)
    JGLK[i].multiplyParentLikelihood(parentGLF[famIdx][i]*parentPrior[i]);
  
  JointGenoLk_denovo jglk;  
  jglk.ResetGenoLk(0.0);
  
  for(int i=0; i<10; i++)
   for(int j=0; j<9; j++)
   {
    jglk.geno[i] += JGLK[j].geno[i];
   }
  return(jglk);
}


//likelihood of reads of all kids conditional on parental genotypes
double FamilyLikelihoodSeq_VCF::likelihoodKids(int fa1, int fa2, int ma1, int ma2, int famIdx)
{
  double lkKids = 1.0;
  double lk=1.0;
  double lk11, lk12, lk22; 
  int famSize = pedGLF->ped->families[famIdx]->count;
  for(int i=2; i<famSize; i++)
    {
      int sex = pedGLF->sexes[famIdx][i];
      lk = likelihoodONEKid(famIdx,i, fa1, fa2, ma1, ma2);
      lkKids *= lk;
    }
  return(lkKids);
}

//likelihood of reads of ONE kid conditional on parental genotypes
//key function for calculating the likelihood of a nuclear family
double FamilyLikelihoodSeq_VCF::likelihoodONEKid(int famIdx, int personIdx, int fa1, int fa2, int ma1, int ma2)
{ 
  double lk=1.0;
  double lk11, lk12, lk22;

   lk11 = fam[famIdx].penetrances[personIdx][geno11];
   lk12 = fam[famIdx].penetrances[personIdx][geno12];
   lk22 = fam[famIdx].penetrances[personIdx][geno22];

      // father 1/1
      if(fa1==allele1 && fa2==allele1 && ma1==allele1 && ma2==allele1)
        { lk = (isChrY && sex==FEMALE) ? 1.0: lk11; }
      else if(fa1==allele1 && fa2==allele1 && ma1==allele1 && ma2==allele2)
        {
        if(isChrX)
         lk = sex==MALE ? 0.5 * (lk11 + lk22) : 0.5 * (lk11 + lk12);
        else if(isChrY)
         lk = sex==MALE ? lk11 : 1.0;
        else if(isMT)
          lk = 0.5 * (lk11 + lk22);
        else lk = 0.5 * (lk11+lk12);
        }
      else if(fa1==allele1 && fa2==allele1 && ma1==allele2 && ma2==allele2)
        if(isChrX)
          lk = sex==MALE ? lk22 : lk12;
        else if(isChrY)
          lk = sex==MALE ? lk11 : 1.0;
        else if(isMT)
         lk = lk22;
        else lk = lk12;
      
      // father 1/2
      else if(fa1==allele1 && fa2==allele2 && ma1==allele1 && ma2==allele1)
        if(isChrX || isChrY || isMT)
          lk = 0.0;
        else lk = 0.5 * (lk11 + lk12);
      else if(fa1==allele1 && fa2==allele2 && ma1==allele1 && ma2==allele2) 
        if(isChrX || isChrY || isMT) lk = 0.0;
        else lk = 0.25*lk11 + 0.5*lk12 + 0.25*lk22;
      else if(fa1==allele1 && fa2==allele2 && ma1==allele2 && ma2==allele2 )
        if(isChrX || isChrY || isMT) lk = 0.0;
        else lk = 0.5 * (lk12 + lk22);       /*kidsCondLikelihood[famIdx][i][5] = lk; */
      
      // father 2/2
      else if(fa1==allele2 && fa2==allele2 && ma1==allele1 && ma2==allele1) 
        if(isChrX)
          lk = sex==MALE ? lk11 : lk12;
        else if(isChrY)
          lk = sex==MALE ? lk22 : 1.0;
        else if(isMT)
          lk = lk11;
        else lk = lk12;      /*kidsCondLikelihood[famIdx][i][6] = lk; */ 
      else if(fa1==allele2 && fa2==allele2 && ma1==allele1 && ma2==allele2)
        if(isChrX)
          lk = sex==MALE ? 0.5*(lk11+lk22) : 0.5 * (lk12 + lk22);
        else if(isChrY)
          lk = sex==MALE ? lk22 : 1.0;
        else if(isMT)
          lk = 0.5 * (lk11 + lk22);
        else lk = 0.5 * (lk12 + lk22);
      else if(fa1==allele2 && fa2==allele2 && ma1==allele2 && ma2==allele2) 
         lk = (isChrY && sex==FEMALE) ? 1.0 : lk22;       /*kidsCondLikelihood[famIdx][i][8] = lk;*/ 
      
  return(lk);
}

//likelihood of the reads and a specific genotype of a kid conditional on parental genotypes
//key function to modify for non-autosomes for calculating posterior genotype likelihoods
JointGenoLk FamilyLikelihoodSeq_VCF::likelihoodKidGenotype(int fa1, int fa2, int ma1, int ma2, int famIdx, int kidIndex)
{
  double lkKidG11 = 1.0;
  double lkKidG12 = 1.0;
  double lkKidG22 = 1.0;
  double lk=0.0;
  double lk11, lk12, lk22;
  double lkg11, lkg12, lkg22;

  int famSize = pedGLF->ped->families[famIdx]->count;

  for(int i=2; i<famSize; i++)
    {
      lk11 = fam[famIdx].penetrances[i][geno11];
      lk12 = fam[famIdx].penetrances[i][geno12];
      lk22 = fam[famIdx].penetrances[i][geno22];

      int sex = pedGLF->sexes[famIdx][i];      
      // father 1/1
      if(fa1==allele1 && fa2==allele1 && ma1==allele1 && ma2==allele1)
        {
        lk = lk11; 
        lkg11 = lk11; lkg12=lkg22=0; 
        }
      if(fa1==allele1 && fa2==allele1 && ma1==allele1 && ma2==allele2) 
        {
        if(isChrX)
         { lk = sex==MALE ? 0.5 * (lk11 + lk22) : 0.5 * (lk11+lk12);
           if(sex==MALE) {lkg11=0.5*lk11; lkg12=0.0; lkg22=0.5*lk22; } 
           else {lkg11=0.5*lk11; lkg12=0.5*lk12; lkg22=0;}
        }
        else if(isChrY)
         { lk = sex==MALE ? lk11 : 1.0; 
           if(sex==MALE){ lkg11=lk11; lkg12=lkg22=0.0;}
           else {lkg11=lkg12=lkg22=0.0;}
         }
        else if(isMT)
         { 
           lk = 0.5 * (lk11 + lk22);
           lkg11 = 0.5*lk11; lkg22 = 0.5*lk22; lkg12=0.0;
         }
        else 
         { lk = 0.5 * (lk11+lk12); lkg11=lk11*0.5; lkg12=lk12*0.5; lkg22=0;
         }
        }
      if(fa1==allele1 && fa2==allele1 && ma1==allele2 && ma2==allele2) 
        {
        if(isChrX){lk = sex==MALE ? lk22 : lk12; if(sex==MALE){lkg11=lkg12=0; lkg22=lk22;} else{lkg11=lkg22=0; lkg12=lk12;}}
        else if(isChrY){lk = sex==MALE ? lk11 : 1.0; if(sex==MALE){lkg11=lk11; lkg12=lkg22=0;} else{lkg11=lkg12=lkg22=0.;}}
        else if(isMT){lk = lk22; lkg11=lkg12=0; lkg22=lk22;}
        else { lk = lk12; lkg11=0; lkg12=lk12; lkg22=0; }
        }
      // father 1/2
      if(fa1==allele1 && fa2==allele2 && ma1==allele1 && ma2==allele1)
        { 
        if(isChrX || isChrY || isMT){lk=0.0; lkg11=lkg12=lkg22=0.0;}
        else {lk = 0.5 * (lk11 + lk12); lkg11=lk11*0.5; lkg12=lk12*0.5; lkg22=0; }
        }
      if(fa1==allele1 && fa2==allele2 && ma1==allele1 && ma2==allele2) 
        {
        if(isChrX || isChrY || isMT){lk=0.0; lkg11=lkg12=lkg22=0.0;}
        else {lk = 0.25*lk11 + 0.5*lk12 + 0.25*lk22; lkg11=lk11*0.25; lkg12=lk12*0.5; lkg22=lk22*0.25;}
        }
      if(fa1==allele1 && fa2==allele2 && ma1==allele2 && ma2==allele2)
        {
        if(isChrX || isChrY || isMT){lk=0.0; lkg11=lkg12=lkg22=0.0;}
        else{ lk = 0.5 * (lk12 + lk22); lkg11=0; lkg12=lk12*0.5; lkg22=lk22*0.5;}
        }
      
      // father 2/2
      if(fa1==allele2 && fa2==allele2 && ma1==allele1 && ma2==allele1) 
        {
        if(isChrX){ lk = sex==MALE ? lk11 : lk12; if(sex==MALE){lkg11=lk11; lkg12=lkg22=0.0;} else{lkg11=lkg22=0.0; lkg12=lk12;}}
        else if(isChrY){lk = sex==MALE ? lk22 : 1.0; if(sex==MALE){lkg11=lkg12=0.0; lkg22=lk22; }else{lkg11=lkg12=lkg22=0.0;}}
        else if(isMT) {lk=lk11; lkg11=lk11; lkg12=lkg22=0.0;}
        else {lk = lk12; lkg11=0; lkg12=lk12; lkg22=0; }
        }
      if(fa1==allele2 && fa2==allele2 && ma1==allele1 && ma2==allele2)
        {
        if(isChrX){lk = sex==MALE ? 0.5*(lk11+lk22) : 0.5*(lk12+lk22); if(sex==MALE){lkg11=0.5*lk11; lkg22=0.5*lk22; lkg12=0.0;} else{lkg11=0.0; lkg12=0.5*lk12; lkg22=0.5*lk22;} }
        else if(isChrY){lk = sex==MALE ? lk22 : 1.0; if(sex==MALE){lkg11=lkg12=0.0; lkg22=lk22; } else {lkg11=lkg12=lkg22=0.0;} }
        else if(isMT){lk=0.5*(lk11+lk22); lkg11=0.5*lk11; lkg22=0.5*lk22; lkg12=0.0;}
        else {lk = 0.5 * (lk12 + lk22); lkg11=0; lkg12=lk12*0.5; lkg22=lk22*0.5;}
        }
      if(fa1==allele2 && fa2==allele2 && ma1==allele2 && ma2==allele2)
        {
        if(isChrX){lk=lk22; lkg11=lkg11=0.0; lkg22=lk22;}
        if(isChrY){lk= sex==MALE ? lk22 : 1.0; if(sex==MALE){lkg11=lkg12=0.0; lkg22=lk22;} else {lkg11=lkg22=lkg12=0.0;}}
        if(isMT){lk=lk22; lkg11=lkg12=0.0; lkg22=lk22;}
        else {lk = lk22; lkg11=0; lkg12=0; lkg22=lk22;}
        }
        
      if(i!=kidIndex)
        {
          lkKidG11 *= lk;
          lkKidG12 *= lk;
          lkKidG22 *= lk;
        }
      else
        {
          lkKidG11*=lkg11;
          lkKidG12*=lkg12;
          lkKidG22*=lkg22;
        }
    }
  JointGenoLk JGLK;
  JGLK.g11 = lkKidG11;
  JGLK.g12 = lkKidG12;
  JGLK.g22 = lkKidG22;

  return(JGLK);
}


