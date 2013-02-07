#include "FamilyLikelihoodSeq.h"

FamilyLikelihoodSeq::FamilyLikelihoodSeq()
{
  fam = NULL;
}

FamilyLikelihoodSeq::~FamilyLikelihoodSeq()
{
  if(fam!=NULL)
    delete [] fam;
}

void FamilyLikelihoodSeq::InitFamilyLikelihoodES()
{
  fam = new FamilyLikelihoodES[nFam];
  for(int i=0; i<nFam; i++)
    {
//      if(ped->families[i]->count==ped->families[i]->founders) continue;

      fam[i].SetPedigree(ped);
      fam[i].SetFamily(ped->families[i]);
      fam[i].SetFamilyIndex(i);
      fam[i].SetTransmissionMatrix();
      fam[i].SetTransmissionMatrix_BA();
      fam[i].SetTransmissionMatrix_BA_CHRX_2Female();
      fam[i].SetTransmissionMatrix_BA_CHRX_2Male();
      fam[i].SetTransmissionMatrix_BA_CHRY();
      fam[i].SetTransmissionMatrix_BA_MITO();
      if(par->denovo){ fam[i].SetGenotypeMutationModel(&gM);
        fam[i].SetTransmissionMatrix_denovo(); }
      if(ped->families[i]->count!=ped->families[i]->founders) fam[i].PreparePeeling();
//fam[i].es.PrintPeelable();
//fam[i].es.PrintPeelingOrder();

    }
}

double FamilyLikelihoodSeq::f(double freq)
{
  return(-(CalcAllFamLogLikelihood(freq)));
}

void FamilyLikelihoodSeq::BackupFounderCount()
{
 if(backupFounderCount.size()==0)
 {
  backupFounderCount.resize(pedGLF->ped->familyCount);
  for(int i=0; i < pedGLF->ped->familyCount; i++)
   backupFounderCount[i] = pedGLF->ped->families[i]->founders;
 }
}

void FamilyLikelihoodSeq::MakeUnrelated()
{
 if(backupFounderCount.size()==0) error("Founder count not backed up yet!\n");
 for(int i=0; i < pedGLF->ped->familyCount; i++)
   pedGLF->ped->families[i]->founders = pedGLF->ped->families[i]->count;
}

void FamilyLikelihoodSeq::RestoreFounderCount()
{
 if(backupFounderCount.size()==0) error("backupFounderCount has not been allocated memory");
 for(int i=0; i < pedGLF->ped->familyCount; i++)
  pedGLF->ped->families[i]->founders = backupFounderCount[i];
}

double FamilyLikelihoodSeq::MonomorphismLogLikelihood_denovo(int refBase, int alt)
{
 SetAlleles(refBase, alt);
 return(CalcAllFamLogLikelihood(1.0)); 
}

void FamilyLikelihoodSeq::CalcPostProb(double freq)
{
  for(int i=0; i<nFam; i++)
    {
      if(ped->families[i]->isNuclear() || ped->families[i]->count==ped->families[i]->founders)
	{
	  if(par->denovo) CalcPostProb_SingleNucFam_denovo(i, freq);
	  else CalcPostProb_SingleNucFam(i, freq);
	}
      else 
	{
	  if(par->denovo) CalcPostProb_SingleExtendedPed_denovo(i, freq);
	  else CalcPostProb_SingleExtendedPed_BA(i, freq);
	}
    }
}

double FamilyLikelihoodSeq::PolymorphismLogLikelihood(int  a1, int a2)
{
  SetAlleles(a1, a2);
  if(nFam>1 || nFam==1 && !ped->families[0]->isNuclear())
    {
      OptimizeFrequency();
      return(GetMaxLogLikelihood());
    }
  else // Single nuclear family and prior is fixed for two parents
  {
     double dummy = 0.5;
     return(CalcAllFamLogLikelihood(dummy));
  }
}
                                    
                                    
void FamilyLikelihoodSeq::CalcPostProb_SingleExtendedPed(int i, double freq)
{
  double lk[10];
  
  for(int j=0; j<pedGLF->ped->families[i]->count; j++)
    {
      for(int k=0; k<10; k++)
	{
	  FillZeroPenetrance(&fam[i], pedGLF, j, k);
	  lk[k] = CalcSingleFamLikelihood(i, freq);
	}
      
      double sum = 0.0;
      for(int k=0; k<10; k++)
        sum += lk[k];
      
      if(sum==0)
	for(int k=0; k<10; k++) postProb[i][j][k] = 0;
      else
	for(int k=0; k<10; k++) postProb[i][j][k] = lk[k]/sum;
      
      int best = 0; double max = 0.0;
      for(int k=0; k<10; k++)
	if(max<lk[k]) { max=lk[k]; best=k; }

      bestGenoIdx[i][j] = best;       
      bestGenoLabel[i][j] = GetBestGenoLabel_denovo(best); //This is correct for the 10 genotpye case
      
    }
    

}

void FamilyLikelihoodSeq::CalcPostProb_SingleExtendedPed_denovo(int i, double freq)
{
  double lk[10];
  
  for(int j=0; j<pedGLF->ped->families[i]->count; j++)
    {
      for(int k=0; k<10; k++)
	{
	  FillZeroPenetrance(&fam[i], pedGLF, j, k);
	  lk[k] = CalcSingleFamLikelihood_denovo(i, freq);
	}
      
      double sum = 0.0;
      for(int k=0; k<10; k++)
	sum += lk[k];
      
      if(sum==0)
	for(int k=0; k<10; k++) postProb[i][j][k] = 0;
      else
	for(int k=0; k<10; k++) postProb[i][j][k] = lk[k]/sum;
      
      int best = 0; double max = 0.0;
      for(int k=0; k<10; k++) {
	if(max<lk[k]) { max=lk[k]; best=k;}
    }
      bestGenoIdx[i][j] = best;      
      bestGenoLabel[i][j] = GetBestGenoLabel_denovo(best);

    }
}

void FamilyLikelihoodSeq::CalcPostProb_SingleExtendedPed_BA(int i, double freq)
{
  double lk11, lk12, lk22;
  double sum;
  int best;
  
  for(int j=0; j<pedGLF->ped->families[i]->count; j++)
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

      FillZeroPenetrance(&fam[i], pedGLF, j, glfHandler::GenotypeIndex(allele1, allele1));
      lk11 = CalcSingleFamLikelihood_BA(i, freq);
      
      FillZeroPenetrance(&fam[i], pedGLF, j, glfHandler::GenotypeIndex(allele1, allele2));
      lk12 = CalcSingleFamLikelihood_BA(i, freq);
      
      FillZeroPenetrance(&fam[i], pedGLF, j, glfHandler::GenotypeIndex(allele2, allele2));
      lk22 = CalcSingleFamLikelihood_BA(i, freq);
      
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
      bestGenoLabel[i][j] = GetBestGenoLabel_vcfv4(best);
      dosage[i][j] = CalcDosage(i,j);
    }
}

double FamilyLikelihoodSeq::CalcAllFamLikelihood(double freq)
{
}

double FamilyLikelihoodSeq::CalcAllFamLogLikelihood(double freq)
{
  double loglk = 0.0;
#pragma omp parallel for reduction(+:loglk)
  for(int i=0; i<nFam; i++)
    {
      if(ped->families[i]->isNuclear() || ped->families[i]->count==ped->families[i]->founders)
{	loglk += par->denovo ? logLkSingleFam_denovo(i, freq) : logLkSingleFam(i, freq); 
//	if(par->denovo) printf("Nuc fam de novo %d lk=%f\n", i, logLkSingleFam_denovo(i, freq) );
//	else printf("Nuc fam %d lk=%f\n", i, logLkSingleFam(i, freq) ); 
}
      else
{	loglk +=  par->denovo ? CalcSingleFamLogLikelihood_denovo(i, freq) :  CalcSingleFamLogLikelihood_BA(i, freq); 
//	if(par->denovo) printf("Ext fam de novo %d lk=%f\n", i, CalcSingleFamLogLikelihood_denovo(i, freq) );
//	else printf("Ext fam %d lk=%f\n", i, CalcSingleFamLogLikelihood_BA(i, freq) );
 }
    }
  return(loglk);
}

// key function and most of calculations are here. Note: this is for the 10 genotype case and not used yet
double FamilyLikelihoodSeq::CalcSingleFamLikelihood(int i, double freq)
{
  double lk=1.0;
  
  fam[i].SetAlleles(allele1,allele2);
  fam[i].SetFounderPriors(freq);
  fam[i].InitializePartials();
  lk = fam[i].CalculateLikelihood();

  return(lk);
}

// for bi-allelic markers for computation speedup
double FamilyLikelihoodSeq::CalcSingleFamLikelihood_BA(int i, double freq)
{
  double lk=1.0;
  
  fam[i].SetAlleles(allele1,allele2);
  fam[i].SetFounderPriors_BA(freq);
  fam[i].InitializePartials_BA();
  lk = fam[i].CalculateLikelihood_BA();

  return(lk);
}

// allowing for de novo mutations
double FamilyLikelihoodSeq::CalcSingleFamLikelihood_denovo(int i, double freq)
{
  double lk=1.0;

  fam[i].SetAlleles(allele1,allele2);
  fam[i].SetFounderPriors(freq); //Note: we still assume bi-allelic variants founders
  fam[i].InitializePartials();
  lk = fam[i].CalculateLikelihood_denovo();

  return(lk);
}

double FamilyLikelihoodSeq::CalcSingleFamLogLikelihood(int i, double freq)
{
  return(log10(CalcSingleFamLikelihood(i, freq)));
}

double FamilyLikelihoodSeq::CalcSingleFamLogLikelihood_BA(int i, double freq)
{
  return(log10(CalcSingleFamLikelihood_BA(i, freq)));
}

double FamilyLikelihoodSeq::CalcSingleFamLogLikelihood_denovo(int i, double freq)
{
  return(log10(CalcSingleFamLikelihood_denovo(i, freq)));
}

void FamilyLikelihoodSeq::FillPenetrance()
{
 for(int i=0; i<nFam; i++)
  FillPenetrance(&fam[i], pedGLF);
}

void  FamilyLikelihoodSeq::FillPenetrance(FamilyLikelihoodES* famlk, PedigreeGLF* pedGLF)
{
  famlk->penetrances.Zero();
  for(int i=0; i<famlk->famSize;i++)
    { 
      if(pedGLF->glf[famlk->famIdx][i].handle==NULL)
	{
	  for(int j=0; j<10; j++)
	    famlk->penetrances[i][j] = 1.0;
	  continue;
	}
      
   for(int j=0; j<10; j++)
	  famlk->penetrances[i][j] = pedGLF->glf[famlk->famIdx][i].GetLikelihoods(pedGLF->currentPos)[j];
    }
}

void  FamilyLikelihoodSeq::FillPenetrance(FamilyLikelihoodSeq& sourceFamLk)
{
  for(int f=0; f<nFam; f++)
   for(int i=0; i<fam[f].famSize;i++)
    for(int j=0; j<10; j++)
	  fam[f].penetrances[i][j] = sourceFamLk.fam[f].penetrances[i][j];
}

void FamilyLikelihoodSeq::FillZeroPenetrance(FamilyLikelihoodES *famlk, PedigreeGLF *pedGLF, int person, int genoIdx)
{
  for(int i=0; i<famlk->famSize;i++)
    { 
      if(pedGLF->glf[famlk->famIdx][i].handle==NULL)
	{
	  if(i!=person)
	    for(int j=0; j<10; j++)
	      famlk->penetrances[i][j] = 1;
	  else
	    for(int j=0; j<10; j++)
	      if(j==genoIdx)
		famlk->penetrances[i][j] = 1.0;
	      else famlk->penetrances[i][j] = 0.0;

	  continue;
	}
      
  if(i!=person)
	for(int j=0; j<10; j++)
	  famlk->penetrances[i][j] = pedGLF->glf[famlk->famIdx][i].GetLikelihoods(pedGLF->currentPos)[j];
  else
  {
	for(int j=0; j<10; j++)
	  if(j==genoIdx)
	    famlk->penetrances[i][j] = pedGLF->glf[famlk->famIdx][i].GetLikelihoods(pedGLF->currentPos)[j];
	  else famlk->penetrances[i][j] = 0.0;
   }   
 }
}

void FamilyLikelihoodSeq::SetNonAutosomeFlags(bool x, bool y, bool mt)
{
 if(x && y || x && mt || y && mt) error("Only one non-autosomal chrom can be specified\n");

 isChrX = x; isChrY=y; isMT=mt;
  for(int i=0; i<nFam; i++) 
  {
   fam[i].isChrX = x;
   fam[i].isChrY=y;
   fam[i].isMT=mt;
  }
}
