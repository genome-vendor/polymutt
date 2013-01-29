#include "NucFamGenotypeLikelihood.h"
#include "StringMap.h"
#include <stdlib.h>
#include <iostream>
#include "CmdLinePar.h"
#include "MutationModel.h"
#include <map>
#include "FamilyLikelihoodES.h"
#include "FamilyLikelihoodSeq.h"
#include "PedVCF.h"

StringArray tokens;
StringArray lists;

int readGLFannoFile(String GLFannoFile, StringMap *glfMap)
{
  IFILE fh = ifopen(GLFannoFile.c_str(), "r");
  if(fh==NULL)
    error("%s open failed\n", GLFannoFile.c_str());
  
  String line;
  StringArray tokens;
  
  int count = 0;
  while(!ifeof(fh))
    {
      line.ReadLine(fh);
      tokens.Clear();
      tokens.ReplaceTokens(line);
      if(tokens.Length()<2) { continue; }
      lists.Add(tokens[1]);
      glfMap->Add(tokens[0],&lists[count]);
      count++;
    }
  ifclose(fh);
  return(count);
}

int LoadPositionFile(String &file, std::map<String, int>& positionMap)
{
  FILE *fh = fopen(file.c_str(), "r");
  if(fh==NULL) error("Open position file %s failed!\n", file.c_str());
  String buffer;
  StringArray tokens;
  while(!feof(fh))
    {
      buffer.ReadLine(fh);
      tokens.ReplaceTokens(buffer);
      if(tokens.Length()==0) continue;
      String chrPos = tokens[0]+":"+tokens[1];
      positionMap[chrPos]++;     
    }
  fclose(fh);
  return(0);
}

int main(int argc, char * argv[])
{
  double posterior = 0.5;
  int minTotalDepth = 0;
  int maxTotalDepth = 0;
  double minPS = 0;
  int minMapQuality = 0;
  String pedFile, datFile, glfListFile;
  String vcfOutFile = "";
  String vcfInFile ="";
  String positionfile;
  double theta = 0.001;
  double theta_indel = 0.0001;
  double tstv_ratio = 2.0;
  double precision = 0.0001;
  int num_threads = 1;
  bool denovo = false;
  double denovo_mut_rate = 1.5e-08;
  double denovo_tstv_ratio = 2.0;
  double denovoLR = 0.01;
  bool gl_off = false;
  String chrs2process;
  String chrX_label("X");
  String chrY_label("Y");
  String MT_label("MT");
  bool quick_call = false;
  bool use_ext = false;
  bool force_call = false;
  bool out_all_sites = false;
  int out_cnt = 0;

  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("Alternative input file")
    LONG_STRINGPARAMETER("in_vcf", &vcfInFile)
    LONG_PARAMETER_GROUP("Scaled mutation rate")
    LONG_DOUBLEPARAMETER("theta", &theta)
    LONG_DOUBLEPARAMETER("indel_theta", &theta_indel)
    LONG_PARAMETER_GROUP("Prior of ts/tv ratio")
    LONG_DOUBLEPARAMETER("poly_tstv", &tstv_ratio)
  LONG_PARAMETER_GROUP("Non-autosome labels")
      LONG_STRINGPARAMETER("chrX", &chrX_label)
      LONG_STRINGPARAMETER("chrY", &chrY_label)
      LONG_STRINGPARAMETER("MT", &MT_label)

  LONG_PARAMETER_GROUP("de novo mutation")
      LONG_PARAMETER("denovo", &denovo)
      LONG_DOUBLEPARAMETER("rate_denovo", &denovo_mut_rate)
      LONG_DOUBLEPARAMETER("tstv_denovo", &denovo_tstv_ratio)
      LONG_DOUBLEPARAMETER("minLLR_denovo", &denovoLR)
    LONG_PARAMETER_GROUP("Optimization precision")
    LONG_DOUBLEPARAMETER("prec", &precision)
  LONG_PARAMETER_GROUP("Multiple threading")
    LONG_INTPARAMETER("nthreads", &num_threads)
//  LONG_PARAMETER_GROUP("Chromosomes to process")
//    LONG_STRINGPARAMETER("chr2process", &chrs2process)
    LONG_PARAMETER_GROUP("Filters")
    LONG_INTPARAMETER("minMapQuality", &minMapQuality)
    LONG_INTPARAMETER("minDepth", &minTotalDepth)
    LONG_INTPARAMETER("maxDepth", &maxTotalDepth)
    LONG_DOUBLEPARAMETER("minPercSampleWithData", &minPS)
  LONG_PARAMETER_GROUP("Output")
    LONG_STRINGPARAMETER("out_vcf", &vcfOutFile)
    LONG_STRINGPARAMETER("pos", &positionfile)
    LONG_PARAMETER("all_sites", &out_all_sites)
    LONG_PARAMETER("gl_off", &gl_off)
    LONG_PARAMETER("quick_call", &quick_call)
    //LONG_PARAMETER("ext", &use_ext)

    END_LONG_PARAMETERS();
  
  pl.Add(new StringParameter('p', "pedfile", pedFile));
  pl.Add(new StringParameter('d', "datfile", datFile));
  pl.Add(new StringParameter('g', "glfIndexFile", glfListFile));
  pl.Add(new DoubleParameter('c', "posterior cutoff", posterior));
  pl.Add(new LongParameters("Additional Options", longParameters));
  pl.Read(argc, argv);

  pl.Status();


  if(vcfInFile==vcfOutFile) error("Input and output VCF files are the same!\n");

  if(pedFile.Length()==0)
    error("pedFile not provided for input!\n");
//  if(datFile.Length()==0)
//    error("datFile not provided for input!\n");
  if(glfListFile.Length()==0 && vcfInFile.Length()==0)
    error("glfListFile or input VCF file not provided for input!\n");
  if(vcfOutFile.Length()==0)
    error("vcfOutFile not provided for output!\n");

  std::map<String, int> positionMap;
  if(positionfile.Length()>0) { LoadPositionFile(positionfile, positionMap); force_call = true; quick_call = false;out_all_sites=false;}

  if(out_all_sites) quick_call = false;

  #ifdef _OPENMP
  if(num_threads>0) omp_set_num_threads(num_threads);
  #endif

  std::string cmd;
  for(int a=0; a<argc; a++)
  {
   cmd += std::string(argv[a]);
   cmd += " ";
  }

  CmdLinePar par;
  par.cmd = cmd;
  par.theta = theta;
  par.theta_indel = theta_indel;
  par.minTotalDepth = minTotalDepth;
  par.maxTotalDepth = maxTotalDepth;
  par.minMapQuality = minMapQuality;
  par.minPS = minPS;
  par.posterior     = posterior;
  par.precision     = precision;
  par.denovo_mut_rate = denovo_mut_rate;
  par.denovo_tstv_ratio = denovo_tstv_ratio;
  par.denovo = denovo;
  par.denovoLR = denovoLR;
  par.gl_off = gl_off;
  par.chrX_label = chrX_label;
  par.chrY_label = chrY_label;
  par.MT_label = MT_label;
  par.vcfInFile = vcfInFile;
  par.vcfOutFile = vcfOutFile;
  par.force_call = force_call;
  par.out_all_sites = out_all_sites;

  if(denovo && denovoLR<0) error("denovo_min_LLR can only be greater than 0 !\n");

  // Prior of ts and tv
  double prior_ts = tstv_ratio/(tstv_ratio + 1);
  double prior_tv = (1-prior_ts)/2;

  double polyPrior, polyPrior_unr;
  double refTransition = .0;
  double refTransvers1 = .0;
  double refTransvers2 = .0;
  double tstvs1=0; double tstvs2=.0; double tvs1tvs2=.0;
  double refTransitionFreq, refTransvers1Freq, refTransvers2Freq, tstvs1Freq, tstvs2Freq, tvs1tvs2Freq;

  StringMap glfMap;
  String glfFileKey;

  if(vcfInFile.Length()==0) readGLFannoFile(glfListFile, &glfMap);

  Pedigree ped;
  PedigreeGLF pedGLF;

  if(vcfInFile.Length()>0) pedGLF.SetVCFInputFlag(true);

  IFILE datFH = ifopen(datFile, "r");
  IFILE pedFH = ifopen(pedFile, "r");
  FILE *vcfFH = fopen(vcfOutFile, "w");
  if(datFH==NULL)
    error("datFile open for input failed!\n");
  if(pedFH==NULL)
    error("pedFile open for input failed!\n");
  if(vcfFH==NULL)
    error("vcfOutFile can not be opened for output!\n");

  ped.Prepare(datFH);
  ped.Load(pedFH);

//  This is for debug purpose which is to compare results of nuclear families
//      using two algorithms. They should give IDENTICAL results
 if(use_ext)
  for(int f=0; f<ped.familyCount; f++)
   ped.families[f]->generations=3;


  if(datFH != NULL) ifclose(datFH);
  if(pedFH != NULL) ifclose(pedFH);

  pedGLF.SetGLFMap(&glfMap);
  pedGLF.SetPedGLF(&ped);

  if(vcfInFile.Length()>0){
   PedVCF pv;
   pv.SetPedigreeGLF(&pedGLF);
   pv.SetPedigreePtr(&ped);
   pv.SetPar(par);
   pv.SetVCFInput(vcfInFile);
   pv.VarCallFromVCF();
   exit(0);
  }

  FamilyLikelihoodSeq famlk[7];
  for(int i=0; i<7; i++)
    {
      famlk[i].SetCmdLinePar(&par);
      if(denovo) famlk[i].SetDenovoMutationModel();
      famlk[i].SetTheta(theta);
      famlk[i].SetTheta_indel(theta_indel);
      famlk[i].SetGLF(&pedGLF);
      famlk[i].InitFamilyLikelihoodES();
      if(quick_call) famlk[i].BackupFounderCount();
    }

  int cnt=0;
  int cntSec=0;
  int totalEntryCnt = 0;

  //sumarry statistics
  uint minTotalDepthFilter = 0;
  uint maxTotalDepthFilter = 0;
  uint maxAvgDepthFilter = 0;
  uint minAvgDepthFilter = 0;
  uint minMapQualFilter = 0;
  uint minAvgMapQualFilter = 0;
  uint minPSFilter = 0;

  int refBaseCounts[5] = {0,0,0,0,0};
  int homoRef = 0;
  int transitions = 0;
  int transversions = 0;
  int otherPolymorphism = 0;
  int tstvs1Cnt = 0;
  int tstvs2Cnt = 0;
  int tvs1tvs2Cnt = 0;
  int nocall = 0;
  int actuaBases = 0;

  time_t t; time(&t);

  std::map<String, int> chrs2process_map;
  StringArray chrs;
  chrs.AddTokens(chrs2process, ',');
  for(int cidx=0; cidx<chrs.Length(); cidx++)
   chrs2process_map[chrs[cidx]]++;

  int chrs2processCount = chrs2process_map.size();

  printf("Analysis started on %s\n", ctime(&t));
  Matrix pen; 
  int maxidx = 0;
  int chrProcessedCount=0;
  String invalid_pid; int invalid_fam_idx, invalid_person_idx;

  while(pedGLF.Move2NextSection())
    {
      if(chrs2process_map.size()>0 && chrProcessedCount >= chrs2processCount) break;

//      if(!pedGLF.CheckSectionLabels(invalid_pid, invalid_fam_idx, invalid_person_idx)) {
//	fprintf(stderr, "Error: GLF of person with PID %s has invalid Section label '%s'\n", invalid_pid.c_str(), pedGLF.glf[invalid_fam_idx][invalid_person_idx].label.c_str());
//	exit(1);}

      if(chrs2process_map.size()>0 && chrs2process_map[pedGLF.GetNonNULLglf()->label]<1) { while(pedGLF.Move2NextBaseEntry()){}; continue;}

      //printf("Processing reference %s ...\n", pedGLF.GetNonNULLglf()->label.c_str());

      if(pedGLF.GetNonNULLglf()->label == par.chrX_label) { for(int i=0; i<7; i++) famlk[i].SetNonAutosomeFlags(true, false, false); }
      else if(pedGLF.GetNonNULLglf()->label == par.chrY_label) { for(int i=0; i<7; i++) famlk[i].SetNonAutosomeFlags(false, true, false); }
      else if(pedGLF.GetNonNULLglf()->label == par.MT_label) {  for(int i=0; i<7; i++)  famlk[i].SetNonAutosomeFlags(false, false, true); }
      else { for(int i=0; i<7; i++) famlk[i].SetNonAutosomeFlags(false, false, false);  }

      homoRef = transitions = transversions = otherPolymorphism = tstvs1Cnt = tstvs2Cnt = tvs1tvs2Cnt = nocall = actuaBases = 0; for(int k=0; k<5; k++) refBaseCounts[k]=0; 
      cnt=cntSec=totalEntryCnt = 0; minTotalDepthFilter = maxTotalDepthFilter = maxAvgDepthFilter = minAvgDepthFilter = minMapQualFilter = minAvgMapQualFilter = 0;

      polyPrior = famlk[0].GetPolyPrior();
      polyPrior_unr = famlk[0].GetPolyPrior_unr();

      chrProcessedCount++;

      while(pedGLF.Move2NextBaseEntry())
	{
	  for(int r=0; r<7; r++)
	    famlk[r].FillPenetrance();

	  if(totalEntryCnt==0) totalEntryCnt = pedGLF.GetNonNULLglf()->maxPosition;

	  if(positionfile.Length()>0)
	    {
	      int pos = pedGLF.currentPos+1;
	      String chrPos = pedGLF.GetNonNULLglf()->label+":"+pos;
	      if(positionMap.count(chrPos)==0) continue;
	    }

	  int refBase  = pedGLF.GetRefBase();
	  if(refBase!=1 && refBase!=2 && refBase!=3 && refBase!=4) continue;
	  refBaseCounts[refBase]++;

	  famlk[0].CalcReadStats();

	  if(famlk[0].totalDepth<minTotalDepth) { minTotalDepthFilter++; continue; }
	  if(maxTotalDepth>0 && famlk[0].totalDepth>maxTotalDepth) { maxTotalDepthFilter++; continue; }
	  if(famlk[0].percSampWithData*100<minPS){ minPSFilter++; continue; }
	  if(famlk[0].avgMapQual<minMapQuality) { minMapQualFilter++; continue; }

	  int ts       = (Poly::ts(refBase));   //transition
	  int tvs1     = (Poly::tvs1(refBase)); //transvertions1
	  int tvs2     = (Poly::tvs2(refBase)); //transvertion2

if(quick_call)
{
for(int r=0; r<7; r++)
   famlk[r].MakeUnrelated();

  double lRef_unr;double refTransition_unr ;double refTransvers1_unr ; double refTransvers2_unr;

# ifdef _OPENMP
# pragma omp parallel sections
# endif
{
# ifdef _OPENMP
#pragma omp section
# endif
	 {
	lRef_unr = log10(1-polyPrior_unr) + famlk[0].MonomorphismLogLikelihood(refBase);
	famlk[0].varllk[0]  = lRef_unr;
	}
# ifdef _OPENMP
# pragma omp section
# endif
	{
	refTransition_unr = log10(polyPrior_unr * prior_ts) + famlk[1].PolymorphismLogLikelihood(refBase, ts);
	famlk[0].varllk[1]  = refTransition_unr;
	}
# ifdef _OPENMP
# pragma omp section
# endif
	{
	refTransvers1_unr = log10(polyPrior_unr * prior_tv) + famlk[2].PolymorphismLogLikelihood(refBase, tvs1);
	famlk[0].varllk[2]  = refTransvers1_unr;
	}
# ifdef _OPENMP
# pragma omp section
# endif
	{
	refTransvers2_unr = log10(polyPrior_unr * prior_tv) + famlk[3].PolymorphismLogLikelihood(refBase, tvs2);
	famlk[0].varllk[3]  = refTransvers2_unr;
	}
}

	maxidx = famlk[0].CalcVarPosterior(4);

       // Calculate likelihoods for less likely SNP configurations
	double tstvs1_unr; double tstvs2_unr ; double tvs1tvs2_unr ;
	if(famlk[0].varPostProb<0.99)
	{
# ifdef _OPENMP
# pragma omp parallel sections
# endif
{
# ifdef _OPENMP
#pragma omp section
# endif
	{
	tstvs1_unr   = log10(polyPrior_unr * 0.001) + famlk[4].PolymorphismLogLikelihood(ts, tvs1);
	famlk[0].varllk[4]  = tstvs1_unr;
	}
# ifdef _OPENMP
#pragma omp section
# endif
	{
	 tstvs2_unr   = log10(polyPrior_unr * 0.001) + famlk[5].PolymorphismLogLikelihood(ts, tvs2);
	 famlk[0].varllk[5]  = tstvs2_unr;
	}
# ifdef _OPENMP
#pragma omp section
# endif
	{
	 tvs1tvs2_unr = log10(polyPrior_unr * 0.001) + famlk[6].PolymorphismLogLikelihood(tvs1, tvs2);
	 famlk[0].varllk[6]  = tvs1tvs2_unr;
	}
}
	 maxidx = famlk[0].CalcVarPosterior(7);
} // end of less likely configurations

//printf("%f - %d: %f %f %f %f %f %f %f : %f : %f\n", famlk[0].varPostProb, pedGLF.currentPos+1, lRef_unr, refTransition_unr, refTransvers1_unr, refTransvers2_unr, tstvs1_unr, tstvs2_unr, tvs1tvs2_unr, polyPrior_unr, tstvs1Freq = famlk[2].GetMinimizer());

	  if(famlk[0].varPostProb<posterior) continue; // no call due to lack of evidence of calling anything
	  if(maxidx==0) continue;  // monomorphism and skip


	for(int r=0; r<7; r++) famlk[r].RestoreFounderCount();
} // end of quick variant calling

# ifdef _OPENMP
# pragma omp parallel sections
# endif
	  {
# ifdef _OPENMP
#pragma omp section
# endif
	    {
	      if(par.denovo==false)
	      {
		//Calculate likelihood assuming everyone is homozygous for the reference
		double lRef = log10(1-polyPrior) + famlk[0].MonomorphismLogLikelihood(refBase);
		famlk[0].varllk[0]  = lRef;
		famlk[0].varllk_noprior[0] = lRef - log10(1-polyPrior);
		famlk[0].varfreq[0] = 1.0;
	      }
	      else 
	      {
		//Calculate likelihood assuming everyone is homozygous for the reference allowing for denovo mutations
		double lRef_denovo = log10(1-polyPrior) + famlk[0].MonomorphismLogLikelihood_denovo(refBase, refBase==4 ? refBase-1 : refBase+1);
		famlk[0].varllk[0] = lRef_denovo;
		famlk[0].varllk_noprior[0] = lRef_denovo - log10(1-polyPrior);
		famlk[0].varfreq[0] = 1.0;
	      }
	    }
# ifdef _OPENMP
# pragma omp section
# endif
	    {
	      // Calculate likelihoods for the most likelily SNP configurations
	      refTransition = log10(polyPrior * prior_ts) + famlk[1].PolymorphismLogLikelihood(refBase, ts);
	      refTransitionFreq = famlk[1].GetMinimizer();
	      famlk[0].varllk[1]  = refTransition;
		   famlk[0].varllk_noprior[1] = refTransition -  log10(polyPrior * 2./3. );
	      famlk[0].varfreq[1] = refTransitionFreq;
	    }
# ifdef _OPENMP
# pragma omp section
# endif
	    {
	      refTransvers1 = log10(polyPrior * prior_tv) + famlk[2].PolymorphismLogLikelihood(refBase, tvs1);
	      refTransvers1Freq = famlk[2].GetMinimizer();
	      famlk[0].varllk[2]  = refTransvers1;
	    	famlk[0].varllk_noprior[2] = refTransvers1 - log10(polyPrior * 1./6.);
	      famlk[0].varfreq[2] = refTransvers1Freq;
	    }
# ifdef _OPENMP
# pragma omp section
# endif
	    {
	      refTransvers2 = log10(polyPrior * prior_tv) + famlk[3].PolymorphismLogLikelihood(refBase, tvs2);
	      refTransvers2Freq = famlk[3].GetMinimizer();
	      famlk[0].varllk[3]  = refTransvers2;
	      famlk[0].varllk_noprior[3] = refTransvers2 - log10(polyPrior * 1./6.);
	      famlk[0].varfreq[3] = refTransvers2Freq;
	    }
	  }
	  maxidx = famlk[0].CalcVarPosterior(4);

	  // Calculate likelihoods for less likely SNP configurations
	  if(famlk[0].varPostProb<0.99)
	    {
# ifdef _OPENMP
#pragma omp parallel sections
# endif
{
# ifdef _OPENMP
#pragma omp section
# endif
		{
		  tstvs1   = log10(polyPrior * 0.001) + famlk[4].PolymorphismLogLikelihood(ts, tvs1);
		  tstvs1Freq = famlk[4].GetMinimizer();
		  famlk[0].varllk[4]  = tstvs1;
		  famlk[0].varllk_noprior[4] = tstvs1 -  log10(polyPrior * 0.001);
		  famlk[0].varfreq[4] = tstvs1Freq;
		}
# ifdef _OPENMP
#pragma omp section
# endif
		{
		  tstvs2   = log10(polyPrior * 0.001) + famlk[5].PolymorphismLogLikelihood(ts, tvs2);
		  tstvs2Freq = famlk[5].GetMinimizer();
		  famlk[0].varllk[5]  = tstvs2;
		  famlk[0].varllk_noprior[5] = tstvs2 - log10(polyPrior * 0.001) ;
		  famlk[0].varfreq[5] = tstvs2Freq;
		}
# ifdef _OPENMP
#pragma omp section
# endif
		{
		  tvs1tvs2 = log10(polyPrior * 0.001) + famlk[6].PolymorphismLogLikelihood(tvs1, tvs2);
		  tvs1tvs2Freq = famlk[6].GetMinimizer();
		  famlk[0].varllk[6]  = tvs1tvs2;
		  famlk[0].varllk_noprior[6] = tvs1tvs2 - log10(polyPrior * 0.001);
		  famlk[0].varfreq[6] = tvs1tvs2Freq;
		}
	      }
	      maxidx = famlk[0].CalcVarPosterior(7);
	    }

	  if(famlk[0].varPostProb<posterior) { nocall++; if(force_call==false && out_all_sites==false) continue; }

	  switch(maxidx)
	    {
	    //case 0: { homoRef++; famlk[0].SetAlleles(refBase, refBase==4?refBase-1:refBase+1); famlk[0].min=1.0; break;}
	    case 0: { homoRef++; if(force_call || out_all_sites) famlk[0].min=1.0; break;}
	    case 1: { transitions++; famlk[0].SetAlleles(refBase, ts); famlk[0].min=famlk[1].min; break;}
	    case 2: { transversions++; famlk[0].SetAlleles(refBase, tvs1);  famlk[0].min=famlk[2].min; break;}
	    case 3: { transversions++; famlk[0].SetAlleles(refBase, tvs2); famlk[0].min=famlk[3].min; break;}
	    case 4: { tstvs1Cnt++; famlk[0].SetAlleles(ts, tvs1); famlk[0].min=famlk[4].min; break; }
	    case 5: { tstvs2Cnt++; famlk[0].SetAlleles(ts, tvs2); famlk[0].min=famlk[5].min; break; }
	    case 6: { tvs1tvs2Cnt++; famlk[0].SetAlleles(tvs1, tvs2); famlk[0].min=famlk[6].min; break;}
	    case -1: { nocall++; break;}
	    default: error("Invalid maxidx!\n");
	    }

	  if(maxidx==-1 || maxidx==0 && par.denovo==false && force_call==false && out_all_sites==false) continue;

	  if(maxidx==0)
	  {
	    if(par.denovo) {
	     double lk_mono = famlk[0].MonomorphismLogLikelihood(refBase);
	     famlk[0].min = 1.0;
	     famlk[0].denovoLR = famlk[0].varllk_noprior[0] - lk_mono;
	     if(famlk[0].denovoLR <= log10(par.denovoLR) && out_all_sites==false && force_call==false) continue;
	    }
	  }
	  else {
	    if(par.denovo)
	    {
	      par.denovo = false;
	      double lk_poly = famlk[0].PolymorphismLogLikelihood(famlk[0].allele1, famlk[0].allele2);
	      famlk[0].denovoLR = famlk[0].varllk_noprior[maxidx] - lk_poly;
	      par.denovo = true;
	    }
	  }

    if(maxidx==0)
    {
      if(par.denovo) {
      famlk[0].denovo_mono = true;
     	famlk[0].CalcPostProb(1.0);
      }
     else { famlk[0].isMono=true; famlk[0].CalcPostProb(1-theta);}
    }
    else {
	famlk[0].isMono=false;
	famlk[0].CalcPostProb(famlk[0].min); //mlefreq is the freq of the reference allele
    }

    denovo ? famlk[0].OutputVCF_denovo(vcfFH) : famlk[0].OutputVCF(vcfFH);
    famlk[0].denovo_mono = false;

    out_cnt++;
    if(force_call && out_cnt >= positionMap.size() ) return(0);
   }

      //output summary statisitcs after the run
      int totalBases = 0;
      for(int i=0; i<5; i++){
	totalBases += refBaseCounts[i];
      }
      printf("Summary of reference -- %s\n", pedGLF.GetNonNULLglf()->label.c_str()); 
      printf("Total Entry Count: %9d\n", totalEntryCnt);
      printf("Total Base Cout: %9d\n", totalBases);
      //printf("Total '0' Base Count: %9d\n", refBaseCounts[0]);
      printf("Non-Polymorphic Count: %9d\n", homoRef);
      printf("Transition Count: %9d\n", transitions);
      printf("Transversion Count: %9d\n", transversions);
      printf("Other Polymorphism Count: %9d\n", tstvs1Cnt+tstvs2Cnt+tvs1tvs2Cnt);
      printf("Filter counts:\n");
      printf("\tminMapQual %u\n", minMapQualFilter);
      printf("\tminTotalDepth %u\n", minTotalDepthFilter);
      printf("\tmaxTotalDepth %u\n", maxTotalDepthFilter);
      printf("Hard to call: %9d\n", nocall);
      printf("Skipped bases: %u\n", totalEntryCnt-homoRef-transitions-transversions-(tstvs1Cnt+tstvs2Cnt+tvs1tvs2Cnt));
      time_t tfinished = t;
      time(&tfinished);
      time_t duration = tfinished - t;
      printf("Analysis ended on %s\n", ctime(&tfinished));
      printf("Running time is %u seconds\n\n", (unsigned int)(duration));

      fflush(vcfFH);

    }

  fclose(vcfFH);
  return(0);
}
