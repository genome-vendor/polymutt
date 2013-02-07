#include "PedVCF.h"
#include "FamilyLikelihoodSeq_VCF.h"

PedVCF::PedVCF(){
 pedGLF = NULL;
 ped = NULL;
 tstv_ratio = 2.0;
 meta_data = "";
}
PedVCF::~PedVCF(){}


double PedVCF::GetPrior_ts(double tstv_ratio)
{  
  return(tstv_ratio/(tstv_ratio + 1));
}

double PedVCF::GetPrior_tv(double tstv_ratio)
{
  return(0.5/(tstv_ratio+1));
}

bool PedVCF::isTs(int a1, int a2)
{
 return ( ( (a1==1 && a2==3) || (a1==2 && a2==4) ) ? true : false );
}

bool PedVCF::isTv(int a1, int a2)
{
 return(!isTs(a1, a2));
}

void PedVCF::SetVCFInput(String& vcfInput){ 
  vcfInFile = std::string(vcfInput.c_str());
  VCFInputFile vin(vcfInFile.c_str());

   std::vector<std::string> pids;
   vin.getVCFHeader()->getPeopleName(&pids);

}
void PedVCF::SetPar(CmdLinePar& cmdpar) { par = cmdpar; }

void PedVCF::VarCallFromVCF()
{
 VCFInputFile vin(vcfInFile.c_str());

 FamilyLikelihoodSeq_VCF famlk_vcf;
 famlk_vcf.SetPedigree(ped);
 famlk_vcf.SetGLF(pedGLF);
 famlk_vcf.SetVinPtr(&vin);
 famlk_vcf.SetCmdLinePar(&par);
 famlk_vcf.SetTheta(par.theta);
 famlk_vcf.InitFamilyLikelihoodES();

if(ped==NULL) error("ped is null\n");

 famlk_vcf.MapPID2Traverse();

 //Set up included PIDs and META data
 std::vector<std::string> pids;
 vin.getVCFHeader()->getPeopleName(&pids);

 famlk_vcf.includedPIDs.clear();
 famlk_vcf.included.clear();

 for(int i=0; i<pids.size(); i++)
 {
	std::map<std::string, std::pair<int, int> >::iterator it;
        it = famlk_vcf.pid2traverse.find(pids[i]);
        if(it!=famlk_vcf.pid2traverse.end() ) {
		famlk_vcf.includedPIDs.push_back(pids[i]);
		famlk_vcf.included[pids[i]]=1;
	}
 }

 vin.includePeople(famlk_vcf.includedPIDs);

 //vin.rewriteVCFHeader();
 //VCFHeader *vh = vin.getVCFHeader();
 //std::string header = (*vh)[vh->size()-1];

  std::string header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  for(int i=0; i<famlk_vcf.includedPIDs.size(); i++) {
	header += "\t";
	header += famlk_vcf.includedPIDs[i];
  }

  meta_data = "##fileformat=VCFv4.1\n";
  meta_data += "##Polymutt=";
  meta_data += famlk_vcf.par->cmd;
  meta_data += "\n";
  meta_data = meta_data +"##Note=VCF file modified by polymutt. Updated fileds include: QUAL, GT and GQ, AF and AC. NOTE: modification was applied only to biallelic variants\n" +
	"##FILTER=<ID=LOWDP,Description=\"Low Depth filter when the average depth per sample is lessn than 1\">\n" +
	"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Read Depth\">\n" +
	"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternative Allele Frequency\">\n" +
	"##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Alternative Allele Count\">\n" +
	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
	"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n" +
	"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n" +
	"##FORMAT=<ID=PL,Number=3,Type=Integer,Description=\"Phred-scaled Genotype Likelihoods\">\n" +
	"##FORMAT=<ID=GL,Number=3,Type=Float,Description=\"Log10 Genotype Likelihoods\">\n" +
	header + "\n";

 double polyPrior = famlk_vcf.GetPolyPrior();
 double polyPrior_indel = famlk_vcf.GetPolyPrior_indel();

 int cnt=0;

 FILE *outVCF = fopen(par.vcfOutFile.c_str(), "w");
 if(outVCF==NULL) error("Open outpuf VCF file %f failed!\n", par.vcfOutFile.c_str() );
 //Note: this is for temp use
 fprintf(outVCF, "%s", meta_data.c_str() );

 double qual, posterior;
 bool flag = false;
 while(vin.readRecord())
 {
  famlk_vcf.FillPenetrance();

  if(flag==false) { printf("Total samples in both VCF and PED files: %d\n\n", famlk_vcf.nSamplesWithData); flag = true; }

  if(famlk_vcf.withdata_cnt==0 || famlk_vcf.isBiallelic==false){ famlk_vcf.OutputVCF(outVCF); continue; }

  if(String((*famlk_vcf.r).getChrom()) == par.chrX_label) famlk_vcf.SetNonAutosomeFlags(true, false, false); 
  else if(String((*famlk_vcf.r).getChrom()) == par.chrY_label) famlk_vcf.SetNonAutosomeFlags(false, true, false); 
  else if(String((*famlk_vcf.r).getChrom()) == par.MT_label) famlk_vcf.SetNonAutosomeFlags(false, false, true); 
  else famlk_vcf.SetNonAutosomeFlags(false, false, false);

  //famlk_vcf.PrintPenetranceMatrix(); printf("--\n");
//  famlk_vcf.PrintLogLKMatrix();

  double mono = famlk_vcf.MonomorphismLogLikelihood(famlk_vcf.ref);

//printf("MIN1=%f\n", famlk_vcf.GetMinimizer() );

  double poly = famlk_vcf.PolymorphismLogLikelihood(famlk_vcf.ref, famlk_vcf.alt);

//printf("MIN2=%f\n", famlk_vcf.GetMinimizer() );

  double llk_alt, llk_ref;

  if(!famlk_vcf.isIndel) {
    llk_alt = log10(polyPrior * isTs(famlk_vcf.ref, famlk_vcf.alt) ? GetPrior_ts(tstv_ratio) : GetPrior_tv(tstv_ratio) ) + poly;
    llk_ref = log10(1-polyPrior) + mono;
  }

  if(famlk_vcf.isIndel){
    llk_alt = log10(polyPrior_indel) + poly;
    llk_ref = log10(1-polyPrior_indel) + mono;
  }

  if(llk_alt-llk_ref>10) qual = 10.0*(llk_alt - llk_ref);
  else {
    posterior = 1/(1+pow(10, llk_ref-llk_alt));
    qual = -10*log10(1-posterior);
  }

  famlk_vcf.CalcPostProb(famlk_vcf.GetMinimizer());
  famlk_vcf.SetQUAL(qual);
  famlk_vcf.OutputVCF(outVCF);

//  printf("POS=%d QUAL=%f poly: %f mono: %f at record %d : %f\n", famlk_vcf.position, qual, poly, mono, ++cnt, famlk_vcf.GetMinimizer());
 }
}
