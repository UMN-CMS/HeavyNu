#include "TFile.h"
#include "TTree.h"
#include "SimpleCardfile.hh"
#include "doLimitElectron.hh"
#include <math.h>


int main(int argc, char* argv[]) {

  if (argc!=2) {
    printf("Usage: [cardfile]\n");
    printf("Cardfile format\n");
    printf(" Required\n");
    printf("  DATA     : <filename>\n");
    printf("  DOSYST   : <yes|no> (default no)\n");
    printf("  SIGNALMC : <filename> (more than one is possible)\n");
    printf("  LUMI     : <lumi in pb^-1>\n");
    printf(" Modes (choose one)\n");
    printf("  POINT    : (list of XSEC points to scan)\n");
    printf("  SCAN     : LINEAR <min> <max> <number of points>\n");
    printf("  SCAN     : LOG <min> <max> <number of points>\n");
    printf(" Optional:\n");
    printf("  CYCLES  : <number of MC cycles>\n");
    return 1;
  }

  LimitPointStruct info;
  SimpleCardfile cardfile(argv[1]);
  TFile* rootOutFile=0;
  TTree* rootOutTree=0;

  if (!cardfile.requireEntry("DATA","Root file containing data histogram")) return 2;
  if (!cardfile.requireEntry("SIGNALMC","Root file containing signal MC histogram")) return 2;
  if (!cardfile.requireEntry("LUMI","Luminosity in pb^-1")) return 2;
  

  std::string fdataFile=cardfile.getItem("DATA");
  printf("Loading data from '%s'\n",fdataFile.c_str());
  TFile fdata(fdataFile.c_str());
  
  info.base.lumi=cardfile.getItemFloat("LUMI");
  int cycles=cardfile.getItemInt("CYCLES",0,500);

  bool doSyst=false;
  if (cardfile.hasEntry("DOSYST")) {
    doSyst=cardfile.getItem("DOSYST")=="yes";
  }
  if (doSyst) printf("Running systematic uncertainties\n");
  else printf("Not running systematic uncertainties\n");

  if (cardfile.hasEntry("TUPLEFILE")) {
    rootOutFile=new TFile(cardfile.getItem("TUPLEFILE").c_str(),"RECREATE");
    rootOutTree=new TTree("WRLimit","WR-HeavyNu Limits");
    rootOutTree->Branch("info",&(info.base),"lumi/F:mwr:mnu:xsec:signal:background:data");
    rootOutTree->Branch("cl_sb",&(info.cl_sb),"obs/F:exp:exp_p1s:exp_m1s:exp_p2s:exp_m2s");
    rootOutTree->Branch("cl_b",&(info.cl_b),"obs/F:exp:exp_p1s:exp_m1s:exp_p2s:exp_m2s");
    rootOutTree->Branch("cls",&(info.cls),"obs/F:exp:exp_p1s:exp_m1s:exp_p2s:exp_m2s");
  }

  std::vector<std::string> sigFiles;
  for (int il=0; il<cardfile.getLineCount("SIGNALMC"); il++) 
    for (int it=0; it<cardfile.getItemCount("SIGNALMC",il); it++) 
      sigFiles.push_back(cardfile.getItem("SIGNALMC",il,it));

  for (std::vector<std::string>::const_iterator ifile=sigFiles.begin(); ifile!=sigFiles.end(); ifile++) {

    std::string sigFile=*ifile;
    TFile fsignal(sigFile.c_str());
		
    int iw=0,inu=0;
    sscanf(strstr(sigFile.c_str(),"WR"),"WR%d_nuRmu%d",&iw,&inu);
    
    info.base.mwr=iw;
    info.base.mnu=inu;
    
    std::vector<float> points;
    
    if (cardfile.hasEntry("POINT")) {
      for (int il=0; il<cardfile.getLineCount("POINT"); il++) {
	for (int it=0; it<cardfile.getItemCount("POINT",il); it++) {
	  points.push_back(cardfile.getItemFloat("POINT",il,it,0));
	}
      }
    } else if (cardfile.hasEntry("SCAN")) {
      if (cardfile.getItem("SCAN",0)=="LINEAR"  && cardfile.getItemCount("SCAN",0)==4) {
	double min=cardfile.getItemFloat("SCAN",1);
	double max=cardfile.getItemFloat("SCAN",2);
	int np=cardfile.getItemInt("SCAN",3);
	for (int i=0; i<np; i++) 
	  points.push_back(min+(max-min)/(np-1)*i);    
      } else if (cardfile.getItem("SCAN",0)=="LOG" && cardfile.getItemCount("SCAN",0)==4) {
	double min=cardfile.getItemFloat("SCAN",1);
	double max=cardfile.getItemFloat("SCAN",2);
	int np=cardfile.getItemInt("SCAN",3);
	double factor=pow(max/min,1.0/(np-1));
	for (int i=0; i<np; i++) 
	  points.push_back(min*pow(factor,i));
      }
    }
    
    int ipoint=1;
    for (std::vector<float>::const_iterator ip=points.begin(); ip!=points.end(); ip++) {
      info.base.xsec=*ip;
      //      printf("%d xsec=%f\n",ipoint,info.xsec);
      doLimitElectron(&fdata,&fsignal,cycles,info,doSyst);
      ipoint++;

      if (rootOutTree!=0) rootOutTree->Fill();
    }
  } 
  if (rootOutFile!=0) {
    rootOutFile->Write();
    rootOutFile->Close();
  }
  return 0;
}
