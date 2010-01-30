#include <iostream>
#include <string>

#include <map>
#include <vector>
#include <functional>

#include "ECALTime/EcalTimePi0/interface/EcalTimePi0TreeContent.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TFile.h"

#define BarrelLimit 1.479
#define EndcapLimit 3.0

// initial authors P. Govoni et al
// authors: S. Cooper and G. Franzoni (UMN)

//! main program
int main (int argc, char** argv)
{

  // default output file
  std::string outputRootName = "outputHistos.root";
  std::string stringGenericOption    = "--";
  std::string stringHelp             = "--help";
  std::string stringInputFileName    = "--i";
  std::string stringOutFileName      = "--o";
  std::string stringETGammaMinEB     = "--eTGammaMinEB";
  std::string strings4s9GammaMinEB   = "--s4s9GammaMinEB";
  std::string stringeTPi0MinEB       = "--stringeTPi0MinEB";
  std::string stringETGammaMinEE     = "--eTGammaMinEE";
  std::string strings4s9GammaMinEE   = "--s4s9GammaMinEE";
  std::string stringeTPi0MinEE       = "--stringeTPi0MinEE";
  std::string stringNumEvents        = "--n";

  std::vector<std::string> listOfFiles;
  int   numEvents    =-1;
  float	eTGammaMinEB   = 0.2;
  float s4s9GammaMinEB = 0.85;
  float eTPi0MinEB     = 0.65;
  float	eTGammaMinEE   = 0.250;
  float s4s9GammaMinEE = 0.85;
  float eTPi0MinEE     = 0.800;

  //gf: support development
  //std::cout << "\nargc:       " << argc << std::endl;
  //for (int v=0; v<argc; v++ ){      std::cout << "argument: " << v << " argv: " << argv[v] << std::endl;    }

  
  // if no arguments are passed, suggest help
  if (argc < 2){
    std::cerr << "\n\tERROR: specify arguments, e.g. --help\n" << std::endl ;
    exit (1) ;  
  }

  // loop over input options
  for (int v=1; v<argc; v++ )
    {
      //std::cout << "argv number " << v << " is: " << argv[v] << std::endl;
      
      if (argv[v] == stringHelp) { // help message
	std::cout << " --help : display help" << std::endl ;
	std::cout << " --o : set name of output root file name (e.g. histograms.root)" << std::endl ;
	std::cout << " --n : number of events" << std::endl ;
	std::cout << " --eTGammaMinEB: min eT for EB gammas" << std::endl;
	std::cout << " --s4s9GammaMinEB: min EB shower shape" << std::endl;
	std::cout << " --eTPi0MinEB min eT for EB pi0 candidate" << std::endl;
	std::cout << " --eTGammaMinEE: min eT for EE gammas" << std::endl;
	std::cout << " --s4s9GammaMinEE: min EE shower shape" << std::endl;
	std::cout << " --eTPi0MinEE min eT for EE pi0 candidate" << std::endl;
	std::cout << " --i <list of strings> list of input files" << std::endl ;     
	exit(1);      }

      
      else if (argv[v] == stringNumEvents) { // set number of events
	std::cout << "events number" << std::endl;
	numEvents=atoi(argv[v+1]);
	v++;
      }

      
      else if (argv[v] == stringETGammaMinEB) { // choose et cut for EB single cluster
	eTGammaMinEB = atof(argv[v+1]);
	v++;
      }
      
      else if (argv[v] == strings4s9GammaMinEB) { // choose cut for EB shower shape
	s4s9GammaMinEB = atof(argv[v+1]);
	v++;
      }
      
      else if (argv[v] == stringeTPi0MinEB) { // choose et cut for EB pi0 candidate
	eTPi0MinEB = atof(argv[v+1]);
	v++;
      }

      
      else if (argv[v] == stringETGammaMinEE) { // choose et cut for EE single cluster
	eTGammaMinEE = atof(argv[v+1]);
	v++;
      }
      
      else if (argv[v] == strings4s9GammaMinEE) { // choose cut for EE shower shape
	s4s9GammaMinEE = atof(argv[v+1]);
	v++;
      }
      
      else if (argv[v] == stringeTPi0MinEE) { // choose et cut for EE pi0 candidate
	eTPi0MinEE = atof(argv[v+1]);
	v++;
      }

      else if (argv[v] == stringOutFileName) { // set output file
	outputRootName = argv[v+1];
	v++;
      }

      // handle here the case of multiple arguments for input files
      else if (argv[v] == stringInputFileName){// && v<(argc-1) ) {

	for (int u=v+1; u<argc; u++) {
	  
	  if ( 0==std::string(argv[u]).find( stringGenericOption ) ){
	    if ( 0==listOfFiles.size())  {std::cout << "no input files listed" << std::cout;}
	    //else  {std::cout << "no more files listed, found: " << argv[u] << std::cout;}
	    break;
	  }

	  else {  listOfFiles.push_back(argv[u]);
	    v++;
	  }

	}// loop on arguments following --i

	continue;

      }//end 'if input files'

      
      else
	{std::cout << "input format unrecognized" << std::endl; exit(1);}

    }// loop over arguments input to the program


  
  if (listOfFiles.size()==0){
    std::cout << "\tno input file found" << std::endl;
    return(1);
  }
  else{
    std::cout << "\tfound " << listOfFiles.size() << " input files: " << std::endl;
    for(std::vector<std::string>::const_iterator  file_itr=listOfFiles.begin(); file_itr!=listOfFiles.end(); file_itr++){
      std::cout << "\t" << (*file_itr) << std::endl;
    }
  }
  




  // Tree construction
  TChain * chain = new TChain ("EcalTimePi0Analysis") ;
  std::vector<std::string>::const_iterator file_itr;
  for(file_itr=listOfFiles.begin(); file_itr!=listOfFiles.end(); file_itr++){
    chain->Add( (*file_itr).c_str() );
  }
  int nEntries = chain->GetEntries () ;
  if (numEvents==-1) numEvents = nEntries;
  std::cout << "\n\tFOUND "         << nEntries << " events" << std::endl ;    
  std::cout << "\tWILL run on: "    <<  numEvents << " events" << std::endl;
  std::cout << "\tOutput file: "    <<  outputRootName << std::endl;
  std::cout << "\teTGammaMinEB: "   <<  eTGammaMinEB << std::endl;
  std::cout << "\ts4s9GammaMinEB: " <<  s4s9GammaMinEB << std::endl;
  std::cout << "\teTPi0MinEB: "     <<  eTPi0MinEB << std::endl;
  std::cout << "\teTGammaMinEE: "   <<  eTGammaMinEE << std::endl;
  std::cout << "\ts4s9GammaMinEE: " <<  s4s9GammaMinEE << std::endl;
  std::cout << "\teTPi0MinEE: "     <<  eTPi0MinEE << std::endl;
	
  EcalTimePi0TreeContent treeVars ; 
  setBranchAddresses (chain, treeVars) ;

  // Initialize output root file
  TFile saving (outputRootName.c_str (),"recreate") ;
  saving.cd () ;

  // Initialize histograms -- xtals
  TH1F* xtalEnergyHist_ = new TH1F("XtalEnergy","Crystal energy;GeV",110,-1,10);
  TH1F* xtalTimeHist_ = new TH1F("XtalTime","Time of all crystals;ns",150,-75,75);
  TH1F* xtalIEtaHist_ = new TH1F("xtalIEta","i#eta of crystal",171,-85,86);
  TH1F* xtalIPhiHist_ = new TH1F("xtalIPhi","i#phi of crystal",361,1,361);
  TH1F* xtalIXHist_ = new TH1F("xtalIX","ix of crystal",101,1,101);
  TH1F* xtalIYHist_ = new TH1F("xtalIY","iy of crystal",101,1,101);
  // TH1F* xtalStatusHist_ = new TH1F("XtalStatus","Crystal status flag",16,0,15);
  // TH2F* xtalOccupancyHistEB_ = new TH2F("XtalOccupancyEB","Crystal occupancy;i#phi;i#eta",360,1.,361.,172,-86,86);

  // Initialize histograms -- BasicClusters
  TH1F* BCNumPerEventHist_ = new TH1F("BCNumPerEvent","Number of BC per event",100,0,100);
  TH1F* BCNumCrysHist_ = new TH1F("BCNumCrys","Number of crystals per BC",10,0,10);
  TH1F* BCEnergyHist_ = new TH1F("BCEnergy","Energy of BCs;GeV",100,0,25);
  TH1F* BCEtHist_ = new TH1F("BCEt","E_{T} of BCs;GeV",100,0,25);
  TH2F* BCOccupancyEBHist_  = new TH2F("BCOccupancyEB","BC occupancy;i#eta;i#phi",171,-85,86,361,1.,361.);
  TH2F* BCOccupancyEEPHist_  = new TH2F("BCOccupancyEEP","BC occupancy;ix;iy",101,1.,101.,101,1,101);
  TH2F* BCOccupancyEEMHist_  = new TH2F("BCOccupancyEEM","BC occupancy;ix;iy",101,1.,101.,101,1,101);
  TH2F* BCOccupancyHistAny_ = new TH2F("BCOccupancyAny","BC occupancy;#eta;#phi",50,-3.5,3.5,50,-1*TMath::Pi(),TMath::Pi());
  TH1F* BCEtaHist_ = new TH1F("Cluster #eta","#eta of cluster",171,-3.5,3.5);
  TH1F* BCPhiHist_ = new TH1F("Cluster #phi","#phi of cluster",50,-1*TMath::Pi(),TMath::Pi());
  TH1F* BCClusterShapeEEPHist_ = new TH1F("EEP cluster shape","e2x2 / e3x3",65,-0.1,1.2);
  TH1F* BCClusterShapeEEMHist_ = new TH1F("EEM cluster shape","e2x2 / e3x3",65,-0.1,1.2);
  TH1F* BCClusterShapeEBHist_  = new TH1F("EB cluster shape","e2x2 / e3x3",65,-0.1,1.2);

  // Initialize histograms -- diphotons
  TH1F* massDiGammaHist_ = new TH1F("massDiGamma","m(#gamma#gamma)",50,0,0.500);
  TH1F* massDiGammaEBHist_ = new TH1F("massDiGamma EB","m(#gamma#gamma) EB",50,0,0.500);
  TH1F* massDiGammaEEPHist_ = new TH1F("massDiGamma EEP","m(#gamma#gamma) EEP",50,0,0.500);
  TH1F* massDiGammaEEMHist_ = new TH1F("massDiGamma EEM","m(#gamma#gamma) EEM",50,0,0.500);


  //loop over entries
  for (int entry = 0 ; (entry < nEntries && entry < numEvents); ++entry)
    {
      chain->GetEntry (entry) ;

      bool speak=false;
      if (entry<10 || entry%100==0) speak=true;

      if (speak) std::cout << "------> reading entry " << entry << " <------\n" ; 


      // loop on calorimetric quantities

      if (speak)  std::cout << "  found " << treeVars.nSuperClusters << " superclusters" << std::endl ;
      if (speak)  std::cout << "  found " << treeVars.nClusters << " basic clusters" << std::endl ;
      BCNumPerEventHist_->Fill(treeVars.nClusters);


        /////////////////////////////////////////////////////
	//loop on basic cluster

	for (int bCluster=0; bCluster < treeVars.nClusters; bCluster++)
	  {

	    float eBC=0; // calculate energy of BC for validation
	    for (int cryInBC=0; cryInBC < treeVars.nXtalsInCluster[bCluster]; cryInBC++){
	      eBC+= treeVars.xtalInBCEnergy[bCluster][cryInBC];}

            BCEnergyHist_->Fill(treeVars.clusterEnergy[bCluster]);
            BCEtHist_->Fill(treeVars.clusterTransverseEnergy[bCluster]);
            BCNumCrysHist_->Fill(treeVars.nXtalsInCluster[bCluster]);

	    // basic cluster occupancy in physics coordinates
	    BCOccupancyHistAny_ -> Fill(treeVars.clusterEta[bCluster],treeVars.clusterPhi[bCluster]);
	    BCEtaHist_ -> Fill(treeVars.clusterEta[bCluster]);
	    BCPhiHist_ -> Fill(treeVars.clusterPhi[bCluster]);

	    //  basic cluster occupancy in detector coordinates, using first cry of BC as a representative
            if(treeVars.xtalInBCIEta[bCluster][0] != -999999)                                        // ieta=-999999 tags EE
              BCOccupancyEBHist_->Fill(treeVars.xtalInBCIEta[bCluster][0],treeVars.xtalInBCIPhi[bCluster][0]);

            else if (treeVars.xtalInBCIx[bCluster][0] != -999999 && treeVars.clusterEta[bCluster]>0) // ix=-999999 tags EB
              BCOccupancyEEPHist_->Fill(treeVars.xtalInBCIx[bCluster][0],treeVars.xtalInBCIy[bCluster][0]);

            else if (treeVars.xtalInBCIx[bCluster][0] != -999999 && treeVars.clusterEta[bCluster]<0) // ix=-999999 tags EB
              BCOccupancyEEMHist_->Fill(treeVars.xtalInBCIx[bCluster][0],treeVars.xtalInBCIy[bCluster][0]);

	    if (speak)  std::cout << "\tbCluster: num"               << bCluster 
				  << "\t eBC: "                      << treeVars.clusterEnergy[bCluster]
				  << "\t eBC_predicted: "            << eBC
				  << "\n\t et: "                     << treeVars.clusterTransverseEnergy[bCluster]
				  << "\t predicted et: "             << treeVars.clusterEnergy[bCluster]*sin(2*atan(exp(-1* treeVars.clusterEta[bCluster] )) )
				  << " eta: "                        << treeVars.clusterEta[bCluster]
				  << "\n\t num crystals: "           << treeVars.nXtalsInCluster[bCluster]
				  << "\n\t\tfirst crystal:  \tieta " << treeVars.xtalInBCIEta[bCluster][0] 
				  << "\teta "                        << treeVars.xtalInBCEta[bCluster][0] 
				  << " \t energy "                   << treeVars.xtalInBCEnergy[bCluster][0] 
				  << " \t ADC "                      << treeVars.xtalInBCAmplitudeADC[bCluster][0] 
				  << " \t time "                     << treeVars.xtalInBCTime[bCluster][0] 
				  << std::endl;

	    for(int thisCry=0; thisCry<treeVars.nXtalsInCluster[bCluster]; thisCry++)
	      {
		if (treeVars.xtalInBCIEta[bCluster][thisCry]!=-999999)  xtalIEtaHist_ -> Fill (treeVars.xtalInBCIEta[bCluster][thisCry]);
		if (treeVars.xtalInBCIPhi[bCluster][thisCry]!=-999999)  xtalIPhiHist_ -> Fill (treeVars.xtalInBCIPhi[bCluster][thisCry]);
		if (treeVars.xtalInBCIx[bCluster][thisCry]  !=-999999)  xtalIXHist_   -> Fill (treeVars.xtalInBCIx[bCluster][thisCry]);
		if (treeVars.xtalInBCIy[bCluster][thisCry]  !=-999999)  xtalIYHist_   -> Fill (treeVars.xtalInBCIy[bCluster][thisCry]);
	      }
	    
	  }
	
	if (speak) std::cout << "  found " << treeVars.nXtals << " crystals\n" ;    
      //PG loop over crystals
      for (int XTLindex = 0 ; XTLindex < treeVars.nXtals ; ++XTLindex)
        {
          //TODO: For this to work, we need another way to tell EB from EE,
          //      as hashedIndex itself is not enough!
          //if(EBDetId::validHashIndex(treeVars.xtalHashedIndex[XTLindex]))
          //{
          //  EBDetId ebDet = EBDetId::unhashIndex(treeVars.xtalHashedIndex[XTLindex]);   
          //  xtalOccupancyHistEB_->Fill(ebDet.iphi(),ebDet.ieta());
          //}
          //else if(EEDetId::validHashIndex(treeVars.xtalHashedIndex[XTLindex]))
          //  EEDetId eeDet = EEDetId::unhashIndex(treeVars.xtalHashedIndex[XTLindex]);
          //else
          //  std::cout << "Crystal with hash: " << 

          xtalEnergyHist_->Fill(treeVars.xtalEnergy[XTLindex]);
          xtalTimeHist_->Fill(treeVars.xtalTime[XTLindex]);
          //TODO
          //xtalStatusHist_->Fill

        } //PG loop over crystals




      float eTA, eTB ;
      float e22A, e33A,    e22B, e33B;
      float eTGammaMinA,   eTGammaMinB;
      float s4s9GammaMinA, s4s9GammaMinB;
      bool  AisEB,         BisEB;
      float eTPi0Min;

      for (int bClusterA=0; bClusterA < treeVars.nClusters; bClusterA++)
	{
	  eTA = treeVars.clusterTransverseEnergy[bClusterA];

	  e22A = treeVars.clusterE2x2[bClusterA];
	  e33A = treeVars.clusterE3x3[bClusterA];

	  // discriminate between EE and EB and set thresholds accordingly
	  if ( fabs(treeVars.clusterEta[bClusterA]) < BarrelLimit) {
	    AisEB         = true;
	    eTGammaMinA   = eTGammaMinEB;
	    s4s9GammaMinA = s4s9GammaMinEB;
	  }
	  else{
	    AisEB         = false;
	    eTGammaMinA   = eTGammaMinEE;
	    s4s9GammaMinA = s4s9GammaMinEB;
	  }


	  if(treeVars.clusterEta[bClusterA]<-1.4)     BCClusterShapeEEMHist_ -> Fill(e22A/e33A);
	  else if(treeVars.clusterEta[bClusterA]>1.4) BCClusterShapeEEPHist_ -> Fill(e22A/e33A);
	  else	                                      BCClusterShapeEBHist_  -> Fill(e22A/e33A);


	  // first selecton cut: photon candidate Et
	  if( eTA < eTGammaMinA ) continue;
	  
	  // second selection cut: cluster shape
	  if ( e22A/e33A < s4s9GammaMinA ) continue;
	  
	  for (int bClusterB=(bClusterA+1); bClusterB < treeVars.nClusters; bClusterB++)
	    {

	      eTB = treeVars.clusterTransverseEnergy[bClusterB];
	      
	      e22B = treeVars.clusterE2x2[bClusterB];
	      e33B = treeVars.clusterE3x3[bClusterB];
	      
	      // discriminate between EE and EB and set thresholds accordingly
	      if ( fabs(treeVars.clusterEta[bClusterB]) < BarrelLimit) {
		BisEB         = true;
		eTGammaMinB   = eTGammaMinEB;
		s4s9GammaMinB = s4s9GammaMinEB;
	      }
	      else{
		BisEB         = false;
		eTGammaMinB   = eTGammaMinEE;
		s4s9GammaMinB = s4s9GammaMinEB;
	      }
	      

	      // first selecton cut: photon candidate Et
	      if( eTB < eTGammaMinB ) continue;

	      // second selection cut: cluster shape
	      if ( e22B/e33B < s4s9GammaMinB ) continue;
	      
	      math::PtEtaPhiMLorentzVectorD gammaA (eTA, treeVars.clusterEta[bClusterA], treeVars.clusterPhi[bClusterA], 0);
	      math::PtEtaPhiMLorentzVectorD gammaB (eTB, treeVars.clusterEta[bClusterB], treeVars.clusterPhi[bClusterB], 0);
	      
	      math::PtEtaPhiMLorentzVectorD pi0Candidate = gammaA + gammaB;
	      
	      // std::cout << "gammaA: " << gammaA << " " << gammaA.M() << "\t\t gammaB: " << gammaB << " " << gammaB.M() << std::endl;
	      // std::cout << "pi0Candidate: " << pi0Candidate << " " << pi0Candidate.M() << std::endl;


	      if ( fabs(pi0Candidate.Eta()) < BarrelLimit) {
		eTPi0Min = eTPi0MinEB;	      }
	      else{		eTPi0Min = eTPi0MinEE;
	      }


	      // third selection cut: pi0 candidate Et
	      if(pi0Candidate.Et() < eTPi0Min ) continue;
	      
	      massDiGammaHist_ -> Fill(pi0Candidate.M());
	      
	      if(treeVars.clusterEta[bClusterA]<-1.4)     massDiGammaEEMHist_ -> Fill(pi0Candidate.M());
	      else if(treeVars.clusterEta[bClusterA]>1.4) massDiGammaEEPHist_ -> Fill(pi0Candidate.M());
	      else	                                  massDiGammaEBHist_  -> Fill(pi0Candidate.M());
	      
	      
	    }//loop on candidateB
	}//loop on candidateA
      
      
      

    } //PG loop over entries


  BCNumPerEventHist_->Write();
  BCEnergyHist_->Write();
  BCEtHist_->Write();
  BCNumCrysHist_->Write();
  BCOccupancyEBHist_->Write();
  BCOccupancyEEPHist_->Write();
  BCOccupancyEEMHist_->Write();
  BCOccupancyHistAny_->Write();
  BCEtaHist_->Write();
  BCPhiHist_->Write();
  BCClusterShapeEEPHist_->Write();
  BCClusterShapeEEMHist_->Write();
  BCClusterShapeEBHist_->Write();

  xtalEnergyHist_->Write(); 
  xtalTimeHist_->Write();
  xtalIEtaHist_->Write();
  xtalIPhiHist_->Write();
  xtalIXHist_ ->Write();
  xtalIYHist_ ->Write();

  // xtalStatusHist_->Write();
  // xtalOccupancyHistEB_->Write();
  massDiGammaHist_-> Write(); 
  massDiGammaEBHist_-> Write(); 
  massDiGammaEEPHist_-> Write(); 
  massDiGammaEEMHist_-> Write(); 

  saving.Close () ;

  delete chain ;
  //delete BCNumPerEventHist_;
  //delete BCEnergyHist_;
  //delete BCEtHist_;
  //delete BCNumCrysHist_;
  //delete BCOccupancyHistEB_;
  //delete xtalEnergyHist_; 
  //delete xtalTimeHist_;
  //delete xtalStatusHist_;
  //delete xtalOccupancyHistEB_;
  
  return 0 ;
}
