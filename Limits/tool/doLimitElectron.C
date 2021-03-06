#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooArgSet.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooExtendPdf.h"
#include "RooConstVar.h"
#include "RooHistPdf.h"
#include "RooGlobalFunc.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/HypoTestPlot.h"
#include "RooStats/HypoTestResult.h"
#include "TCanvas.h"
#include "TFile.h"
#include "doLimitSetting.hh"
#include <algorithm>

using namespace RooFit;
using namespace RooStats;

void doLimitElectron(TFile* dataf, TFile* signal, int ntoys, LimitPointStruct& info, bool doSyst) {
  double lumi=info.base.lumi; // pb^-1
  char temp[1024];
  const char* finalVarHistData="Data";
  const char* finalVarHistSig="Signal";

  const double sig_scale_factor=info.base.xsec;
  RooWorkspace w("eejj");

  w.factory("mwr[520,2000]"); // observable is called mwr
  w.defineSet("obs","mwr");

  TH1* signalShape=(TH1*)(signal->Get(finalVarHistSig)->Clone("ss"));
  //  TH1* signalNorm=(TH1*)(signal->Get("hNu/njet"));

  //  signalShape->Print();
  //  double signal_eff=signalShape->Integral()*1.0/signalNorm->GetEntries();

  double signal_eff=signalShape->Integral()/34.0;
  //  signalShape->Scale(1.0/(34.0));

  TH1* dataHist=(TH1*)(dataf->Get(finalVarHistData)->Clone("dataHist"));

  info.base.data=dataHist->Integral(int(dataHist->FindBin(520)),int(dataHist->FindBin(2001)));
  RooDataHist* data=new RooDataHist("data","data",RooArgList(*w.var("mwr")),dataHist) ;

  //  data->Print();

  RooDataHist* signals=new RooDataHist("signals","signals",RooArgList(*w.var("mwr")),signalShape) ;
  RooHistPdf* signalPdf=new RooHistPdf("sig_pdf","sig_pdf",RooArgList(*w.var("mwr")),*signals);

  // Gaussian signal
  w.import(*signalPdf);
  //  w.import(*datardh);
		       //  w.factory("Gaussian::sig_pdf(mwr,sig_mean[600],sig_sigma[10])");


    // Constant                  =     1.2754          ±     0.13397
  // Slope                     =     -0.005042       ±     0.00013223

  w.factory("Exponential::bkg_all(mwr,-0.005042)");
  const double bkg_yield_total=1.64;

  // total model with signal and background yields as parameters ( *beware when changing luminosity* )
  w.factory("SUM::main_pdf(sig_yield[20,0,3000]*sig_pdf,bkg_yield[0.1,1e-5,10]*bkg_all)");

  if (doSyst) {
    // note the "+1" required for the defintion of K! ( K - 1   = [fractional error] )
    //    sprintf(temp,"Lognormal::ttb_xsec(ttb_estimate[%f],ttb_yield,%f)",ttb_yield_mc,ttb_yield_ferror+1);
    //  w.factory(temp);
    // sprintf(temp,"Lognormal::zj_xsec(zj_estimate[%f],zj_yield,%f)",zj_yield_mc,zj_yield_ferror+1);
    //w.factory(temp);
  }

  w.var("bkg_yield")->setVal(bkg_yield_total);
  //  w.var("zj_yield")->setVal(zj_yield_mc);
  //  w.var("qcd_yield")->setVal(qcd_yield_data);
  double sig_yield_expected=signal_eff*sig_scale_factor*lumi;
  //  printf("Components: %f %f %f \n",signal_eff,sig_scale_factor,lumi);
  info.base.background=bkg_yield_total;
  //  w.var("bkg_scale")->setVal(info.background);

  if (!doSyst) {
    w.var("bkg_yield")->setConstant();
    //    w.var("ttb_yield")->setConstant();
    //    w.var("zj_yield")->setConstant();
  }
  //  w.var("qcd_yield")->setConstant();
  //  w.var("bkg_scale")->setConstant();

  info.base.signal=sig_yield_expected;
  // The model for the control sample that constrains the background.
  //w.factory("Gaussian::control_pdf(control_meas[50],bkg_yield,10.)");
  
  if (doSyst) {
    w.factory("PROD::control_pdf(ttb_xsec,zj_xsec)");
  //sprintf(temp,"Gaussian::control_pdf(ttb_estimate[%f],ttb_yield,%f)",ttb_yield_mc,ttb_yield_mc*ttb_yield_ferror);
  //w.factory(temp);
  // The total model including the main measurement and the control sample

    w.factory("PROD::main_with_control(main_pdf,control_pdf)");
  }

  // choose which pdf you want to use
  RooAbsPdf* pdfToUse=w.pdf("main_pdf"); // only use main measurement

  if (doSyst) pdfToUse=w.pdf("main_with_control");
    
  // define sets for reference later
  w.defineSet("poi","sig_yield");
  if (doSyst) {
    w.defineSet("nuis","ttb_yield,zj_yield");
  }
  //  RooDataSet* data = w.pdf("main_pdf")->generate(*w.set("obs"),RooFit::Extended());
  //  data->Print();

  // D E F I N E  N U L L  &  A L T E R N A T I V E   H Y P O T H E S E S  
  ModelConfig b_model("B_model", &w);
  b_model.SetPdf(*pdfToUse);
  b_model.SetObservables(*w.set("obs"));
  b_model.SetParametersOfInterest(*w.set("poi"));
  if (doSyst)
    b_model.SetNuisanceParameters(*w.set("nuis")); 
  w.var("sig_yield")->setVal(0.0);
  b_model.SetSnapshot(*w.set("poi"));  

  ModelConfig sb_model("S+B_model", &w);
  sb_model.SetPdf(*pdfToUse);
  sb_model.SetObservables(*w.set("obs"));
  sb_model.SetParametersOfInterest(*w.set("poi"));
  if (doSyst)
    sb_model.SetNuisanceParameters(*w.set("nuis")); 
  w.var("sig_yield")->setVal(sig_yield_expected);
  sb_model.SetSnapshot(*w.set("poi"));

  
  SimpleLikelihoodRatioTestStat slrts(*b_model.GetPdf(),*sb_model.GetPdf());
  slrts.SetNullParameters(*b_model.GetSnapshot());
  slrts.SetAltParameters(*sb_model.GetSnapshot());

  ToyMCSampler toymcsampler(slrts, ntoys);
  toymcsampler.SetGenerateBinned(true);

  HybridCalculator myH2(*data,sb_model, b_model, &toymcsampler);

  if (doSyst) {
    myH2.ForcePriorNuisanceNull(*w.pdf("control_pdf")); 
    myH2.ForcePriorNuisanceAlt(*w.pdf("control_pdf")); 
  }

  HypoTestResult *res = myH2.GetHypoTest();
  //  res->Print(0);

  // here we get the expected limit
  SamplingDistribution* sd_sb=res->GetAltDistribution();
  SamplingDistribution* sd_b=res->GetNullDistribution();

  //  printf("Sizes: %d %d %f\n",sd_b->GetSize(),sd_sb->GetSize());

  std::vector<double> sampl_sb=sd_sb->GetSamplingDistribution();
  std::vector<double> sampl_b=sd_b->GetSamplingDistribution();

  std::sort(sampl_sb.begin(),sampl_sb.end(),std::less<double>());
  std::sort(sampl_b.begin(),sampl_b.end(),std::less<double>());
  int lim=int(sampl_b.size());

  //  double median=sampl_b[sampl_b.size()/2];
  double clmed=-1,clmed_m2s=-1,clmed_m1s=-1,clmed_p1s=-1,clmed_p2s=-1;

  int m1s_pos = int(lim*(1.0-0.3173/2.0)+0.5)+1;
  int p1s_pos = int(lim*(    0.3173/2.0)+0.5)+1;
  int m2s_pos = int(lim*(1.0-0.0455/2.0)+0.5)+1;
  int p2s_pos = int(lim*(    0.0455/2.0)+0.5)+1;

  float val=sampl_b[lim/2];
  for (int i=0; i<lim; i++) {
    if (sampl_sb[i]>=val && clmed<0) {
      clmed=std::max(0.0f,float(i)/(lim/2));
      //      printf("CLMED %d %f %f %f\n",i,clmed,sampl_sb[i],val);
    }
    if (sampl_sb[i]>=sampl_b[m1s_pos] && clmed_m1s<0) clmed_m1s=std::max(0.0f,float(i)/(lim/2));
    if (sampl_sb[i]>=sampl_b[m2s_pos] && clmed_m2s<0) clmed_m2s=std::max(0.0f,float(i)/(lim/2));
    if (sampl_sb[i]>=sampl_b[p1s_pos] && clmed_p1s<0) clmed_p1s=std::max(0.0f,float(i)/(lim/2));
    if (sampl_sb[i]>=sampl_b[p2s_pos] && clmed_p2s<0) clmed_p2s=std::max(0.0f,float(i)/(lim/2));
  }
  //  printf("%f %f %f %f %f %f \n",val,sampl_b[m2s_pos],sampl_b[m1s_pos],sampl_b[p1s_pos],sampl_b[p2s_pos]);
  if (clmed==-1.0 || clmed>1) clmed=1.0;
  if (clmed_m1s==-1.0 || clmed_m1s > 1) clmed_m1s=1.0;
  if (clmed_m2s==-1.0 || clmed_m2s > 1) clmed_m2s=1.0;
  if (clmed_p1s==-1.0 || clmed_p1s > 1) clmed_p1s=1.0;
  if (clmed_p2s==-1.0 || clmed_p2s > 1) clmed_p2s=1.0;
  
  if (clmed>1.0 && clmed<1.001) clmed=1.0;
  if (clmed<0.0) clmed=0.0;

  info.cls.obs=res->CLs();
  info.cls.exp=clmed;
  info.cls.exp_m2s=clmed_m2s;
  info.cls.exp_m1s=clmed_m1s;
  info.cls.exp_p1s=clmed_p1s;
  info.cls.exp_p2s=clmed_p2s;

  
  // for (int i=0; i<20; i++) 
  //  printf("%3d %f %f \n",i,sampl_sb[i],sampl_b[i]);

  double xobs=res->GetTestStatisticData();
  xobs+=0.0001; // important for handling the "exactly zero data" case, irrelevant for all other purposes
  double alt_obs_sb=sd_sb->Integral(-10000,xobs);
  double alt_obs_b=sd_b->Integral(xobs,10000);

  double mean_sb=sd_sb->InverseCDFInterpolate(0.5);
  double mean_b=sd_b->InverseCDFInterpolate(0.5);

  double expected_sb=sd_sb->Integral(-10000,mean_b);
  double expected_b=sd_b->Integral(mean_sb,10000);
  
  //  printf("DETAIL : mean_sb=%f  mean_b=%f\n", mean_sb,mean_b);
  printf("SUMMARY  : LUMI=%.1f ipb MW=%.0f GeV  MNu=%.0f GeV  XSEC=%.2e pb \n",info.base.lumi,info.base.mwr,info.base.mnu,sig_scale_factor);
  //  printf("SUMMARY2 : EXP(s+b)=%.3f OBS(s+b)=%.3f EXP(b)=%.3f OBS(b)=%.3f\n",expected_sb,res->CLsplusb(),expected_b,1.0-res->CLb());
  printf("SUMMARY2 : EXP(s+b)=%.3f OBS(s+b)=%.3f EXP(b)=%.3f OBS(b)=%.3f\n",expected_sb,alt_obs_sb,expected_b,alt_obs_b);
  //  info.cl_sb_obs = res->CLsplusb();
  info.cl_sb.obs = alt_obs_sb;
  info.cl_sb.exp = expected_sb;
  //  info.cl_b_obs = 1.0-res->CLb();
  info.cl_b.obs = alt_obs_b;
  info.cl_b.exp = expected_b;

  info.cl_sb.exp_p1s = sd_sb->Integral(-10000,sd_b->InverseCDFInterpolate(0.5+0.3413));
  info.cl_sb.exp_m1s = sd_sb->Integral(-10000,sd_b->InverseCDFInterpolate(0.5-0.3413)+0.0001);
  info.cl_sb.exp_p2s = sd_sb->Integral(-10000,sd_b->InverseCDFInterpolate(0.5+0.4773));
  info.cl_sb.exp_m2s = sd_sb->Integral(-10000,sd_b->InverseCDFInterpolate(0.5-0.4773)+0.0001);

  printf("SUMMARY3 : EXP(s+b) (-2s,-1s,1s,2s) %.3f %.3f %.3f %.3f\n",info.cl_sb.exp_m2s,info.cl_sb.exp_m1s,info.cl_sb.exp_p1s,info.cl_sb.exp_p2s);

  info.cl_b.exp_p1s = sd_b->Integral(sd_sb->InverseCDFInterpolate(0.5+0.3413),10000);
  info.cl_b.exp_m1s = sd_b->Integral(sd_sb->InverseCDFInterpolate(0.5-0.3413),10000);
  info.cl_b.exp_p2s = sd_b->Integral(sd_sb->InverseCDFInterpolate(0.5+0.4773),10000);
  info.cl_b.exp_m2s = sd_b->Integral(sd_sb->InverseCDFInterpolate(0.5-0.4773),10000);
  printf("SUMMARY3 : EXP(b) (-2s,-1s,1s,2s) %.3f %.3f %.3f %.3f\n",info.cl_b.exp_m2s,info.cl_b.exp_m1s,info.cl_b.exp_p1s,info.cl_b.exp_p2s);
  printf("SUMMARY3 : CLS %f | %f %f %f %f %f\n",res->CLs(),clmed,clmed_p2s,clmed_p1s,clmed_m1s,clmed_m2s);


  printf("SUMMARY4: %.2f background / %.2f (%.2f%%) signal / %d data\n",info.base.background,info.base.signal,signal_eff*100, int(info.base.data));


  
#ifndef MAKE_MAIN
  TCanvas* c2=new TCanvas("c2","c2",900,900);
  c2->Divide(1,2);

  c2->cd(1);
  // and make a plot. number of bins is optional (default: 100)
  HypoTestPlot *plot = new HypoTestPlot(*res, 80); 
  plot->Draw();

  c2->cd(2);
  RooPlot* mesframe = w.var("mwr")->frame() ;

  data->plotOn(mesframe);
  //  w.pdf("bkg_pdf")->plotOn(mesframe);
  //w.pdf("main_pdf")->plotOn(mesframe);
  w.pdf("main_pdf")->plotOn(mesframe,Normalization(sig_yield_expected+info.base.background));
  w.pdf("main_pdf")->plotOn(mesframe,Normalization(sig_yield_expected+info.base.background),Components("bkg_all"),LineStyle(kDashed));

  mesframe->Draw();

#endif

  delete sd_sb;
  delete sd_b;
  delete res;

  delete dataHist;
  delete signalShape;
  delete signals;
  delete data;

}

void testLimitElectron(TFile* f0, TFile* f1, double xsec=1.0, bool doSyst=true) {
  LimitPointStruct lps;
  lps.base.xsec=xsec;
  lps.base.lumi=34.0;

  doLimitElectron(f0,f1,2000,lps,doSyst);
}
