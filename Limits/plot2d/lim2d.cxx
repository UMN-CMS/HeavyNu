// utility to draw exclusion plot
// to run do:
//   $ cd CMSSW/src/heavynu/datafit
//   $ root -b -q -l 'lim2d.cxx+'

#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TAxis.h"
#include "TList.h"
#include "TGraph.h"
#include "TLine.h"
#include "TBox.h"
#include "TH2.h"
#include "TLatex.h"
#include "TLegend.h"

#include "sampledb.h"
#include "tdrstyle.C"

#include <vector>
#include <algorithm>
using namespace std;

//const double xmax =1500;
const double xmax =1800;
const double ymax =1100;

//const double luminvpb=36;
//const double luminvpb=1000;
const double luminvpb=191;

struct XYZ {
  int x, y;
  float z;
  
  bool operator == (const XYZ& v) const
  {
    return (x == v.x) && (y == v.y);
  }
};

struct LimitsResults {
  TGraph2D* limit_gr;
  TGraph2D* xs_gr;
  TGraph2D* diff_gr;
  TGraph* cont0;
};

static bool sortByXval(std::pair<double,double> u, std::pair<double,double> v)
{
   return u.first < v.first ;
}

LimitsResults calcLimits(const char* fname,
			 const char* scanfmt = "%lg %lg %lg")
{
  // === read limits:
  std::cout << "Checkpoint 1.0" << std::endl ; 
  TGraph2D* limit_gr = new TGraph2D(fname,scanfmt);
  if (limit_gr->IsZombie()) exit(-1);
  limit_gr->SetNameTitle("limit", "limit;WR;nuR;#sigma");

  limit_gr->GetYaxis()->SetRangeUser(0, xmax);
  
  std::cout << "Checkpoint 1.1" << std::endl ; 
  // === read cross-sections:
  // select latest signal samples with ECM = 7 TeV:
  const SampleDB sig7 = SampleDB().find("type", 1).find("channel", 2).find("ECM", 7) ; 
  
  // remove duplicates:
  vector<XYZ> ps;
  for (int i = 0; i < sig7.size(); ++i) {
    const XYZ p = {sig7[i].MW, sig7[i].MNu, sig7[i].CS};
    if (find(ps.begin(), ps.end(), p) == ps.end())
      ps.push_back(p);
  }
  std::cout << "Checkpoint 1.2" << std::endl ; 
  
  // and fill graph:
  TGraph2D* xs_gr = new TGraph2D(ps.size());
  xs_gr->SetNameTitle("xs", "xs;WR;nuR;#sigma");
  
  for (size_t i = 0; i < ps.size(); ++i) {
    xs_gr->SetPoint(i, ps[i].x, ps[i].y, ps[i].z);
  }

  xs_gr->GetYaxis()->SetRangeUser(0, xmax);
  
  cout << "xs: WR = " << xs_gr->GetXmin() << " - " << xs_gr->GetXmax()
       << " NuR = "   << xs_gr->GetYmin() << " - " << xs_gr->GetYmax()
       << endl;
  
  cout << "lim: WR = " << limit_gr->GetXmin() << " - " << limit_gr->GetXmax()
       << " NuR = "    << limit_gr->GetYmin() << " - " << limit_gr->GetYmax()
       << endl;
  
  std::cout << "Checkpoint 1.3" << std::endl ; 
  // === calc difference:
  TGraph2D* diff_gr = new TGraph2D();
  diff_gr->SetNameTitle("diff", "diff;WR;nuR;#sigma (pb)");
  
  const float xmin = max(limit_gr->GetXmin(), xs_gr->GetXmin());
  const float xmax = min(limit_gr->GetXmax(), xs_gr->GetXmax());
  const float ymin = max(limit_gr->GetYmin(), xs_gr->GetYmin());
  const float ymax = min(limit_gr->GetYmax(), xs_gr->GetYmax());
  
  int i = 0;
  for (float wr = xmin; wr < xmax + 1; wr += (xmax - xmin)/40.) {
    for (float nur = ymin; nur < ymax + 1; nur += (ymax - ymin)/40.) {
      double xs = xs_gr->Interpolate(wr,nur);
      double lim = limit_gr->Interpolate(wr,nur);
      if( (xs != 0.0 ) || (lim != 0.0) ) {
	const float diff = xs - lim;
	diff_gr->SetPoint(i, wr, nur, diff);
      }
      i++;
    }
  }

  diff_gr->GetYaxis()->SetRangeUser(0, xmax);

  std::cout << "Checkpoint 1.4" << std::endl ; 
  
  //diff_gr->Draw("surf1");
  //TList* list = el.diff_gr->GetContourList(0);
  //TGraph* cont0 = (TGraph*) list->At(0);
  
  const LimitsResults res = {limit_gr, xs_gr, diff_gr, 0};
  return res;
}                                                          // calcLimits

//======================================================================

TGraph *genContour(const char* fname,
		   const char* scanfmt,
		   const char* tag)
{
  TCanvas* c = new TCanvas(tag, tag, 800, 700);
  c->Divide(2, 2);
  c->SetLogy(0) ; 

  std::cout << "Checkpoint 1" << std::endl ; 
  
  //const LimitsResults el = calcLimits("Data_CS_Up_Lim_E_All.txt");
  const LimitsResults el = calcLimits(fname,scanfmt);
  
  std::cout << "Checkpoint 2" << std::endl ; 
  // =========== plot limits:
  c->cd(1);
  gPad->SetLogy(0);
  gPad->SetLogz(0);

  el.limit_gr->Draw("colz");
  
  // =========== plot crossection:
  c->cd(2);
  gPad->SetLogy(0);
  gPad->SetLogz(1);

  el.xs_gr->DrawClone("colz");
  //el.xs_gr->Draw("pcol");
  
  c->cd(3);
  gPad->SetLogy(0);
  gPad->SetLogz(1);

  el.xs_gr->GetXaxis()->SetRangeUser(700, 2400);
  el.xs_gr->Draw("colz");
  
  
  // =========== plot difference:
  c->cd(4);
  gPad->SetLogy(0);
  gPad->SetLogz(0);

  el.diff_gr->SetMinimum(-0.1) ; 
  el.diff_gr->SetMaximum(0.1) ; 
  el.diff_gr->Draw("colz");
  // gStyle->SetNumberContours(1) ; 
  // el.diff_gr->Draw("CONTsame");
  c->Print("lim2d.png");
  
  std::cout << "Checkpoint 3" << std::endl ; 
  
  TList* list = el.diff_gr->GetContourList(0);
  std::cout << list->GetSize() << std::endl ; 
  TGraph* cont0 = (TGraph*) list->At(0) ;
  cont0->SetNameTitle("limit-cont", "Exclusion mass region;M_{W}, GeV;M_{#nu}, GeV");
  // cont0->GetXaxis()->SetLimits(650, xmax);
  cont0->GetXaxis()->SetRangeUser(650, xmax);
  cont0->GetYaxis()->SetRangeUser(0, ymax); // 1000);

  TGraph* contourPlot = new TGraph() ; 
  std::vector< std::pair<double,double> > pts ; 
  for (int i=0; i<list->GetSize(); i++) { 
    // std::cout << "*** " << i << " ***" << std::endl ;
    TObjLink* lnk = list->FirstLink() ; 
    while (lnk) { 
      int nPts = ((TGraph*)lnk->GetObject())->GetN() ; 
      bool nonNan = false ; 
      for (int j=0; j<nPts; j++) { 
	Double_t xval, yval ; 
	((TGraph*)lnk->GetObject())->GetPoint(j,xval,yval) ; 
	if ( xval == xval ) { // nan fix
	  if ( xval < xmax ) { // kluge...very afraid of this
	    if ( xval - yval > 100 ) { 
	      if ( xval > 700 ) { // Tevatron exclusion region
		nonNan = true ; 
		std::pair<double,double> thePoint ; 
		thePoint.first = xval ; 
		thePoint.second = yval ; 
		bool duplicate = false ; 
		for (unsigned int pItr=0; pItr<pts.size(); pItr++) 
		  if (thePoint == pts.at(pItr)) duplicate = true ;
		if ( !duplicate ) pts.push_back( thePoint ) ; 
		
		// if ( xval < 800 ) std::cout << "(" << xval << "," << yval << ") " ; 
		// contourPlot->SetPoint(pointCounter++,xval,yval) ; 
	      }
	    }
	  }
	}
      }
      // if ( nonNan ) std::cout << std::endl ; 
      lnk = lnk->Next() ; 
    }
  }
  std::cout << "There are " << pts.size() << " points to consider" << std::endl ; 
  unsigned int iMax = 0 ; 
  double maxXval = 0 ; 
  for (unsigned int i=0; i<pts.size(); i++) {
    std::cout << pts.at(i).first << "\t" << pts.at(i).second << std::endl ; 
    if ( pts.at(i).first > maxXval ) { maxXval = pts.at(i).first ; iMax = i ; } 
  }
  std::cout << "max x-val: " << maxXval << std::endl ; 

  // Create two sorted vectors 
  std::vector< std::pair<double,double> > highPts ; 
  std::vector< std::pair<double,double> > lowPts ; 
  for (unsigned int i=0; i<pts.size(); i++) {
    double yMax = pts.at(iMax).second ; 
    if ( pts.at(i).second >= yMax ) highPts.push_back( pts.at(i) ) ; 
    else lowPts.push_back( pts.at(i) ) ; 
  }

  std::cout << "Sorting" << std::endl ; 
  std::sort( highPts.begin(),highPts.end(),sortByXval ) ; 
  std::sort( lowPts.begin(),lowPts.end(),sortByXval ) ; 
  std::cout << "Done sorting" << std::endl ; 

  // Make sure you start exactly at 740 (Tevatron limit)
  // Hide the line at 740, though
  double newx = 735 ; 
#if 0
  std::cout << highPts.at(0).first << " " << lowPts.at(0).first << " " << newx << std::endl ;
  if ( highPts.at(0).first > newx ) { 
    double slope = ( highPts.at(0).second - highPts.at(1).second ) / 
      ( highPts.at(0).first - highPts.at(1).first ) ; 
    double newy = highPts.at(1).second - slope * ( highPts.at(1).first - newx ) ;
    highPts.at(0).first = newx ; 
    highPts.at(0).second = newy ; 
    std::cout << "high: " << slope << " " << newy << std::endl ; 
  }
#endif
  std::cout << "Checking low: " << (lowPts.at(0).first < newx) << std::endl ; 
  if ( lowPts.at(0).first < newx ) { 
    std::cout << "In here" << std::endl ;
    std::cout << lowPts.size() << std::endl ; 
    std::cout << lowPts.at(0).first << " " << lowPts.at(0).second << std::endl ; 
    std::cout << lowPts.at(1).first << " " << lowPts.at(1).second << std::endl ; 
    std::cout << lowPts.at(0).second - lowPts.at(1).second << std::endl ; 
    std::cout << lowPts.at(0).first - lowPts.at(1).first << std::endl ; 
    double slope = ( lowPts.at(0).second - lowPts.at(1).second ) / 
      ( lowPts.at(0).first - lowPts.at(1).first ) ; 
    std::cout << "slope: " << slope << std::endl ;
    double newy = lowPts.at(1).second - slope * ( lowPts.at(1).first - newx ) ;
    std::cout << slope << " " << newy << std::endl ;
    lowPts.at(0).first = newx ; 
    lowPts.at(0).second = newy ; 
    std::cout << "low: " << slope << " " << newy << std::endl ; 
  }

  std::cout << "Still here" << std::endl ; 
  int pointCounter = 0 ; 
  std::cout << "Adding pts to graph" << std::endl ; 

  contourPlot->SetPoint(pointCounter++,740,640) ;
  for (unsigned int i=0; i<highPts.size(); i++) {
    double x=highPts.at(i).first;
    double y=highPts.at(i).second;
    if (y > x - 100) y = x - 100;
    contourPlot->SetPoint(pointCounter++,x,y) ;
    std::cout << x << "\t" << y << std::endl ; 
  }
  std::cout << "High pts added" << std::endl ; 

  for (unsigned int i=0 ; i<lowPts.size(); i++) {
    unsigned int j = lowPts.size() - i - 1 ; 
    contourPlot->SetPoint(pointCounter++,lowPts.at(j).first,lowPts.at(j).second) ;
    std::cout << lowPts.at(i).first << "\t" << lowPts.at(i).second << std::endl ; 
  }
  contourPlot->SetPoint(pointCounter++,740,100) ;
  std::cout << "Low pts added" << std::endl ; 

  // list->At(0)->Print() ; 

  // show exclusion zone
  // cont0->SetLineColor(kBlue - 3);

  return contourPlot;
}                                                          // genContour

//======================================================================

void lim2d(const char *filename = "../10kcycles.csv",
	   //const char *filename = "Data_CS_Up_Lim.txt",
	   const char *obsscanfmt  = "%lg %lg %lg", // to plot observed limit
	   const char *expscanfmt  = "%lg %lg %*lg %lg") // to plot expected limit
{
  //gROOT->SetStyle("Plain");
  setTDRStyle();
  //gStyle->SetTitleOffset(1.2,"Y");
  //gStyle->SetPadRightMargin(0.2);
  //gStyle->SetOptTitle(1);
  //gStyle->SetTitleX(.1);
  gStyle->SetPalette(1);
  // fixOverlay() ; 
  
  // ============ plot contour(s):
  TCanvas *c2 = new TCanvas("c2","c2",800,800);

  c2->SetLeftMargin(0.20);
  c2->SetRightMargin(0.05);
  c2->SetLogy(0);
  gStyle->SetTitleOffset(1.5, "Y");
  //gStyle->SetOptTitle(0);

  TH2D* grid = new TH2D("grid","grid",100,650,xmax,100,0,ymax); //1000) ; 
  grid->GetXaxis()->SetTitle("M_{W_{R}} [GeV]") ; 
  grid->GetYaxis()->SetTitle("M_{N_{#mu}} [GeV]") ; 
  grid->Draw() ; 

  TGraph *contourObs = genContour(filename,obsscanfmt,"obs");

  c2->cd();

  contourObs->SetLineColor(kBlue);
  contourObs->SetFillColor(kBlue-4);
  contourObs->SetFillStyle(3005);
  contourObs->SetLineWidth(5);
  contourObs->Draw("F SAME") ;
  contourObs->Draw("L SAME") ;

  if( expscanfmt && strlen(expscanfmt) ) {
    TGraph *contourExp = genContour(filename,expscanfmt,"exp");

    c2->cd();

    contourExp->SetLineColor(kBlue);
    contourExp->SetLineStyle(2);
    contourExp->SetLineWidth(5);
    contourExp->Draw("L SAME") ;

    TLegend *leg = new TLegend(.8,.8,.95,.93);
    leg->SetTextFont(42);
    leg->SetBorderSize(1);
    leg->SetFillColor(10);
    leg->SetLineWidth(1);
    leg->SetNColumns(1);
    leg->AddEntry(contourObs,"Observed","L");
    leg->AddEntry(contourExp,"Expected","L");

    leg->Draw();
  }  
  // std::cout << "Hello" << std::endl ; 

  // TBox* box2 = new TBox(760,750,900,850) ; 
  // box2->SetLineColor(kWhite) ; 
  // box2->SetFillColor(kWhite) ; 
  // box2->Draw() ; 

  // std::cout << "Hello" << std::endl ; 

  // TBox* box3 = new TBox(1200,100,1400,200) ; 
  // box3->SetLineColor(kWhite) ; 
  // box3->SetFillColor(kWhite) ; 
  // box3->Draw() ; 

  // std::cout << "Hello" << std::endl ; 

  // TBox* box4 = new TBox(1320,100,1400,320) ; 
  // box4->SetLineColor(kWhite) ; 
  // box4->SetFillColor(kWhite) ; 
  // box4->Draw() ; 

  // std::cout << "Hello" << std::endl ; 

  // TBox* box5 = new TBox(1350,700,1400,900) ; 
  // box5->SetLineColor(kWhite) ; 
  // box5->SetFillColor(kWhite) ; 
  // box5->Draw() ; 

  // std::cout << "Hello" << std::endl ; 

  // contourPlot->SetLineColor(kBlue);
  // contourPlot->SetFillColor(kBlue);
  // contourPlot->SetLineWidth(3);
  // contourPlot->DrawCopy("L SAME") ;


  // grid->Draw("A SAME") ; 

  // cont0->Draw("AL");
  // contourPlot->SetLineWidth(4) ; 
  // cont0->Draw("A");
  // el.diff_gr->GetContourList(0)->Draw("same") ; 

  std::cout << "Checkpoint 4" << std::endl ; 

  std::cout << "Checkpoint 5" << std::endl ; 
  // ============== plot Tevatron limit region:

  float x[] = {650, 740, 740, 650};
  float y[] = {0, 0, 740, 650};
  TGraph* Tevatron = new TGraph(4, x, y);
  Tevatron->SetFillColor(kGray+1) ; 
  Tevatron->SetLineColor(kGray) ; 
  Tevatron->Draw("F SAME") ; 
  
  // ============== plot region M_nuR > M_WR:
  float x2[] = {650, ymax, 650} ;
  float y2[] = {ymax, ymax, 650} ;
  TGraph* wrnu = new TGraph(3, x2, y2);
  wrnu->SetLineWidth(3);
  wrnu->SetLineColor(kYellow-4);
  wrnu->SetFillColor(kYellow);
  wrnu->Draw("F SAME") ; 

  // txt_gr->SetLineWidth(2 + 15*100);
  // txt_gr->SetFillStyle(3004);
  // txt_gr->Draw("L");
  
  std::cout << "Checkpoint 6" << std::endl ; 
  // ============= plot labels:

  TLatex text;
  text.SetTextSize(0.04);
  text.DrawLatex(680, .94*ymax, "M_{N_{#mu}} > M_{W_{R}}");
  
  text.SetTextSize(0.035);
  text.SetTextAngle(90);
  text.DrawLatex(730, 30, "Excluded by Tevatron");

  // watermark
  TLatex* mark  = new TLatex( 0.67,0.245, "CMS Preliminary" ) ;
  TLatex* mark2 = new TLatex( 0.625,0.21, Form("#int Ldt = %.0f pb^{-1} at 7 TeV",luminvpb) ); 
  mark->SetNDC(kTRUE) ;
  mark->SetTextSize(0.03) ;
  mark->Draw() ;
  mark2->SetNDC(kTRUE) ;
  mark2->SetTextSize(0.03) ;
  mark2->Draw() ;

  //fixOverlay() ; 
  
  std::cout << "Checkpoint 7" << std::endl ; 
  c2->Print("lim2d-cont0.pdf");
  c2->Print("lim2d-cont0.png");
  c2->Print("lim2d-cont0.eps");
  c2->Print("lim2d-cont0.C");
}
