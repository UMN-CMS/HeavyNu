#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"

#include <iostream>

double accRatio(int mwr) {
  return 1.0;
  // (from the acceptance db) for generated versus MN=MW/2
  int roundmwr=(mwr/100)*100;
  switch (mwr) {
  case (1000) : return 0.605/0.630;
  case (1100) : return 0.637/0.657;
  case (1200) : return 0.678/0.691;
  case (1300) : return 0.696/0.708;
  case (1400) : return 0.713/0.729;
  case (1500) : return 0.729/0.742;
  case (1600) : return 0.740/0.759;
  case (1700) : return 0.749/0.767;
  case (1800) : return 0.766/0.774;
  case (1900) : return 0.772/0.774;
  case (2000) : return 0.778/0.789;
  case (2100) : return 0.791/0.794;
  case (2200) : return 0.794/0.793;
  case (2300) : return 0.787/0.800;
  case (2400) : return 0.799/0.808;
  case (2500) : return 0.803/0.805;
  default: 
    return accRatio(roundmwr)+(mwr-(roundmwr+100))*(accRatio(roundmwr+100)-accRatio(roundmwr))/100.0;
    //return 1.0;
  }
}

void smoothSG(float pts[], int n, int mode=4) {
  // − 3y1 + 12y2 + 17y3 + 12y4 − 3y5) / 35

  float ptsOut[1000];
 
  //  const int mode=19;
  mode+=10;

  if (mode==1) {
    for (int i=0; i<n; i++) {
      double sum=0;
      
      if (i>1) sum+=-3*pts[i-2];
      else if (i>0) sum+=-3*pts[i-1];
      else sum+=-3*pts[i];
      
      if (i>0) sum+=12*pts[i-1];
      else sum+=12*pts[i];
      
      sum+=17*pts[i];
      
      if (i<n-1) sum+=12*pts[i+1];
      else sum+=12*pts[i];
      
      if (i<n-2) sum+=-3*pts[i+2];
      else if (i<n-1) sum+=-3*pts[i+1];
      else sum+=-3*pts[i];
      
      sum/=35;
      ptsOut[i]=sum;
    }
  }
  if (mode==2) {
    for (int i=0; i<n; i++) {
      double sum=0;
      
      if (i>1) sum+=pts[i-2];
      else if (i>0) sum+=pts[i-1];
      else sum+=pts[i];
      
      if (i>0) sum+=pts[i-1];
      else sum+=pts[i];
      
      sum+=pts[i];
      
      if (i<n-1) sum+=pts[i+1];
      else sum+=pts[i];

      
      if (i<n-2) sum+=pts[i+2];
      else if (i<n-1) sum+=pts[i+1];
      else sum+=pts[i];
      
      sum/=5;
      ptsOut[i]=sum;  
    }
  }
  if ((mode/10)==1) {
    int window=mode%10;


    float x[100],y[100];
    for (int i=0; i<=2*window; i++) x[i]=i-window;

    TF1 f1("f1","pol1",-window+0.5,window-0.5);

    for (int i=0; i<n; i++) {

      for (int j=-window; j<=window; j++) {
	if ((i+j)<0) y[j+window]=pts[0];
	else if ((i+j)>=n) y[j+window]=pts[n-1];
	else y[j+window]=pts[i+j];
      }

      TGraph tg(window*2+1,x,y);

      tg.Fit(&f1,"QR","",-window+0.5,+window-0.5);
      
      ptsOut[i]=f1.Eval(0);
    }
  }


  for (int i=0; i<n; i++) pts[i]=ptsOut[i];
}

#include "tdrstyle.C"

void limWR(const char* fname,int which, int smooth=0,const char* asUsed=0,double lumi=19.7) {
  float mw[1000], mn[1000], mwt[1000];
  float obs[1000], exp[1000], expm2s[1000],expm1s[1000],expp1s[1000],expp2s[1000],xsec[1000];
  float obs_ns[1000],exp_ns[1000];
  int n=0,nx=0;
  char buffer[1024];
  //  double pbr=1e3; // pb ratio
  double pbr=1e3; // pb ratio
  double corr=1.0;

  gStyle->SetHatchesLineWidth(3);
  setTDRStyle();
  gStyle->SetHatchesLineWidth(3);

  int ecm;
  if (which==0) { // 2011 Muons
    ecm=7;
  } else if (which==1) { // 2012 Muons
    ecm=8;
  } else if (which==2) { // 2012 Electrons
    ecm=8;
  } else if (which==3) { // 2012 Muons+Electrons
    ecm=8;
  } else if (which==4) { // USe the 8 TeV xsec as master
    ecm=8;
  }


  bool correctToMWMN2=true;

  FILE* f1=fopen(fname,"rt");
  int dummyi;

    while (!feof(f1))
    {
        buffer[0] = 0;
        fgets(buffer, 1000, f1);
        if (sscanf(buffer, "%f %f %f %f %f %f %f %f ", mw + n, mn + n,
                   obs + n, exp + n, expm2s + n, expm1s + n, expp1s + n, expp2s + n) == 8)
        {
            mw[n] /= 1000.0;
            if (correctToMWMN2)
            {
                corr = 1.0 / accRatio(mw[n]);
                mn[n] = mw[n] / 2.0;
            }

            if (which == 3)
            {
                corr = 1.0;// we now apply this to the xsec 3.0 / 2.0; // correct to all leptons
            }


            obs[n] *= pbr*corr;
            exp[n] *= pbr*corr;
            expm2s[n] *= pbr*corr;
            expm1s[n] *= pbr*corr;
            expp1s[n] *= pbr*corr;
            expp2s[n] *= pbr*corr;
            n++;
        }
    }
    fclose(f1);

  FILE* fx=fopen("cs.txt","rt");
  float rawx, kf,xmw,xmn;
  int ecmx;
  while (!feof(fx)) {
    buffer[0]=0;
    fgets(buffer,1000,fx);
    if (sscanf(buffer,"%d %f %f %f %f",&ecmx, &xmw,&xmn,&rawx,&kf)==5) {
        xmw /= 1000.0;
        xmn /= 1000.0;
      if (ecmx!=ecm) continue;
      for (int search=0; search<n; search++) {
	if (fabs(xmw-mw[search])<1.0/1000.0 && fabs(xmn-mn[search])<1.0/1000.0) {
	  xsec[nx]=rawx*kf*pbr;
	  mwt[nx]=xmw;
      std::cout << mwt[nx] << std::endl;

	  if (which==3) {
	    xsec[nx]*=8.0/3.0 * 2.0/3; // from Pythia (2/3 is to remove tau contribution)
	  }

	  //      xsec[nx]/=0.75;
	  nx++;
	  break;
	}
      }
    }
  }
  fclose(fx);

  if (which==4) {

    for (int i=0; i<n; i++) {

      double ax=1;
      for (int search=0; search<nx; search++) {
	if (fabs(mw[i]-mwt[search])<1.0) {
	  ax=xsec[search];
	  xsec[search]=1.0;
	  break;
	}
      }
  
      obs[i]/=ax;
      exp[i]/=ax;
      expm2s[i]/=ax;
      expm1s[i]/=ax;
      expp1s[i]/=ax;
      expp2s[i]/=ax;
      
    }

  }
  
  //  gROOT->SetStyle("Plain");
  TCanvas* c1=new TCanvas("c1","c1",800,800);
  c1->SetTopMargin(0.05);
  c1->SetLeftMargin(0.15);
  c1->SetRightMargin(0.05);
  c1->SetBottomMargin(0.13);
  c1->SetTicks(1,1);
  c1->SetLogy();

  //  TH1* dummy=new TH1F("dummy","",30,1000,3000);
  TH1* dummy;
  if (which!=0) dummy=new TH1F("dummy","",30,1.0,3.5);
  else dummy=new TH1F("dummy","",30,1000,2600);

  //  dummy->SetMinimum(0.5);
  dummy->SetMinimum(0.2);
  if (which==4) {
   dummy->SetMaximum(10);
   dummy->SetMinimum(0.005);
  } else if (which==3) {
    dummy->SetMaximum(1000);
  } else {
    dummy->SetMaximum(50);
  }
  dummy->SetStats(0);
  dummy->GetXaxis()->SetNdivisions(507);
  dummy->GetXaxis()->SetTitle("M_{W_{R}} [TeV]");
  
  if (which==0 || which==1) {
    dummy->GetYaxis()->SetTitle("#sigma(pp#rightarrow W_{R}) #times BR(W_{R}#rightarrow #mu#mujj) [fb]");
    dummy->GetYaxis()->SetRangeUser(0.2, 100);
  } else if (which==2) {
    dummy->GetYaxis()->SetTitle("#sigma(pp#rightarrow W_{R}) #times BR(W_{R}#rightarrow eejj) [fb]");
    dummy->GetYaxis()->SetRangeUser(0.2, 100);
  } else if (which==3) {
    dummy->GetYaxis()->SetTitle("#sigma(pp#rightarrow W_{R}) #times BR(W_{R}#rightarrow (ee+#mu#mu)jj) [fb]");
    dummy->GetYaxis()->SetTitleSize(0.055);
    dummy->GetYaxis()->SetRangeUser(0.2, 100);
  } else if (which==4) {
    dummy->GetYaxis()->SetTitle("#sigma(pp#rightarrow W_{R}) #times BR(W_{R}#rightarrow #mu#mujj)/#sigma_{g_{R}=g_{L}}");
  }
  dummy->GetYaxis()->SetTitleOffset(1.1);
  dummy->Draw("HIST");

  TLegend* tl;
  if (which==0 || which==1) {
    tl=new TLegend(0.55,0.66,0.99,0.87);//,"CL_{S} Method          M_{N_{#mu}}=M_{W_{R}}/2");
    //tl->SetHeader("CL_{S} Method    95% CL");
    // tl->SetHeader("CL_{S} Method    M_{#scale[1.25]{N_{#scale[1.5]{#mu}}}}= M_{#scale[1.25]{W_{R}}}/2");
    //tl->SetHeader("M_{#scale[1.25]{N_{#scale[1.5]{#mu}}}}= M_{#scale[1.25]{W_{R}}}/2");
    TLatex *latex = new TLatex(0.57, 0.89, "M_{#scale[1.25]{N_{#scale[1.5]{#mu}}}}= M_{#scale[1.25]{W_{R}}}/2");
    latex->SetNDC();
    latex->SetTextSize(0.032);
    latex->SetTextFont(42);
    latex->Draw();
  } else if (which==3) {
    //tl=new TLegend(0.60,0.55,0.99,0.90,"CL_{S} Method          M_{#scale[1.25]{N}}=M_{#scale[1.25]{W_{R}}}/2");
      tl=new TLegend(0.535,0.66,0.975,0.88);//,"M_{#scale[1.25]{N}}=M_{#scale[1.25]{W_{R}}}/2");
      TLatex *latex = new TLatex(0.555, 0.895, "M_{#scale[1.25]{N}} = M_{#scale[1.25]{W_{R}}}/2");
      latex->SetNDC();
      latex->SetTextSize(0.032);
      latex->SetTextFont(42);
      latex->Draw();
  } else if (which==2) {
    //tl=new TLegend(0.55,0.55,0.99,0.90,"CL_{S} Method          M_{#scale[1.25]{N_{#scale[1.5]{e}}}}=M_{#scale[1.25]{W_{R}}}/2");
      tl = new TLegend(0.55,0.66,0.99,0.87);//,"M_{#scale[1.25]{N_{#scale[1.5]{e}}}}=M_{#scale[1.25]{W_{R}}}/2");
      TLatex *latex = new TLatex(0.57, 0.89, "M_{#scale[1.25]{N_{#scale[1.5]{e}}}}= M_{#scale[1.25]{W_{R}}}/2");
      latex->SetNDC();
      latex->SetTextSize(0.032);
      latex->SetTextFont(42);
      latex->Draw();
  } else if (which==4) {
    tl=new TLegend(0.16,0.90,0.90,0.72,"CL_{S} Method          M_{N_{#mu}}= M_{#scale[1.25]{W_{R}}}/2");
    tl=new TLegend(0.16,0.90,0.90,0.72,"M_{N_{#mu}}= M_{#scale[1.25]{W_{R}}}/2");
    tl->SetNColumns(2);
  }
  tl->SetTextFont(42);
  tl->SetTextSize(0.032);
  tl->SetFillStyle(0);
  tl->SetBorderSize(0);


  if ((smooth/100)%10) {
    smoothSG(expm2s,n);
    smoothSG(expp2s,n);
  }

  if ((smooth/10)%10) {
    smoothSG(expm1s,n);
    smoothSG(expp1s,n);
  }

  if ((smooth/1)%10) {
    smoothSG(exp,n);
  }

  if ((smooth/1000)%10) {
    smoothSG(obs,n,(smooth/1000)%10);
  }

  float mw2[200],e2s[200],e1s[200];
  for (int i=0; i<n; i++) {
    mw2[i]=mw[i];
    mw2[n+i]=mw[n-i-1];
    e2s[i]=expm2s[i];
    e2s[n+i]=expp2s[n-i-1];
    e1s[i]=expm1s[i];
    e1s[n+i]=expp1s[n-i-1];
    //    xsec[i]*=0.75; // 0.75 is average reduced xsection for a single channel * k factor (determined at MW=1.3 TeV)
  }

  TGraph* tg_e2s=new TGraph(n*2,mw2,e2s);
  tg_e2s->SetLineWidth(0);
  tg_e2s->SetFillColor(kYellow);
  tg_e2s->Draw("F SAME");

  TGraph* tg_e1s=new TGraph(n*2,mw2,e1s);
  tg_e1s->SetLineWidth(0);
  tg_e1s->SetFillColor(kGreen);
  tg_e1s->Draw("F SAME");

  TGraph* tg_theory=new TGraph(nx,mwt,xsec);
  tg_theory->SetLineWidth(3);
  tg_theory->SetLineColor(kRed+2);
  tg_theory->SetLineStyle(0);
  tg_theory->SetFillStyle(3002);
  tg_theory->SetFillColor(kRed);
  
  
  double tesx[] = {1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9};
  //double tesy[] = {0.072374,0.069254,0.067429,0.066900,0.067665,0.069726,0.073082,0.077733,0.083679,0.090920,0.099457,0.109289,0.120415,0.132838,0.146555,0.161567,0.177875,0.195478,0.214375,0.234569,0.256057,0.278840,0.302919,0.328293,0.354962,0.382926,0.412185,0.442740,0.474589};
  double tesy[] = {0.051550, 0.056700, 0.061725, 0.066400, 0.071525, 0.077125, 0.082925, 0.089725, 0.097000, 0.103725, 0.112475, 0.121300, 0.131475, 0.142200, 0.151650, 0.163925, 0.178325, 0.192050, 0.211700, 0.235250, 0.267625, 0.292425, 0.328975, 0.34375, 0.3775, 0.4152, 0.4496, 0.0, 0.0, 0.0};
  TGraph* tg_theoryerror=new TGraph();
  int iValidPoint = 0, ivps1m;
  for(int i = 0; i < 29; i++)
  {
      //if(fabs(tesx[i] - mwt[iValidPoint]) > 1.0/1000.0) continue;
      tg_theoryerror->SetPoint(iValidPoint, tesx[i], xsec[iValidPoint]*(1-tesy[i]));
      std::cout << tesx[i] << "\t" << mwt[iValidPoint] << std::endl;
      iValidPoint++;
  }
  ivps1m = iValidPoint;
  for(int i = 28; i >= 0; i--)
  {
      if(fabs(tesx[i] - mwt[2*ivps1m-iValidPoint - 1]) > 1.0/1000.0) continue;
      tg_theoryerror->SetPoint(iValidPoint, tesx[i], xsec[2*ivps1m-iValidPoint - 1]*(1+tesy[i]));
      iValidPoint++;
  }
  tg_theoryerror->SetLineWidth(3);
  tg_theoryerror->SetLineColor(kRed+2);
  tg_theoryerror->SetFillColor(kRed);
  tg_theoryerror->SetFillStyle(3002);
  tg_theoryerror->Draw("F SAME");
  
  tg_theory->Draw("L SAME");
  
  TGraph* tg_exp=new TGraph(n,mw,exp);
  tg_exp->SetLineWidth(2);
  tg_exp->SetLineColor(kBlue);
  tg_exp->SetLineStyle(2);
  tg_exp->Draw("L SAME");

  TGraph* tg_obs=new TGraph(n,mw,obs);
  tg_obs->SetLineWidth(2);
  tg_obs->Draw("L SAME");

  tl->AddEntry(tg_obs,"Observed limit","L");
  tl->AddEntry(tg_exp,"Expected limit","L");
  tl->AddEntry(tg_e1s,"Expected #pm 1 #sigma","F");
  tl->AddEntry(tg_e2s,"Expected #pm 2 #sigma","F");
  
  if (which==3) {
    tl->AddEntry(tg_theory,"Theory","LF");
    tl->AddEntry((TObject*)0,"g_{R}= g_{L}, M_{#scale[1.25]{N}_{#scale[1.5]{e}}}= M_{#scale[1.25]{N}_{#scale[1.5]{#mu}}}= M_{#scale[1.25]{N}_{#scale[1.5]{#tau}}}","");
  } else {
    tl->AddEntry(tg_theoryerror,"Theory (g_{R}= g_{L})","LF");
  }

  TText *text;
  if (which==0) {
    
    text = new TText(0.89,0.97,"CMS");
  } else {
    text = new TText(0.72,0.97,"CMS Preliminary");
  }
  //TText *text = new TText(0.89,0.94,"CMS");
   text->SetNDC();
   text->SetTextFont(42);
   //   text->SetTextSize(0.03);
   text->SetTextSize(0.05);
   //   text->Draw();

   //text = new TLatex(0.2,0.2,"M_{N_{#mu}}= M_{W_{R}}/2");
   //text->SetNDC();
   //text->SetTextFont(42);
   //text->SetTextSize(0.05);
   //text->Draw();

   /*
   TDatime now;
   sprintf(buffer,"%d/%d/%d",now.GetDay(),now.GetMonth(),now.GetYear());

   text = new TText(0.80,0.88,buffer);
   text->SetNDC();
   text->SetTextSize(0.03);
   text->Draw();
   */
   /*
   if (which==3) {
     text = new TText(0.16,0.16,"No systematics");
     text->SetNDC();
     text->SetTextColor(kBlue);
     text->SetTextSize(0.04);
     text->Draw();
   }
   */

   

   
   if (which==4) {
     //TLatex *   texa = new TLatex(0.152,0.972,"#int");
     //texa->SetTextFont(42);
     //texa->SetNDC();
     //texa->SetTextSize(0.017);
     //texa->SetLineWidth(2);
     //texa->Draw();
     sprintf(buffer,"CMS Preliminary %.1f fb^{-1} at #sqrt{s}=%d TeV + %.1f fb^{-1} at #sqrt{s}=%d TeV",lumi,ecm,5.0,7);
     TLatex *   tex = new TLatex(0.150,0.950,buffer);
     tex->SetNDC();
     tex->SetTextSize(0.032);
     tex->SetLineWidth(2);
     tex->Draw();

     //sprintf(buffer,"%.1f fb^{-1} at #sqrt{s}=%d TeV",5.0,7);
     //tex = new TLatex(0.47,0.970,buffer);
     //tex->SetNDC();
     //tex->SetTextSize(0.026);
     //tex->SetLineWidth(2);
     //tex->Draw();
   } else {
     /*
     TLatex *   texa = new TLatex(0.152,0.972,"#int");
     texa->SetNDC();
     texa->SetTextFont(42);
     texa->SetTextSize(0.017);
     texa->SetLineWidth(2);
     texa->Draw();
     */
     //sprintf(buffer,"CMS    #sqrt{s} = %d TeV    %.1f fb^{-1}",ecm,lumi);
     //sprintf(buffer,"        CMS    #sqrt{s} = %d TeV    %.1f fb^{-1}",ecm,lumi);
     //sprintf(buffer, "CMS    #sqrt{s} = %d TeV    L = %0.1f fb^{-1}",ecm,lumi);
     sprintf(buffer, "%0.1f fb^{-1} (%d TeV)", lumi, ecm);
     //mark->DrawLatex(0.220, 0.9575, "        CMS    #sqrt{s} = 8 TeV    19.7 fb^{-1}");
     TLatex *   tex = new TLatex(c1->GetLeftMargin()+(1-c1->GetLeftMargin()-c1->GetRightMargin())/2, 0.9575,buffer);
     tex->SetNDC();
     tex->SetTextFont(42);
     tex->SetTextAlign(31);
     tex->SetTextSize(0.04 * 1.1);
     tex->DrawLatex(1 - c1->GetRightMargin(), 0.9575, buffer);
     tex->SetTextAlign(13);
     tex->SetTextSize(0.04 * 1.1 * 1.25);
     tex->SetTextFont(61);
     //tex->DrawLatex(c1->GetLeftMargin() + 0.025, 1 - (c1->GetTopMargin() + 0.025), "#splitline{CMS}{#scale[0.7]{#it{Preliminary}}}");
     tex->DrawLatex(c1->GetLeftMargin() + 0.025, 1 - (c1->GetTopMargin() + 0.025), "CMS");
     //tex->SetNDC();
     //tex->SetTextAlign(21);
     //tex->SetTextFont(42);
     //tex->SetTextSize(0.04 * 1.1);
     //tex->SetLineWidth(2);
     //tex->Draw();
   }

   /*

   tex = new TLatex(0.75,0.5,"M(N_{R})=#frac{1}{2} M(W_{R})");
   tex->SetNDC();
   tex->SetTextSize(0.03);
   tex->Draw();
   */

     
   tl->Draw("SAME");
  c1->RedrawAxis();
  char feps[100];
  char fpng[100];
  char* fbase="limWR";

  if (which==4) fbase="limWR_2y";
  if (which==3) fbase="limWR_emu";
  if (which==2) fbase="limWR_e12";
  if (which==1) fbase="limWR_mu12";
  if (which==0) fbase="limWR_mu11";

  sprintf(feps,"%s.eps",fbase);
  sprintf(fpng,"%s.png",fbase);
  c1->Print(feps);
  c1->Print(fpng);
  char cmd[100];
  sprintf(cmd,"epstopdf %s",feps);
  system(cmd);

  if (asUsed!=0) {
    FILE* f=fopen(asUsed,"w");

    for (int i=0; i<n; i++) {
      fprintf(f,"%4d %4d %7f %7f %7f %7f %7f %7f 0 \n",
	      int(tg_obs->GetX()[i]),int(tg_obs->GetX()[i])/2,
	      tg_obs->GetY()[i],tg_exp->GetY()[i],
	      expm2s[i],expm1s[i],expp1s[i],expp2s[i]);

    }
    fclose(f);
  }

}
