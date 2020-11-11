
TString graphOpt = "";
TString graphTitle = "";
bool drawText = true;
int lineColor = kBlue;
int markerType = 20; 
double markerSize = 1.2; 

TGraph *  PlotDeviations(int N=256, int type = 0, int nskip = 0) {
   // plot the deviations for the different engine

   //dirName = "Deviations_Random1";
   double xfmin = 0;
   double xfmax = 10;

   double xgmax = 10.5;  // max of graphs
   
   double lambda_max = 1; 
   int ntrials = 10000;

   //int type = 0; // for different types giving N
   
   //int N = 10;
   if (N==10) { 
      lambda_max = TMath::Log2(21.3497);
      ntrials = 8000;
      nskip = 0;
      xfmin=0;
      xfmax=14.5;
      xgmax = 20.5;
   }
   if (N==16) { 
      lambda_max = TMath::Log2(41.697);
      ntrials = 20000;
      nskip = 0;
      xfmin=0;
      xfmax=11.5;
      xgmax = 20.5;
   }
   if (N== 8) {
      //lambda_max = TMath::Log2( TMath::Power(2.,36)+1);
      // obtained from Mathematica
      lambda_max = TMath::Log2( 2.042572776E9);
      ntrials = 3281;
      //xfmin=0.5;
      xfmax=4.5;
   }
   if (N== 8 && type == 1) {
      // with m = 2^30+1
      //lambda_max = TMath::Log2( TMath::Power(2.,36)+1);
      // obtained from ROOT
      // entropy is 125
      lambda_max = TMath::Log2( 6.01919e+07);
      ntrials = 3281;
      //xfmin=0.5;
      xfmax=2.5;
   }

   if (N== 8 && type == 2) {
      // m = 2^24+1
      // entropy = 99.8
      lambda_max = TMath::Log2( 1.83548e+06);
      ntrials = 3281;
      //xfmin=0.5;
      xfmax=2.5;
   }
   if (N== 8 && type == 3) {
      // m = 2^18+1
      // entropy = 74.9
      lambda_max = TMath::Log2( 59766.6);
      ntrials = 3280;
      //xfmin=0.5;
      xfmax=3.5;
   }
   if (N== 8 && type == 4) {
      // m = 2^12+1
      // entropy = 49.9
      lambda_max = TMath::Log2( 2223.45);
      ntrials = 3280;
      //xfmin=0.5;
      xfmax=4.5;
   }
      
   if (N== 17) {
      //lambda_max = TMath::Log2( TMath::Power(2.,36)+1);
      lambda_max = TMath::Log2(2.226219317E10 );
      ntrials = 20000;
      //xfmin=0.5;
      xfmax=3.5;
   }
  if (N==240) {
     //lambda_max = TMath::Log2( TMath::Power(2.,32)+1);
      lambda_max = TMath::Log2( 3.882098668E11 );
      ntrials = 20000;
      //xfmin=0.5;
      xfmax=2.5;
   }
  if (N== 256) {
      lambda_max = TMath::Log2( 3002);
      ntrials = 2000;
      xfmin=0.;
      xfmax=5.5;
   }
   if (N== 44) {
      lambda_max = TMath::Log2( 183.671);
      ntrials = 20000;
      xfmin=0.;
      xfmax=8.;
   }
   if (N== 88) {
      lambda_max = TMath::Log2( 537.652);
      ntrials = 20000;
      xfmin=0.;
      xfmax=6.5;
   }
   if (N== 1000) {
      lambda_max = TMath::Log2( 29450.2);
      ntrials = 20000;
      xfmin=0.;
      xfmax=4.5;
   }
   if (N== 44851) {
      lambda_max = TMath::Log2( 29450.2);
      ntrials = 20000;
      xfmin=0.;
      xfmax=2.5;
   }



  
   
   TString fname;
   fname = TString::Format("Deviations_trials%d_Test::TestEngine<ROOT::Math::MixMaxEngine<%d,%d> >_N%d_test.root",ntrials,N,nskip,N);
   if (type != 0 && N == 8)  fname = TString::Format("Deviations_trials%d_Test::TestEngine<ROOT::Math::MixMaxEngine<%d,%d> >_N%d-%d_test.root",ntrials,N,nskip,N,type);

   std::cout << fname << std::endl;
      
   auto file = TFile::Open(fname);
   if (!file) return 0; 

   TString gname = TString::Format("deviations_N%d",N);

   TString title = TString::Format("Deviations for MIXMAX %d;Iteration Number;Distance in log2",N);

   auto g = (TGraphErrors*) file->Get(gname); 

   if (!g) return 0;

   std::cout << "Study " << g->GetTitle() << std::endl;

   if (graphTitle != "") title = graphTitle + ";Iteration Number;Distance in log2";
                                             
   g->SetTitle(title);
   g->SetMarkerStyle(markerType);                                               
   g->SetMarkerSize(markerSize);
                                             

   auto f1 = new TF1("f1","[Slope]*x + [Constant]",xfmin, xfmax);

   f1->SetParameters(1,1);
   f1->SetLineWidth(3);
   f1->SetLineColor(lineColor);
   f1->SetLineStyle(2);
   g->SetLineColor(lineColor+1);

                
   g->Fit(f1,"R");
   g->Draw(graphOpt);

   g->GetXaxis()->SetRangeUser(0,xgmax);
   g->GetXaxis()->CenterTitle(true);
   g->GetXaxis()->SetTitleSize(0.045);
   g->GetXaxis()->SetTitleOffset(0.95);
   g->GetXaxis()->SetNdivisions(10);

   g->GetYaxis()->CenterTitle(true);
   g->GetYaxis()->SetTitleSize(0.045);
   g->GetXaxis()->SetTitleOffset(0.8);
  
   
   double x1,y1;
   g->GetPoint(1,x1,y1);
   double e1 = g->GetErrorY(1);

   if (drawText) { 
      TObject * o = nullptr; 
      // auto legend  = new TLegend(0.5,0.3,0.87,0.5);
      //auto legend = new TPaveText(0.5,0.3,0.87,0.5,"NDC");
      auto legend = new TPaveText(0.4,0.3,0.87,0.5,"NDC");  // for 3 plots
      legend->AddText(TString::Format("log2(#lambda_{max} )= %8.5f",lambda_max));
      legend->AddText(TString::Format("Slope = %5.2f #pm %4.2f",f1->GetParameter(1),f1->GetParError(1)));
      legend->AddText(TString::Format("point1 = %5.2f #pm %4.2f",y1,e1));
      legend->Draw();
   }

   return g; 
}


void PlotHor(const std::vector<int>  vlist) {
   auto c1 = new TCanvas("mixmax","mixmax",1400,800);
   c1->Divide(vlist.size()/2,2);
   for (size_t i = 0; i < vlist.size(); ++i) { 
      c1->cd(i+1);
      PlotDeviations(vlist[i] );
   }
}


void PlotMany(const std::vector<int>  vlist) {

   drawText = false;
   graphTitle = "Deviations for Original MIXMAX";
   std::vector<int> lineColors = {kBlue, kBlack, kGreen, kRed, kCyan, kMagenta};
   std::vector<int> mTypes = { 20, 21, 22, 23, 29, 47};
   std::vector<double> mSizes = { 1.3, 1.3, 1.35, 1.45, 1.5, 1.5};
   std::vector<TGraph *> glist(vlist.size() ); 
   auto legend  = new TLegend(0.6,0.2,0.87,0.6);
   
   for (size_t i = 0; i < vlist.size(); ++i) {
      if (i > 0) graphOpt = "PL same";
      markerType = mTypes[i];
      markerSize = mSizes[i];
      lineColor = lineColors[i];
      auto g = PlotDeviations(vlist[i] );
      legend->AddEntry(g,TString::Format("MIXMAX %d",vlist[i]),"lp");
   }
   legend->Draw();
}

