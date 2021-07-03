/*--------------------------------------------------------------------------------------------------*\
|                                                                                                    |
|    Macro that plotes an arbitrary look-up table                                                    |
|                                                                                                    |
| Author: 10/2016, Benjamin Brudnyj                                                                  | 
\*--------------------------------------------------------------------------------------------------*/

void LUTToFunction() {

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Create variables
  ///________________________________
  static const Int_t Bins     = 2008;
  static const Int_t Register = 502;   // 2008/4 = 502
  static const Int_t n        = 4;     // 4 entries (PID-values) in each register
  ULong64_t ScaleRange        = 0;    
  
  char *LUT_file_name         = "lhc11dv3en";//"lhc11dv2en";
  char *file_type             = "tcs"; // the dot is set below
  
  ULong64_t NoBins            = Bins;
  Int_t scaleq0               = 0;    
  Int_t scaleq1               = 0;    

  Int_t addrPID[Bins]         = {0};  
  Int_t addr[Register]        = {0};   // 4 PID values in each register

  Double_t ChargeDep[Bins]    = {0};
  Int_t gChargeDep[Bins]      = {0};
  Int_t Iter                  = 0;  
  Double_t IterCD             = 0.;    

  Int_t Arbitrary[2]          = {0};
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Fill variables
  ///_______________________________________
  ifstream LUTfile;  
  LUTfile.open(Form("%sDraw.%s",LUT_file_name,file_type));
  if (LUTfile.is_open()) {

    LUTfile >> scaleq0;
    LUTfile >> scaleq1;
    LUTfile >> Arbitrary[0];
    LUTfile >> Arbitrary[1];

    ScaleRange = 20000;//(NoBins<<32) / scaleq0;
    //if((scaleq0 * ScaleRange)>>32 != NoBins) ScaleRange = ScaleRange + 1;

    IterCD = (Double_t)ScaleRange / (Double_t)NoBins;

    for(Int_t i = 0; i < Register; i++) {

      LUTfile >> addr[i];
      
      for(Int_t j = 0; j < n; j++) { 

        addrPID[Iter] = (addr[i] >> 8*j) & 0xFF;      

        if(i == (Register - 1) && j == 3) continue;
	ChargeDep[Iter+1] = ChargeDep[Iter] + IterCD; 

	Iter++;	
      }
    }
      
    LUTfile.close();
  }
  else cout << "Unable to open file";
  
  for(Int_t i = 0; i < Bins; i++) {
    gChargeDep[i] = (Int_t)ChargeDep[i];
  }
  
  TGraph *gLUT = new TGraph(Bins,gChargeDep,addrPID);
  gLUT->SetTitle(Form("LUT %s",LUT_file_name));
  gLUT->GetXaxis()->SetTitle("Charge deposition (a.u.)");
  gLUT->GetYaxis()->SetTitle("PID (a.u.)");
  gLUT->SetLineColor(4);
  gLUT->SetMaximum(280.);
  gLUT->SetMinimum(0.);

  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Draw and save graph
  ///_______________________________________
  TCanvas *c = new TCanvas("c",Form("LUT from file %s.%s",LUT_file_name,file_type));
  c->Divide(1,1);
  c->cd(1);
  gLUT->Draw("AL");

  
}
