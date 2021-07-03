/*--------------------------------------------------------------------------------------------------*\
|                                                                                                    |
|    Macro that translates an arbitrary function and sacle factor into a look-up table (LUT)         |
|                                                                                                    |
| Author: 10/2016, Benjamin Brudnyj                                                                  |
\*--------------------------------------------------------------------------------------------------*/


void FunctionToLUT() {

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Create variables
  ///________________________________
  static const Int_t Bins     = 2008;
  static const Int_t Register = 502; // 2008/4 = 502
  static const Int_t n        = 4;   // 4 entries (PID-values) in each register
  Int_t ScaleRange            = 20000;//20080;

  char *LUT_file_name         = "lhc11dv3en";//"lhc11dv2en";
  char *file_type             = "tcs"; // the dot is set below

  ULong64_t NoBins            = Bins;

  ULong64_t scaleq0           = 429496730;//(NoBins << 32) / ScaleRange;
  //if(((scaleq0*ScaleRange)>>32) != NoBins) scaleq0 = scaleq0 + 1;
  ULong64_t scaleq1           = 0;

  Int_t addrPID[Bins]         = {0};
  Int_t addr[Register]        = {0}; // 4 PID values in each register

  Int_t Kink                  = 2988;
  Double_t ChargeDep          = 0.;
  Int_t Iter                  = 0;
  Double_t IterCD             = (Double_t)ScaleRange/(Double_t)Bins;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Fill variables
  ///_______________________________________
  TF1* elLUT = new TF1("elLUT","pol1",0,Kink);
  elLUT->SetParameters(-101.892,0.0889764);

  TF1* nuclLUT = new TF1("nuclLUT","pol1",Kink,ScaleRange);
  nuclLUT->SetParameters(1.06657899999999998e+02,0.0193284);


  // fill each addrPID with values from 0 to 255
  for(Int_t i = 0; i < Bins; i++) {
    if(ChargeDep <= Kink)     addrPID[i] = elLUT  ->Eval(ChargeDep);
    else if(ChargeDep > Kink) addrPID[i] = nuclLUT->Eval(ChargeDep);

    if(addrPID[i] < 0)   addrPID[i] = 0;
    if(addrPID[i] > 255) addrPID[i] = 255;

    ChargeDep = ChargeDep + IterCD;
  }

  // put always 4 addrPID (8bit) in 1 addr (32bit)
  for(Int_t i = 0; i < Register; i++) {
    for(Int_t j = 0; j < n; j++) {

        addr[i] = (addr[i] << 8) + (addrPID[Iter+3-j] & 0xFF);

    }
    Iter = Iter + 4;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Write and save LUT in pwd in a .txt file
  ///______________
  ofstream LUTfile;
  LUTfile.open(Form("%s.%s",LUT_file_name,file_type));

  LUTfile << setw(8) << "tracklet" << setw(8) << "scaleq0" << setw(14) << Form("%u;\n",scaleq0);
  LUTfile << setw(8) << "tracklet" << setw(8) << "scaleq1" << setw(14) << Form("%u;\n",scaleq1);
  LUTfile << setw(8) << "write   " << setw(8) << "127, "   << setw(14) << "0xc029,  " << setw(13) << "2008;\n";
  LUTfile << setw(8) << "write   " << setw(8) << "127, "   << setw(14) << "0xc02b,  " << setw(13) << "2008;\n";

  for(Int_t i = 5; i <= 506; i++) {
    LUTfile << setw(8) << "write   " << setw(8) << "127, " << setw(14) << Form("0x%x,  ",49408+(i-5)) << setw(13) << Form("%i;\n",addr[i-5]);
  }
  LUTfile.close();


  //___________________
  ofstream LUTDrawfile;
  LUTDrawfile.open(Form("%sDraw.%s",LUT_file_name,file_type));

  LUTDrawfile << Form("%u\n",scaleq0);
  LUTDrawfile << Form("%u\n",scaleq1);
  LUTDrawfile << "2008\n";
  LUTDrawfile << "2008\n";

  for(Int_t i = 5; i <= 506; i++) {
    LUTDrawfile << Form("%i\n",addr[i-5]);
  }
  LUTDrawfile.close();
}
