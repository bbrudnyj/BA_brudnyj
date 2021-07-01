/*--------------------------------------------------------------------------------------------------*\
|                                                                                                    |
|    Macro that makes projection histograms with branches of "nmartin_TestTriggerTRDTree.root"       |
|                                                                                                    |
|    - ERP: Efficiency, Rejection and Purity for several cuts in <Q_s> for                           |
|      Deuteron, Triton, Helium3, Alpha, Helium3 & Alpha, d & t & He3 & Alpha and their Antis        |
|                                                                                                    |
| Author: Benjamin Brudnyj (2016)                                                                    |
\*--------------------------------------------------------------------------------------------------*/
void bbrudnyj_ReadTreeTestTRDTriggerERP() {

  static const Int_t n = 26;               // for purity
  static const Int_t m = 26;               // for efficiency and rejection

  // number of i=n,m equals to PID cut starting point in the graphs
  // i=26 eq.to <Q_s> = 0    i=19 eq.to <Q_s> = 70    i=12 eq.to <Q_s> = 140   i=5  eq.to <Q_s> = 210
  // i=25 eq.to <Q_s> = 10   i=18 eq.to <Q_s> = 80    i=11 eq.to <Q_s> = 150   i=4  eq.to <Q_s> = 220
  // i=24 eq.to <Q_s> = 20   i=17 eq.to <Q_s> = 90    i=10 eq.to <Q_s> = 160   i=3  eq.to <Q_s> = 230
  // i=23 eq.to <Q_s> = 30   i=16 eq.to <Q_s> = 100   i=9  eq.to <Q_s> = 170   i=2  eq.to <Q_s> = 240
  // i=22 eq.to <Q_s> = 40   i=15 eq.to <Q_s> = 110   i=8  eq.to <Q_s> = 180   i=1  eq.to <Q_s> = 250
  // i=21 eq.to <Q_s> = 50   i=14 eq.to <Q_s> = 120   i=7  eq.to <Q_s> = 190
  // i=20 eq.to <Q_s> = 60   i=13 eq.to <Q_s> = 130   i=6  eq.to <Q_s> = 200

  gStyle->SetTitleW(.8f);
  gStyle->SetTitleH(.08f);
  gStyle->SetTitleSize(.045,"xy");
  //gStyle->SetTitleOffset(.8,"x");
  gROOT->ForceStyle();

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Define Variables
  ///__________________
  Int_t Sigma      = 2;                    // 2 or 3 is testing in this analysis

  char *AllOrAnti = "All";                 // which particles shall be analyse
  //char *AllOrAnti = "Anti";

  char *Collisions = "PbPb_11h";
  //char *Collisions = "PbPb_11hStd";
  //char *Collisions = "pPb_13b";
  //char *Collisions = "pPb_13c";
  //char *Collisions = "pp_15f";

  if(Sigma == 2) {
    Double_t ClearCutDeuteronTPC    = 1.5; // Rigidity cuts on TPCnSigma histograms
    Double_t ClearCutDeuteronTPCTOF = 1.7;
    Double_t ClearCutTritonTPC      = 1.8;
    Double_t ClearCutTritonTPCTOF   = 2.;
  }
  else if(Sigma == 3) {
    Double_t ClearCutDeuteronTPC    = 1.4;
    Double_t ClearCutDeuteronTPCTOF = 1.7;
    Double_t ClearCutTritonTPC      = 1.7;
    Double_t ClearCutTritonTPCTOF   = 2.;
  }
  else {
    Double_t ClearCutDeuteronTPC    = 20.;
    Double_t ClearCutDeuteronTPCTOF = 20.;
    Double_t ClearCutTritonTPC      = 20.;
    Double_t ClearCutTritonTPCTOF   = 20.;
  }

  Int_t nArray[n] = {0};
  Int_t k = 0;
  for(Int_t i = n; i > 0; i--) {
    nArray[i-1] = 250 - k*10;
    k++;
  }
  Int_t mArray[m] = {0};
  Int_t l = 0;
  for(Int_t i = m; i > 0; i--) {
    nArray[i-1] = 250 - l*10;
    l++;
  }

  Int_t nopc          = 5;                 // "n"umber "o"f "p"ossible "c"andidates are shown in legends (nopc == 3 or 5)
  Float_t y1          = .0;
  if(nopc == 3) y1    = .57;
  if(nopc == 5) y1    = .45;

  Int_t Deu           = 90;               // beginning of the possible cut candidates in the legends in steps of 10
  Int_t Tri           = 90;
  Int_t He3           = 90;
  Int_t Alp           = 90;
  Int_t He3Alp        = 90;
  Int_t DeuTriHe3Alp  = 90;

  Int_t nDeu          = (Deu          - nArray[26-n]) / 10;
  Int_t nTri          = (Tri          - nArray[26-n]) / 10;
  Int_t nHe3          = (He3          - nArray[26-n]) / 10;
  Int_t nAlp          = (Alp          - nArray[26-n]) / 10;
  Int_t nHe3Alp       = (He3Alp       - nArray[26-n]) / 10;
  Int_t nDeuTriHe3Alp = (DeuTriHe3Alp - nArray[26-n]) / 10;

  Int_t mDeu          = (Deu          - mArray[26-m]) / 10;
  Int_t mTri          = (Tri          - mArray[26-m]) / 10;
  Int_t mHe3          = (He3          - mArray[26-m]) / 10;
  Int_t mAlp          = (Alp          - mArray[26-m]) / 10;
  Int_t mHe3Alp       = (He3Alp       - mArray[26-m]) / 10;
  Int_t mDeuTriHe3Alp = (DeuTriHe3Alp - mArray[26-m]) / 10;


  printf(Form("\nCreating and saving 12 canvases with Sigma = %i for %s particles from %s...\n\n",Sigma,AllOrAnti,Collisions));


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Get data from file /lustre/nyx/alice/users/bbrudnyj/trunk/nmartin_nuclei/AliAnalysisTaskTestTriggerTRD.cxx
  /// Get histograms from /lustre/nyx/alice/users/bbrudnyj/codes/bbrudnyj_ReadTreeTestTRDTrigger2PID.C
  ///                 and                                    .../bbrudnyj_ReadTreeTestTRDTrigger3.C
  ///_________________________________________________________________________________________________
  TFile *cinput1     = TFile::Open(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTrigger2/bbrudnyj_TestTRDTrigger2%s%iSigmaPID.root",Collisions,Collisions,Sigma));
  TList *inlistPID   = (TList*)cinput1->Get(Form("bbrudnyj_TestTRDTrigger2%s%iSigmaPID;1",Collisions,Sigma));
  TList *inlistProj  = (TList*)cinput1->Get(Form("bbrudnyj_TestTRDTrigger2%s%iSigmaProj;1",Collisions,Sigma));

  TFile *cinput2     = TFile::Open(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTrigger3/bbrudnyj_TestTRDTrigger3%s%iSigma.root",Collisions,Collisions,Sigma));
  TList *inlistSigma = (TList*)cinput2->Get(Form("bbrudnyj_TestTRDTrigger3%s%iSigma;1",Collisions,Sigma));

  // PID histograms
    TH2D *fHistTRDpTvPID                                  = (TH2D*)inlistPID->FindObject("fHistTRDpTvPID");

    if(AllOrAnti == "All") {
      TH2D *fHistTRDpTvPIDvDeuteronTPC_clear                = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvDeuteronTPC_clear");
      TH2D *fHistTRDpTvPIDvTritonTPC_clear                  = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvTritonTPC_clear");
      TH2D *fHistTRDpTvPIDvHelium3TPC                       = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvHelium3TPC");
      TH2D *fHistTRDpTvPIDvAlphaTPC                         = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAlphaTPC");
      TH2D *fHistTRDpTvPIDvHelium3AlphaTPC                  = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvHelium3AlphaTPC");
      TH2D *fHistTRDpTvPIDvDeuTriHe3AlpTPC                  = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvDeuTriHe3AlpTPC");

      TH2D *fHistTRDpTvPIDvAllexceptDeuteronTPC_clear       = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptDeuteronTPC_clear");
      TH2D *fHistTRDpTvPIDvAllexceptTritonTPC_clear         = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptTritonTPC_clear");
      TH2D *fHistTRDpTvPIDvAllexceptHelium3TPC              = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptHelium3TPC");
      TH2D *fHistTRDpTvPIDvAllexceptAlphaTPC                = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAlphaTPC");
      TH2D *fHistTRDpTvPIDvAllexceptHelium3AlphaTPC         = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptHelium3AlphaTPC");
      TH2D *fHistTRDpTvPIDvAllexceptDeuTriHe3AlpTPC         = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptDeuTriHe3AlpTPC");

      TH2D *fHistTRDpTvPIDvDeuteronTPCTOF_clear              = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvDeuteronTPCTOF_clear");
      TH2D *fHistTRDpTvPIDvTritonTPCTOF_clear                = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvTritonTPCTOF_clear");
      TH2D *fHistTRDpTvPIDvHelium3TPCTOF                     = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvHelium3TPCTOF");
      TH2D *fHistTRDpTvPIDvAlphaTPCTOF                       = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAlphaTPCTOF");
      TH2D *fHistTRDpTvPIDvHelium3AlphaTPCTOF                = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvHelium3AlphaTPCTOF");
      TH2D *fHistTRDpTvPIDvDeuTriHe3AlpTPCTOF                = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvDeuTriHe3AlpTPCTOF");

      TH2D *fHistTRDpTvPIDvAllexceptDeuteronTPCTOF_clear     = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptDeuteronTPCTOF_clear");
      TH2D *fHistTRDpTvPIDvAllexceptTritonTPCTOF_clear       = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptTritonTPCTOF_clear");
      TH2D *fHistTRDpTvPIDvAllexceptHelium3TPCTOF            = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptHelium3TPCTOF");
      TH2D *fHistTRDpTvPIDvAllexceptAlphaTPCTOF              = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAlphaTPCTOF");
      TH2D *fHistTRDpTvPIDvAllexceptHelium3AlphaTPCTOF       = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptHelium3AlphaTPCTOF");
      TH2D *fHistTRDpTvPIDvAllexceptDeuTriHe3AlpTPCTOF       = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptDeuTriHe3AlpTPCTOF");

    // Projection histograms
      // TPC preselection
      TH1D *fHistTRDPIDDeuteronTPC_clear                     = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPC_clear");
      TH1D *fHistTRDPIDTritonTPC_clear                       = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPC_clear");
      TH1D *fHistTRDPIDHelium3TPC                            = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3TPC");
      TH1D *fHistTRDPIDAlphaTPC                              = (TH1D*)inlistProj->FindObject("fHistTRDPIDAlphaTPC");
      TH1D *fHistTRDPIDHelium3AlphaTPC                       = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3AlphaTPC");
      TH1D *fHistTRDPIDDeuTriHe3AlpTPC                       = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuTriHe3AlpTPC");

      TH1D *fHistTRDPIDAllexceptDeuteronTPC_clear            = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptDeuteronTPC_clear");
      TH1D *fHistTRDPIDAllexceptTritonTPC_clear              = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptTritonTPC_clear");
      TH1D *fHistTRDPIDAllexceptHelium3TPC                   = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptHelium3TPC");
      TH1D *fHistTRDPIDAllexceptAlphaTPC                     = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAlphaTPC");
      TH1D *fHistTRDPIDAllexceptHelium3AlphaTPC              = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptHelium3AlphaTPC");
      TH1D *fHistTRDPIDAllexceptDeuTriHe3AlpTPC              = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptDeuTriHe3AlpTPC");

      TH1D *fHistTRDPIDallTPC                                = (TH1D*)inlistProj->FindObject("fHistTRDPIDallTPC");
      TH1D *fHistTRDPIDallDeuteronTPC                        = (TH1D*)inlistProj->FindObject("fHistTRDPIDallDeuteronTPC");
      TH1D *fHistTRDPIDallTritonTPC                          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallTritonTPC");
      TH1D *fHistTRDPIDallHelium3TPC                         = (TH1D*)inlistProj->FindObject("fHistTRDPIDallHelium3TPC");
      TH1D *fHistTRDPIDallAlphaTPC                           = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAlphaTPC");
      TH1D *fHistTRDPIDallHelium3AlphaTPC                    = (TH1D*)inlistProj->FindObject("fHistTRDPIDallHelium3AlphaTPC");
      TH1D *fHistTRDPIDallDeuTriHe3AlpTPC                    = (TH1D*)inlistProj->FindObject("fHistTRDPIDallDeuTriHe3AlpTPC");

      TH1D *fHistTRDPIDDeuteronTPC_clear_Clone               = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPC_clear_Clone");
      TH1D *fHistTRDPIDDeuteronTPC_clear_Clone2              = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPC_clear_Clone2");
      TH1D *fHistTRDPIDTritonTPC_clear_Clone                 = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPC_clear_Clone");
      TH1D *fHistTRDPIDTritonTPC_clear_Clone2                = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPC_clear_Clone2");
      TH1D *fHistTRDPIDHelium3TPC_Clone                      = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3TPC_Clone");
      TH1D *fHistTRDPIDHelium3TPC_Clone2                     = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3TPC_Clone2");
      TH1D *fHistTRDPIDAlphaTPC_Clone                        = (TH1D*)inlistProj->FindObject("fHistTRDPIDAlphaTPC_Clone");
      TH1D *fHistTRDPIDAlphaTPC_Clone2                       = (TH1D*)inlistProj->FindObject("fHistTRDPIDAlphaTPC_Clone2");
      TH1D *fHistTRDPIDHelium3AlphaTPC_Clone                 = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3AlphaTPC_Clone");
      TH1D *fHistTRDPIDDeuTriHe3AlpTPC_Clone                 = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuTriHe3AlpTPC_Clone");

      // TPC & TOF preselection
      TH1D *fHistTRDPIDDeuteronTPCTOF_clear                  = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPCTOF_clear");
      TH1D *fHistTRDPIDTritonTPCTOF_clear                    = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPCTOF_clear");
      TH1D *fHistTRDPIDHelium3TPCTOF                         = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3TPCTOF");
      TH1D *fHistTRDPIDAlphaTPCTOF                           = (TH1D*)inlistProj->FindObject("fHistTRDPIDAlphaTPCTOF");
      TH1D *fHistTRDPIDHelium3AlphaTPCTOF                    = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3AlphaTPCTOF");
      TH1D *fHistTRDPIDDeuTriHe3AlpTPCTOF                    = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuTriHe3AlpTPCTOF");

      TH1D *fHistTRDPIDAllexceptDeuteronTPCTOF_clear         = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptDeuteronTPCTOF_clear");
      TH1D *fHistTRDPIDAllexceptTritonTPCTOF_clear           = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptTritonTPCTOF_clear");
      TH1D *fHistTRDPIDAllexceptHelium3TPCTOF                = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptHelium3TPCTOF");
      TH1D *fHistTRDPIDAllexceptAlphaTPCTOF                  = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAlphaTPCTOF");
      TH1D *fHistTRDPIDAllexceptHelium3AlphaTPCTOF           = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptHelium3AlphaTPCTOF");
      TH1D *fHistTRDPIDAllexceptDeuTriHe3AlpTPCTOF           = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptDeuTriHe3AlpTPCTOF");

      TH1D *fHistTRDPIDallTPCTOF                             = (TH1D*)inlistProj->FindObject("fHistTRDPIDallTPCTOF");
      TH1D *fHistTRDPIDallDeuteronTPCTOF                     = (TH1D*)inlistProj->FindObject("fHistTRDPIDallDeuteronTPCTOF");
      TH1D *fHistTRDPIDallTritonTPCTOF                       = (TH1D*)inlistProj->FindObject("fHistTRDPIDallTritonTPCTOF");
      TH1D *fHistTRDPIDallHelium3TPCTOF                      = (TH1D*)inlistProj->FindObject("fHistTRDPIDallHelium3TPCTOF");
      TH1D *fHistTRDPIDallAlphaTPCTOF                        = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAlphaTPCTOF");
      TH1D *fHistTRDPIDallHelium3AlphaTPCTOF                 = (TH1D*)inlistProj->FindObject("fHistTRDPIDallHelium3AlphaTPCTOF");
      TH1D *fHistTRDPIDallDeuTriHe3AlpTPCTOF                 = (TH1D*)inlistProj->FindObject("fHistTRDPIDallDeuTriHe3AlpTPCTOF");

      TH1D *fHistTRDPIDDeuteronTPCTOF_clear_Clone            = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPCTOF_clear_Clone");
      TH1D *fHistTRDPIDDeuteronTPCTOF_clear_Clone2           = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPCTOF_clear_Clone2");
      TH1D *fHistTRDPIDTritonTPCTOF_clear_Clone              = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPCTOF_clear_Clone");
      TH1D *fHistTRDPIDTritonTPCTOF_clear_Clone2             = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPCTOF_clear_Clone2");
      TH1D *fHistTRDPIDHelium3TPCTOF_Clone                   = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3TPCTOF_Clone");
      TH1D *fHistTRDPIDHelium3TPCTOF_Clone2                  = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3TPCTOF_Clone2");
      TH1D *fHistTRDPIDAlphaTPCTOF_Clone                     = (TH1D*)inlistProj->FindObject("fHistTRDPIDAlphaTPCTOF_Clone");
      TH1D *fHistTRDPIDAlphaTPCTOF_Clone2                    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAlphaTPCTOF_Clone2");
      TH1D *fHistTRDPIDHelium3AlphaTPCTOF_Clone              = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3AlphaTPCTOF_Clone");
      TH1D *fHistTRDPIDDeuTriHe3AlpTPCTOF_Clone              = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuTriHe3AlpTPCTOF_Clone");

    // nSigma histograms
      TH2D *fHistTPCnSigmavRigvDeuteronTPC_clear             = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvDeuteronTPC_clear");
      TH2D *fHistTPCnSigmavRigvTritonTPC_clear               = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvTritonTPC_clear");
      TH2D *fHistTPCnSigmavRigvHelium3TPC                    = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvHelium3TPC");
      TH2D *fHistTPCnSigmavRigvAlphaTPC                      = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvAlphaTPC");
      //TH2D *fHistTPCnSigmavRigvHelium3AlphaTPC               = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvHelium3AlphaTPC");
      //TH2D *fHistTPCnSigmavRigvDeuTriHe3AlpTPC               = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvDeuTriHe3AlpTPC");

      //TH2D *fHistTPCnSigmavRigvDeuteronTPC_clear             = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvDeuteronTPC_clear");
      //TH2D *fHistTPCnSigmavRigvTritonTPC_clear               = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvTritonTPC_clear");
      //TH2D *fHistTPCnSigmavRigvHelium3TPC                    = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvHelium3TPC");
      //TH2D *fHistTPCnSigmavRigvAlphaTPC                      = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvAlphaTPC");
      //TH2D *fHistTPCnSigmavRigvHelium3AlphaTPC               = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvHelium3AlphaTPC");
      //TH2D *fHistTPCnSigmavRigvDeuTriHe3AlpTPC               = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvDeuTriHe3AlpTPC");

      TH2D *fHistTPCnSigmavRigvDeuteronTPCTOF_clear          = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvDeuteronTPCTOF_clear");
      TH2D *fHistTPCnSigmavRigvTritonTPCTOF_clear            = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvTritonTPCTOF_clear");
      TH2D *fHistTPCnSigmavRigvHelium3TPCTOF                 = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvHelium3TPCTOF");
      TH2D *fHistTPCnSigmavRigvAlphaTPCTOF                   = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvAlphaTPCTOF");
      //TH2D *fHistTPCnSigmavRigvHelium3AlphaTPCTOF            = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvHelium3AlphaTPCTOF");
      //TH2D *fHistTPCnSigmavRigvDeuTriHe3AlpTPCTOF            = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvDeuTriHe3AlpTPCTOF");

      //TH2D *fHistTPCnSigmavRigvDeuteronTPCTOF_clear          = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvDeuteronTPCTOF_clear");
      //TH2D *fHistTPCnSigmavRigvTritonTPCTOF_clear            = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvTritonTPCTOF_clear");
      //TH2D *fHistTPCnSigmavRigvHelium3TPCTOF                 = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvHelium3TPCTOF");
      //TH2D *fHistTPCnSigmavRigvAlphaTPCTOF                   = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvAlphaTPCTOF");
      //TH2D *fHistTPCnSigmavRigvHelium3AlphaTPCTOF            = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvHelium3AlphaTPCTOF");
      //TH2D *fHistTPCnSigmavRigvDeuTriHe3AlpTPCTOF            = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvDeuTriHe3AlpTPCTOF");
    }



    if(AllOrAnti == "Anti") {
      TH2D *fHistTRDpTvPIDvDeuteronTPC_clear                = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiDeuteronTPC_clear");
      TH2D *fHistTRDpTvPIDvTritonTPC_clear                  = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiTritonTPC_clear");
      TH2D *fHistTRDpTvPIDvHelium3TPC                       = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiHelium3TPC");
      TH2D *fHistTRDpTvPIDvAlphaTPC                         = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiAlphaTPC");
      TH2D *fHistTRDpTvPIDvHelium3AlphaTPC                  = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiHelium3AlphaTPC");
      TH2D *fHistTRDpTvPIDvDeuTriHe3AlpTPC                  = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiDeuTriHe3AlpTPC");

      TH2D *fHistTRDpTvPIDvAllexceptDeuteronTPC_clear       = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiDeuteronTPC_clear");
      TH2D *fHistTRDpTvPIDvAllexceptTritonTPC_clear         = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiTritonTPC_clear");
      TH2D *fHistTRDpTvPIDvAllexceptHelium3TPC              = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiHelium3TPC");
      TH2D *fHistTRDpTvPIDvAllexceptAlphaTPC                = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiAlphaTPC");
      TH2D *fHistTRDpTvPIDvAllexceptHelium3AlphaTPC         = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiHelium3AlphaTPC");
      TH2D *fHistTRDpTvPIDvAllexceptDeuTriHe3AlpTPC         = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiDeuTriHe3AlpTPC");

      TH2D *fHistTRDpTvPIDvDeuteronTPCTOF_clear              = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiDeuteronTPCTOF_clear");
      TH2D *fHistTRDpTvPIDvTritonTPCTOF_clear                = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiTritonTPCTOF_clear");
      TH2D *fHistTRDpTvPIDvHelium3TPCTOF                     = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiHelium3TPCTOF");
      TH2D *fHistTRDpTvPIDvAlphaTPCTOF                       = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiAlphaTPCTOF");
      TH2D *fHistTRDpTvPIDvHelium3AlphaTPCTOF                = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiHelium3AlphaTPCTOF");
      TH2D *fHistTRDpTvPIDvDeuTriHe3AlpTPCTOF                = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiDeuTriHe3AlpTPCTOF");

      TH2D *fHistTRDpTvPIDvAllexceptDeuteronTPCTOF_clear     = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiDeuteronTPCTOF_clear");
      TH2D *fHistTRDpTvPIDvAllexceptTritonTPCTOF_clear       = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiTritonTPCTOF_clear");
      TH2D *fHistTRDpTvPIDvAllexceptHelium3TPCTOF            = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiHelium3TPCTOF");
      TH2D *fHistTRDpTvPIDvAllexceptAlphaTPCTOF              = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiAlphaTPCTOF");
      TH2D *fHistTRDpTvPIDvAllexceptHelium3AlphaTPCTOF       = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiHelium3AlphaTPCTOF");
      TH2D *fHistTRDpTvPIDvAllexceptDeuTriHe3AlpTPCTOF       = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiDeuTriHe3AlpTPCTOF");

    // Projection histograms
      // TPC preselection
      TH1D *fHistTRDPIDDeuteronTPC_clear                     = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuteronTPC_clear");
      TH1D *fHistTRDPIDTritonTPC_clear                       = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiTritonTPC_clear");
      TH1D *fHistTRDPIDHelium3TPC                            = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3TPC");
      TH1D *fHistTRDPIDAlphaTPC                              = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiAlphaTPC");
      TH1D *fHistTRDPIDHelium3AlphaTPC                       = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3AlphaTPC");
      TH1D *fHistTRDPIDDeuTriHe3AlpTPC                       = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuTriHe3AlpTPC");

      TH1D *fHistTRDPIDAllexceptDeuteronTPC_clear            = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiDeuteronTPC_clear");
      TH1D *fHistTRDPIDAllexceptTritonTPC_clear              = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiTritonTPC_clear");
      TH1D *fHistTRDPIDAllexceptHelium3TPC                   = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiHelium3TPC");
      TH1D *fHistTRDPIDAllexceptAlphaTPC                     = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiAlphaTPC");
      TH1D *fHistTRDPIDAllexceptHelium3AlphaTPC              = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiHelium3AlphaTPC");
      TH1D *fHistTRDPIDAllexceptDeuTriHe3AlpTPC              = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiDeuTriHe3AlpTPC");

      TH1D *fHistTRDPIDallDeuteronTPC                        = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiDeuteronTPC");
      TH1D *fHistTRDPIDallTritonTPC                          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiTritonTPC");
      TH1D *fHistTRDPIDallHelium3TPC                         = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiHelium3TPC");
      TH1D *fHistTRDPIDallAlphaTPC                           = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiAlphaTPC");
      TH1D *fHistTRDPIDallHelium3AlphaTPC                    = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiHelium3AlphaTPC");
      TH1D *fHistTRDPIDallDeuTriHe3AlpTPC                    = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiDeuTriHe3AlpTPC");

      TH1D *fHistTRDPIDDeuteronTPC_clear_Clone               = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuteronTPC_clear_Clone");
      TH1D *fHistTRDPIDDeuteronTPC_clear_Clone2              = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuteronTPC_clear_Clone2");
      TH1D *fHistTRDPIDTritonTPC_clear_Clone                 = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiTritonTPC_clear_Clone");
      TH1D *fHistTRDPIDTritonTPC_clear_Clone2                = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiTritonTPC_clear_Clone2");
      TH1D *fHistTRDPIDHelium3TPC_Clone                      = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3TPC_Clone");
      TH1D *fHistTRDPIDHelium3TPC_Clone2                     = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3TPC_Clone2");
      TH1D *fHistTRDPIDAlphaTPC_Clone                        = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiAlphaTPC_Clone");
      TH1D *fHistTRDPIDAlphaTPC_Clone2                       = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiAlphaTPC_Clone2");
      TH1D *fHistTRDPIDHelium3AlphaTPC_Clone                 = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3AlphaTPC_Clone");
      TH1D *fHistTRDPIDDeuTriHe3AlpTPC_Clone                 = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuTriHe3AlpTPC_Clone");

      // TPC & TOF preselection
      TH1D *fHistTRDPIDDeuteronTPCTOF_clear                  = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuteronTPCTOF_clear");
      TH1D *fHistTRDPIDTritonTPCTOF_clear                    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiTritonTPCTOF_clear");
      TH1D *fHistTRDPIDHelium3TPCTOF                         = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3TPCTOF");
      TH1D *fHistTRDPIDAlphaTPCTOF                           = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiAlphaTPCTOF");
      TH1D *fHistTRDPIDHelium3AlphaTPCTOF                    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3AlphaTPCTOF");
      TH1D *fHistTRDPIDDeuTriHe3AlpTPCTOF                    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuTriHe3AlpTPCTOF");

      TH1D *fHistTRDPIDAllexceptDeuteronTPCTOF_clear         = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiDeuteronTPCTOF_clear");
      TH1D *fHistTRDPIDAllexceptTritonTPCTOF_clear           = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiTritonTPCTOF_clear");
      TH1D *fHistTRDPIDAllexceptHelium3TPCTOF                = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiHelium3TPCTOF");
      TH1D *fHistTRDPIDAllexceptAlphaTPCTOF                  = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiAlphaTPCTOF");
      TH1D *fHistTRDPIDAllexceptHelium3AlphaTPCTOF           = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiHelium3AlphaTPCTOF");
      TH1D *fHistTRDPIDAllexceptDeuTriHe3AlpTPCTOF           = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiDeuTriHe3AlpTPCTOF");

      TH1D *fHistTRDPIDallDeuteronTPCTOF                     = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiDeuteronTPCTOF");
      TH1D *fHistTRDPIDallTritonTPCTOF                       = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiTritonTPCTOF");
      TH1D *fHistTRDPIDallHelium3TPCTOF                      = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiHelium3TPCTOF");
      TH1D *fHistTRDPIDallAlphaTPCTOF                        = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiAlphaTPCTOF");
      TH1D *fHistTRDPIDallHelium3AlphaTPCTOF                 = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiHelium3AlphaTPCTOF");
      TH1D *fHistTRDPIDallDeuTriHe3AlpTPCTOF                 = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiDeuTriHe3AlpTPCTOF");

      TH1D *fHistTRDPIDDeuteronTPCTOF_clear_Clone            = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuteronTPCTOF_clear_Clone");
      TH1D *fHistTRDPIDDeuteronTPCTOF_clear_Clone2           = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuteronTPCTOF_clear_Clone2");
      TH1D *fHistTRDPIDTritonTPCTOF_clear_Clone              = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiTritonTPCTOF_clear_Clone");
      TH1D *fHistTRDPIDTritonTPCTOF_clear_Clone2             = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiTritonTPCTOF_clear_Clone2");
      TH1D *fHistTRDPIDHelium3TPCTOF_Clone                   = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3TPCTOF_Clone");
      TH1D *fHistTRDPIDHelium3TPCTOF_Clone2                  = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3TPCTOF_Clone2");
      TH1D *fHistTRDPIDAlphaTPCTOF_Clone                     = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiAlphaTPCTOF_Clone");
      TH1D *fHistTRDPIDAlphaTPCTOF_Clone2                    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiAlphaTPCTOF_Clone2");
      TH1D *fHistTRDPIDHelium3AlphaTPCTOF_Clone              = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3AlphaTPCTOF_Clone");
      TH1D *fHistTRDPIDDeuTriHe3AlpTPCTOF_Clone              = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuTriHe3AlpTPCTOF_Clone");

    // nSigma histograms
      //TH2D *fHistTPCnSigmavRigvDeuteronTPC_clear             = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvAntiDeuteronTPC_clear"); // have to be also created in nicoles code
      //TH2D *fHistTPCnSigmavRigvTritonTPC_clear               = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvAntiTritonTPC_clear");
      //TH2D *fHistTPCnSigmavRigvHelium3TPC                    = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvAntiHelium3TPC");
      //TH2D *fHistTPCnSigmavRigvAlphaTPC                      = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvAntiAlphaTPC");
      //TH2D *fHistTPCnSigmavRigvHelium3AlphaTPC               = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvAntiHelium3AlphaTPC");
      //TH2D *fHistTPCnSigmavRigvDeuTriHe3AlpTPC               = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvAntiDeuTriHe3AlpTPC");

      //TH2D *fHistTPCnSigmavRigvDeuteronTPCTOF_clear          = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvAntiDeuteronTPCTOF_clear");
      //TH2D *fHistTPCnSigmavRigvTritonTPCTOF_clear            = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvAntiTritonTPCTOF_clear");
      //TH2D *fHistTPCnSigmavRigvHelium3TPCTOF                 = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvAntiHelium3TPCTOF");
      //TH2D *fHistTPCnSigmavRigvAlphaTPCTOF                   = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvAntiAlphaTPCTOF");
      //TH2D *fHistTPCnSigmavRigvHelium3AlphaTPCTOF            = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvAntiHelium3AlphaTPCTOF");
      //TH2D *fHistTPCnSigmavRigvDeuTriHe3AlpTPCTOF            = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvAntiDeuTriHe3AlpTPCTOF");
    }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Create graphs
  ///______
  // Purity
    Double_t nMinBin                      = 260. - n*10;
    Double_t nMeanQRange[n]               = {0};

    Int_t nPIDall[n]                      = {0};

    Int_t nDeuteronTPC[n]                 = {0};
    Int_t nDeuteronTPCTOF[n]              = {0};
    Int_t nTritonTPC[n]                   = {0};
    Int_t nTritonTPCTOF[n]                = {0};
    Int_t nHelium3TPC[n]                  = {0};
    Int_t nHelium3TPCTOF[n]               = {0};
    Int_t nAlphaTPC[n]                    = {0};
    Int_t nAlphaTPCTOF[n]                 = {0};
    Int_t nHelium3AlphaTPC[n]             = {0};
    Int_t nHelium3AlphaTPCTOF[n]          = {0};
    Int_t nDeuTriHe3AlpTPC[n]             = {0};
    Int_t nDeuTriHe3AlpTPCTOF[n]          = {0};

    Double_t PurityDeuteronTPC[n]         = {0};
    Double_t PurityDeuteronTPCTOF[n]      = {0};
    Double_t PurityTritonTPC[n]           = {0};
    Double_t PurityTritonTPCTOF[n]        = {0};
    Double_t PurityHelium3TPC[n]          = {0};
    Double_t PurityHelium3TPCTOF[n]       = {0};
    Double_t PurityAlphaTPC[n]            = {0};
    Double_t PurityAlphaTPCTOF[n]         = {0};
    Double_t PurityHelium3AlphaTPC[n]     = {0};
    Double_t PurityHelium3AlphaTPCTOF[n]  = {0};
    Double_t PurityDeuTriHe3AlpTPC[n]     = {0};
    Double_t PurityDeuTriHe3AlpTPCTOF[n]  = {0};

    Double_t ex[n]                        = {0};

    Double_t ePIDall[n]                   = {0};

    Double_t eDeuteronTPC[n]              = {0};
    Double_t eDeuteronTPCTOF[n]           = {0};
    Double_t eTritonTPC[n]                = {0};
    Double_t eTritonTPCTOF[n]             = {0};
    Double_t eHelium3TPC[n]               = {0};
    Double_t eHelium3TPCTOF[n]            = {0};
    Double_t eAlphaTPC[n]                 = {0};
    Double_t eAlphaTPCTOF[n]              = {0};
    Double_t eHelium3AlphaTPC[n]          = {0};
    Double_t eHelium3AlphaTPCTOF[n]       = {0};
    Double_t eDeuTriHe3AlpTPC[n]          = {0};
    Double_t eDeuTriHe3AlpTPCTOF[n]       = {0};

    Double_t ePurityDeuteronTPC[n]        = {0};
    Double_t ePurityDeuteronTPCTOF[n]     = {0};
    Double_t ePurityTritonTPC[n]          = {0};
    Double_t ePurityTritonTPCTOF[n]       = {0};
    Double_t ePurityHelium3TPC[n]         = {0};
    Double_t ePurityHelium3TPCTOF[n]      = {0};
    Double_t ePurityAlphaTPC[n]           = {0};
    Double_t ePurityAlphaTPCTOF[n]        = {0};
    Double_t ePurityHelium3AlphaTPC[n]    = {0};
    Double_t ePurityHelium3AlphaTPCTOF[n] = {0};
    Double_t ePurityDeuTriHe3AlpTPC[n]    = {0};
    Double_t ePurityDeuTriHe3AlpTPCTOF[n] = {0};


    for(Int_t i = 0; i < n; i++) {
      nMeanQRange[i]               = nMinBin + i*10;

      nPIDall[i]                   = fHistTRDpTvPID                     ->IntegralAndError(1,500,nMeanQRange[i],260,ePIDall[i]);

      nDeuteronTPC[i]              = fHistTRDpTvPIDvDeuteronTPC_clear   ->IntegralAndError(1,160,nMeanQRange[i],260,eDeuteronTPC[i]);
      nDeuteronTPCTOF[i]           = fHistTRDpTvPIDvDeuteronTPCTOF_clear->IntegralAndError(1,160,nMeanQRange[i],260,eDeuteronTPCTOF[i]);
      nTritonTPC[i]                = fHistTRDpTvPIDvTritonTPC_clear     ->IntegralAndError(1,160,nMeanQRange[i],260,eTritonTPC[i]);
      nTritonTPCTOF[i]             = fHistTRDpTvPIDvTritonTPCTOF_clear  ->IntegralAndError(1,160,nMeanQRange[i],260,eTritonTPCTOF[i]);
      nHelium3TPC[i]               = fHistTRDpTvPIDvHelium3TPC          ->IntegralAndError(1,160,nMeanQRange[i],260,eHelium3TPC[i]);
      nHelium3TPCTOF[i]            = fHistTRDpTvPIDvHelium3TPCTOF       ->IntegralAndError(1,160,nMeanQRange[i],260,eHelium3TPCTOF[i]);
      nAlphaTPC[i]                 = fHistTRDpTvPIDvAlphaTPC            ->IntegralAndError(1,160,nMeanQRange[i],260,eAlphaTPC[i]);
      if(AllOrAnti == "All")nAlphaTPCTOF[i] = fHistTRDpTvPIDvAlphaTPCTOF->IntegralAndError(1,160,nMeanQRange[i],260,eAlphaTPCTOF[i]);
      nHelium3AlphaTPC[i]          = fHistTRDpTvPIDvHelium3AlphaTPC     ->IntegralAndError(1,160,nMeanQRange[i],260,eHelium3AlphaTPC[i]);
      nHelium3AlphaTPCTOF[i]       = fHistTRDpTvPIDvHelium3AlphaTPCTOF  ->IntegralAndError(1,160,nMeanQRange[i],260,eHelium3AlphaTPCTOF[i]);
      nDeuTriHe3AlpTPC[i]          = fHistTRDpTvPIDvDeuTriHe3AlpTPC     ->IntegralAndError(1,160,nMeanQRange[i],260,eDeuTriHe3AlpTPC[i]);
      nDeuTriHe3AlpTPCTOF[i]       = fHistTRDpTvPIDvDeuTriHe3AlpTPCTOF  ->IntegralAndError(1,160,nMeanQRange[i],260,eDeuTriHe3AlpTPCTOF[i]);

      PurityDeuteronTPC[i]         = (Double_t)nDeuteronTPC[i]        / nPIDall[i];
      PurityDeuteronTPCTOF[i]      = (Double_t)nDeuteronTPCTOF[i]     / nPIDall[i];
      PurityTritonTPC[i]           = (Double_t)nTritonTPC[i]          / nPIDall[i];
      PurityTritonTPCTOF[i]        = (Double_t)nTritonTPCTOF[i]       / nPIDall[i];
      PurityHelium3TPC[i]          = (Double_t)nHelium3TPC[i]         / nPIDall[i];
      PurityHelium3TPCTOF[i]       = (Double_t)nHelium3TPCTOF[i]      / nPIDall[i];
      PurityAlphaTPC[i]            = (Double_t)nAlphaTPC[i]           / nPIDall[i];
      if(AllOrAnti == "All")PurityAlphaTPCTOF[i] = (Double_t)nAlphaTPCTOF[i] / nPIDall[i];
      PurityHelium3AlphaTPC[i]     = (Double_t)nHelium3AlphaTPC[i]    / nPIDall[i];
      PurityHelium3AlphaTPCTOF[i]  = (Double_t)nHelium3AlphaTPCTOF[i] / nPIDall[i];
      PurityDeuTriHe3AlpTPC[i]     = (Double_t)nDeuTriHe3AlpTPC[i]    / nPIDall[i];
      PurityDeuTriHe3AlpTPCTOF[i]  = (Double_t)nDeuTriHe3AlpTPCTOF[i] / nPIDall[i];

      ePurityDeuteronTPC[i]        = GaussError(nDeuteronTPC[i]       ,eDeuteronTPC[i]       ,nPIDall[i],ePIDall[i]); //own function, see below
      ePurityDeuteronTPCTOF[i]     = GaussError(nDeuteronTPCTOF[i]    ,eDeuteronTPCTOF[i]    ,nPIDall[i],ePIDall[i]);
      ePurityTritonTPC[i]          = GaussError(nTritonTPC[i]         ,eTritonTPC[i]         ,nPIDall[i],ePIDall[i]);
      ePurityTritonTPCTOF[i]       = GaussError(nTritonTPCTOF[i]      ,eTritonTPCTOF[i]      ,nPIDall[i],ePIDall[i]);
      ePurityHelium3TPC[i]         = GaussError(nHelium3TPC[i]        ,eHelium3TPC[i]        ,nPIDall[i],ePIDall[i]);
      ePurityHelium3TPCTOF[i]      = GaussError(nHelium3TPCTOF[i]     ,eHelium3TPCTOF[i]     ,nPIDall[i],ePIDall[i]);
      ePurityAlphaTPC[i]           = GaussError(nAlphaTPC[i]          ,eAlphaTPC[i]          ,nPIDall[i],ePIDall[i]);
      if(AllOrAnti == "All")ePurityAlphaTPCTOF[i] = GaussError(nAlphaTPCTOF[i],eAlphaTPCTOF[i],nPIDall[i],ePIDall[i]);
      ePurityHelium3AlphaTPC[i]    = GaussError(nHelium3AlphaTPC[i]   ,eHelium3AlphaTPC[i]   ,nPIDall[i],ePIDall[i]);
      ePurityHelium3AlphaTPCTOF[i] = GaussError(nHelium3AlphaTPCTOF[i],eHelium3AlphaTPCTOF[i],nPIDall[i],ePIDall[i]);
      ePurityDeuTriHe3AlpTPC[i]    = GaussError(nDeuTriHe3AlpTPC[i]   ,eDeuTriHe3AlpTPC[i]   ,nPIDall[i],ePIDall[i]);
      ePurityDeuTriHe3AlpTPCTOF[i] = GaussError(nDeuTriHe3AlpTPCTOF[i],eDeuTriHe3AlpTPCTOF[i],nPIDall[i],ePIDall[i]);

      ex[i]                        = 0.;

      //cout<<nMeanQRange[i]<<"\t"<<nDeuteronTPC[i]       <<"\t"<<eDeuteronTPC[i]       <<"\t"<<nPIDall[i]<<"\t"<<ePIDall[i]<<"\t("<<PurityDeuteronTPC[i]       <<"\t+- "<<ePurityDeuteronTPC[i]    <<")"<<endl;
      //cout<<nMeanQRange[i]<<"\t"<<nDeuteronTPCTOF[i]    <<"\t"<<eDeuteronTPCTOF[i]    <<"\t"<<nPIDall[i]<<"\t"<<ePIDall[i]<<"\t("<<PurityDeuteronTPCTOF[i]    <<"\t+- "<<ePurityDeuteronTPCTOF[i] <<")"<<endl;
      //cout<<nMeanQRange[i]<<"\t"<<nTritonTPC[i]         <<"\t"<<eTritonTPC[i]         <<"\t"<<nPIDall[i]<<"\t"<<ePIDall[i]<<"\t("<<PurityTritonTPC[i]         <<"\t+- "<<ePurityTritonTPC[i]      <<")"<<endl;
      //cout<<nMeanQRange[i]<<"\t"<<nTritonTPCTOF[i]      <<"\t"<<eTritonTPCTOF[i]      <<"\t"<<nPIDall[i]<<"\t"<<ePIDall[i]<<"\t("<<PurityTritonTPCTOF[i]      <<"\t+- "<<ePurityTritonTPCTOF[i]   <<")"<<endl;
      //cout<<nMeanQRange[i]<<"\t"<<nHelium3TPC[i]        <<"\t"<<eHelium3TPC[i]        <<"\t"<<nPIDall[i]<<"\t"<<ePIDall[i]<<"\t("<<PurityHelium3TPC[i]        <<"\t+- "<<ePurityHelium3TPC[i]     <<")"<<endl;
      //cout<<nMeanQRange[i]<<"\t"<<nHelium3TPCTOF[i]     <<"\t"<<eHelium3TPCTOF[i]     <<"\t"<<nPIDall[i]<<"\t"<<ePIDall[i]<<"\t("<<PurityHelium3TPCTOF[i]     <<"\t+- "<<ePurityHelium3TPCTOF[i]  <<")"<<endl;
      //cout<<nMeanQRange[i]<<"\t"<<nAlphaTPC[i]          <<"\t"<<eAlphaTPC[i]          <<"\t"<<nPIDall[i]<<"\t"<<ePIDall[i]<<"\t("<<PurityAlphaTPC[i]          <<"\t+- "<<ePurityAlphaTPC[i]       <<")"<<endl;
      //cout<<nMeanQRange[i]<<"\t"<<nAlphaTPCTOF[i]       <<"\t"<<eAlphaTPCTOF[i]       <<"\t"<<nPIDall[i]<<"\t"<<ePIDall[i]<<"\t("<<PurityAlphaTPCTOF[i]       <<"\t+- "<<ePurityAlphaTPCTOF[i]    <<")"<<endl;
      //cout<<nMeanQRange[i]<<"\t"<<nHelium3AlphaTPC[i]   <<"\t"<<eHelium3AlphaTPC[i]   <<"\t"<<nPIDall[i]<<"\t"<<ePIDall[i]<<"\t("<<PurityHelium3AlphaTPC[i]   <<"\t+- "<<ePurityHelium3AlphaTPC[i]<<")"<<endl;
      //cout<<nMeanQRange[i]<<"\t"<<nHelium3AlphaTPCTOF[i]<<"\t"<<eHelium3AlphaTPCTOF[i]<<"\t"<<nPIDall[i]<<"\t"<<ePIDall[i]<<"\t("<<PurityHelium3AlphaTPCTOF[i]<<"\t+- "<<ePurityHelium3AlphaTPCTOF[i]<<")"<<endl;
      //cout<<nMeanQRange[i]<<"\t"<<nDeuTriHe3AlpTPC[i]   <<"\t"<<eDeuTriHe3AlpTPC[i]   <<"\t"<<nPIDall[i]<<"\t"<<ePIDall[i]<<"\t("<<PurityDeuTriHe3AlpTPC[i]   <<"\t+- "<<ePurityDeuTriHe3AlpTPC[i]<<")"<<endl;
      //cout<<nMeanQRange[i]<<"\t"<<nDeuTriHe3AlpTPCTOF[i]<<"\t"<<eDeuTriHe3AlpTPCTOF[i]<<"\t"<<nPIDall[i]<<"\t"<<ePIDall[i]<<"\t("<<PurityDeuTriHe3AlpTPCTOF[i]<<"\t+- "<<ePurityDeuTriHe3AlpTPCTOF[i]<<")"<<endl;

    }

    TGraphErrors *fPurityDeuteronTPC = new TGraphErrors(n,nMeanQRange,PurityDeuteronTPC,ex,ePurityDeuteronTPC);
    if(AllOrAnti == "All") fPurityDeuteronTPC->SetTitle("Deuteron Purity (TPC ps);<Q_{s}> (a.u.);Purity");
    if(AllOrAnti == "Anti")fPurityDeuteronTPC->SetTitle("Anti-Deuteron Purity (TPC ps);<Q_{s}> (a.u.);Purity");
    fPurityDeuteronTPC->SetMarkerStyle(7);
    fPurityDeuteronTPC->SetMarkerColor(4);
    TGraphErrors *fPurityDeuteronTPCTOF = new TGraphErrors(n,nMeanQRange,PurityDeuteronTPCTOF,ex,ePurityDeuteronTPCTOF);
    if(AllOrAnti == "All") fPurityDeuteronTPCTOF->SetTitle("Deuteron Purity (TPC&TOF ps);<Q_{s}> (a.u.);Purity");
    if(AllOrAnti == "Anti")fPurityDeuteronTPCTOF->SetTitle("Anti-Deuteron Purity (TPC&TOF ps);<Q_{s}> (a.u.);Purity");
    fPurityDeuteronTPCTOF->SetMarkerStyle(7);
    fPurityDeuteronTPCTOF->SetMarkerColor(4);

    TGraphErrors *fPurityTritonTPC = new TGraphErrors(n,nMeanQRange,PurityTritonTPC,ex,ePurityTritonTPC);
    if(AllOrAnti == "All") fPurityTritonTPC->SetTitle("Triton Purity (TPC ps);<Q_{s}> (a.u.);Purity");
    if(AllOrAnti == "Anti")fPurityTritonTPC->SetTitle("Anti-Triton Purity (TPC ps);<Q_{s}> (a.u.);Purity");
    fPurityTritonTPC->SetMarkerStyle(7);
    fPurityTritonTPC->SetMarkerColor(4);
    TGraphErrors *fPurityTritonTPCTOF = new TGraphErrors(n,nMeanQRange,PurityTritonTPCTOF,ex,ePurityTritonTPCTOF);
    if(AllOrAnti == "All") fPurityTritonTPCTOF->SetTitle("Triton Purity (TPC&TOF ps);<Q_{s}> (a.u.);Purity");
    if(AllOrAnti == "Anti")fPurityTritonTPCTOF->SetTitle("Anti-Triton Purity (TPC&TOF ps);<Q_{s}> (a.u.);Purity");
    fPurityTritonTPCTOF->SetMarkerStyle(7);
    fPurityTritonTPCTOF->SetMarkerColor(4);

    TGraphErrors *fPurityHelium3TPC = new TGraphErrors(n,nMeanQRange,PurityHelium3TPC,ex,ePurityHelium3TPC);
    if(AllOrAnti == "All") fPurityHelium3TPC->SetTitle("Helium3 Purity (TPC ps);<Q_{s}> (a.u.);Purity");
    if(AllOrAnti == "Anti")fPurityHelium3TPC->SetTitle("Anti-Helium3 Purity (TPC ps);<Q_{s}> (a.u.);Purity");
    fPurityHelium3TPC->SetMarkerStyle(7);
    fPurityHelium3TPC->SetMarkerColor(4);
    TGraphErrors *fPurityHelium3TPCTOF = new TGraphErrors(n,nMeanQRange,PurityHelium3TPCTOF,ex,ePurityHelium3TPCTOF);
    if(AllOrAnti == "All") fPurityHelium3TPCTOF->SetTitle("Helium3 Purity (TPC&TOF ps);<Q_{s}> (a.u.);Purity");
    if(AllOrAnti == "Anti")fPurityHelium3TPCTOF->SetTitle("Anti-Helium3 Purity (TPC&TOF ps);<Q_{s}> (a.u.);Purity");
    fPurityHelium3TPCTOF->SetMarkerStyle(7);
    fPurityHelium3TPCTOF->SetMarkerColor(4);

    TGraphErrors *fPurityAlphaTPC = new TGraphErrors(n,nMeanQRange,PurityAlphaTPC,ex,ePurityAlphaTPC);
    if(AllOrAnti == "All") fPurityAlphaTPC->SetTitle("Alpha Purity (TPC ps);<Q_{s}> (a.u.);Purity");
    if(AllOrAnti == "Anti")fPurityAlphaTPC->SetTitle("Anti-Alpha Purity (TPC ps);<Q_{s}> (a.u.);Purity");
    fPurityAlphaTPC->SetMarkerStyle(7);
    fPurityAlphaTPC->SetMarkerColor(4);
    TGraphErrors *fPurityAlphaTPCTOF = new TGraphErrors(n,nMeanQRange,PurityAlphaTPCTOF,ex,ePurityAlphaTPCTOF);
    if(AllOrAnti == "All") fPurityAlphaTPCTOF->SetTitle("Alpha Purity (TPC&TOF ps);<Q_{s}> (a.u.);Purity");
    if(AllOrAnti == "Anti")fPurityAlphaTPCTOF->SetTitle("Anti-Alpha Purity (TPC&TOF ps);<Q_{s}> (a.u.);Purity");
    fPurityAlphaTPCTOF->SetMarkerStyle(7);
    fPurityAlphaTPCTOF->SetMarkerColor(4);

    TGraphErrors *fPurityHelium3AlphaTPC = new TGraphErrors(n,nMeanQRange,PurityHelium3AlphaTPC,ex,ePurityHelium3AlphaTPC);
    if(AllOrAnti == "All") fPurityHelium3AlphaTPC->SetTitle("Z=2 Particle Purity (TPC ps);<Q_{s}> (a.u.);Purity");
    if(AllOrAnti == "Anti")fPurityHelium3AlphaTPC->SetTitle("Z=2 Anti-Particle Purity (TPC ps);<Q_{s}> (a.u.);Purity");
    fPurityHelium3AlphaTPC->SetMarkerStyle(7);
    fPurityHelium3AlphaTPC->SetMarkerColor(4);
    TGraphErrors *fPurityHelium3AlphaTPCTOF = new TGraphErrors(n,nMeanQRange,PurityHelium3AlphaTPCTOF,ex,ePurityHelium3AlphaTPCTOF);
    if(AllOrAnti == "All") fPurityHelium3AlphaTPCTOF->SetTitle("Z=2 Particle Purity (TPC&TOF ps);<Q_{s}> (a.u.);Purity");
    if(AllOrAnti == "Anti")fPurityHelium3AlphaTPCTOF->SetTitle("Z=2 Anti-Particle Purity (TPC&TOF ps);<Q_{s}> (a.u.);Purity");
    fPurityHelium3AlphaTPCTOF->SetMarkerStyle(7);
    fPurityHelium3AlphaTPCTOF->SetMarkerColor(4);

    TGraphErrors *fPurityDeuTriHe3AlpTPC = new TGraphErrors(n,nMeanQRange,PurityDeuTriHe3AlpTPC,ex,ePurityDeuTriHe3AlpTPC);
    if(AllOrAnti == "All") fPurityDeuTriHe3AlpTPC->SetTitle("d & t & He3 & Alpha Purity (TPC ps);<Q_{s}> (a.u.);Purity");
    if(AllOrAnti == "Anti")fPurityDeuTriHe3AlpTPC->SetTitle("Anti-d & -t & -He3 & -Alpha Purity (TPC ps);<Q_{s}> (a.u.);Purity");
    fPurityDeuTriHe3AlpTPC->SetMarkerStyle(7);
    fPurityDeuTriHe3AlpTPC->SetMarkerColor(4);
    TGraphErrors *fPurityDeuTriHe3AlpTPCTOF = new TGraphErrors(n,nMeanQRange,PurityDeuTriHe3AlpTPCTOF,ex,ePurityDeuTriHe3AlpTPCTOF);
    if(AllOrAnti == "All") fPurityDeuTriHe3AlpTPCTOF->SetTitle("d & t & He3 & Alpha Purity (TPC&TOF ps);<Q_{s}> (a.u.);Purity");
    if(AllOrAnti == "Anti")fPurityDeuTriHe3AlpTPCTOF->SetTitle("Anti-d & -t & -He3 & -Alpha Purity (TPC&TOF ps);<Q_{s}> (a.u.);Purity");
    fPurityDeuTriHe3AlpTPCTOF->SetMarkerStyle(7);
    fPurityDeuTriHe3AlpTPCTOF->SetMarkerColor(4);



  // efficiency
    Double_t mMinBin                          = 260. - m*10;
    Double_t mMeanQRange[m]                   = {0};

    // for the denominator
    Double_t eSigmaEntriesDeuteronTPC         = 0.;
    Double_t eSigmaEntriesDeuteronTPCTOF      = 0.;
    Double_t eSigmaEntriesTritonTPC           = 0.;
    Double_t eSigmaEntriesTritonTPCTOF        = 0.;
    Double_t eSigmaEntriesHelium3TPC          = 0.;
    Double_t eSigmaEntriesHelium3TPCTOF       = 0.;
    Double_t eSigmaEntriesAlphaTPC            = 0.;
    Double_t eSigmaEntriesAlphaTPCTOF         = 0.;
    Double_t eSigmaEntriesHelium3AlphaTPC     = 0.;
    Double_t eSigmaEntriesHelium3AlphaTPCTOF  = 0.;
    Double_t eSigmaEntriesDeuTriHe3AlpTPC     = 0.;
    Double_t eSigmaEntriesDeuTriHe3AlpTPCTOF  = 0.;

    Int_t nSigmaEntriesDeuteronTPC            = fHistTRDpTvPIDvDeuteronTPC_clear       ->IntegralAndError(1,160 ,1,260,eSigmaEntriesDeuteronTPC);        // not correct
    Int_t nSigmaEntriesDeuteronTPCTOF         = fHistTRDpTvPIDvDeuteronTPCTOF_clear    ->IntegralAndError(1,160 ,1,260,eSigmaEntriesDeuteronTPCTOF);     // not correct
    //Int_t nSigmaEntriesDeuteronTPC            = fHistTPCnSigmavRigvDeuteronTPC_clear   ->IntegralAndError(1,2000,1,400,eSigmaEntriesDeuteronTPC);
    //Int_t nSigmaEntriesDeuteronTPCTOF         = fHistTPCnSigmavRigvDeuteronTPCTOF_clear->IntegralAndError(1,2000,1,400,eSigmaEntriesDeuteronTPCTOF);
    Int_t nSigmaEntriesTritonTPC              = fHistTRDpTvPIDvTritonTPC_clear         ->IntegralAndError(1,160 ,1,260,eSigmaEntriesTritonTPC);          // not correct
    Int_t nSigmaEntriesTritonTPCTOF           = fHistTRDpTvPIDvTritonTPCTOF_clear      ->IntegralAndError(1,160 ,1,260,eSigmaEntriesTritonTPCTOF);       // not correct
    //Int_t nSigmaEntriesTritonTPC              = fHistTPCnSigmavRigvTritonTPC_clear     ->IntegralAndError(1,2000,1,400,eSigmaEntriesTritonTPC);
    //Int_t nSigmaEntriesTritonTPCTOF           = fHistTPCnSigmavRigvTritonTPCTOF_clear  ->IntegralAndError(1,2000,1,400,eSigmaEntriesTritonTPCTOF);
    Int_t nSigmaEntriesHelium3TPC             = fHistTRDpTvPIDvHelium3TPC              ->IntegralAndError(1,160 ,1,260,eSigmaEntriesHelium3TPC);         // not correct
    Int_t nSigmaEntriesHelium3TPCTOF          = fHistTRDpTvPIDvHelium3TPCTOF           ->IntegralAndError(1,160 ,1,260,eSigmaEntriesHelium3TPCTOF);      // not correct
    //Int_t nSigmaEntriesHelium3TPC             = fHistTPCnSigmavRigvHelium3TPC          ->IntegralAndError(1,2000,1,400,eSigmaEntriesHelium3TPC);
    //Int_t nSigmaEntriesHelium3TPCTOF          = fHistTPCnSigmavRigvHelium3TPCTOF       ->IntegralAndError(1,2000,1,400,eSigmaEntriesHelium3TPCTOF);
    Int_t nSigmaEntriesAlphaTPC               = fHistTRDpTvPIDvAlphaTPC                ->IntegralAndError(1,160 ,1,260,eSigmaEntriesAlphaTPC);           // not correct
    Int_t nSigmaEntriesAlphaTPCTOF            = fHistTRDpTvPIDvAlphaTPCTOF             ->IntegralAndError(1,160 ,1,260,eSigmaEntriesAlphaTPCTOF);        // not correct
    //Int_t nSigmaEntriesAlphaTPC               = fHistTPCnSigmavRigvAlphaTPC            ->IntegralAndError(1,2000,1,400,eSigmaEntriesAlphaTPC);
    //Int_t nSigmaEntriesAlphaTPCTOF            = fHistTPCnSigmavRigvAlphaTPCTOF         ->IntegralAndError(1,2000,1,400,eSigmaEntriesAlphaTPCTOF);
    Int_t nSigmaEntriesHelium3AlphaTPC        = fHistTRDpTvPIDvHelium3AlphaTPC         ->IntegralAndError(1,160 ,1,260,eSigmaEntriesHelium3AlphaTPC);    // not correct
    Int_t nSigmaEntriesHelium3AlphaTPCTOF     = fHistTRDpTvPIDvHelium3AlphaTPCTOF      ->IntegralAndError(1,160 ,1,260,eSigmaEntriesHelium3AlphaTPCTOF); // not correct
    //Int_t nSigmaEntriesHelium3AlphaTPC        = fHistTPCnSigmavRigvHelium3AlphaTPC     ->IntegralAndError(1,2000,1,400,eSigmaEntriesHelium3AlphaTPC);
    //Int_t nSigmaEntriesHelium3AlphaTPCTOF     = fHistTPCnSigmavRigvHelium3AlphaTPCTOF  ->IntegralAndError(1,2000,1,400,eSigmaEntriesHelium3AlphaTPCTOF);
    Int_t nSigmaEntriesDeuTriHe3AlpTPC        = fHistTRDpTvPIDvDeuTriHe3AlpTPC         ->IntegralAndError(1,160 ,1,260,eSigmaEntriesDeuTriHe3AlpTPC);    // not correct
    Int_t nSigmaEntriesDeuTriHe3AlpTPCTOF     = fHistTRDpTvPIDvDeuTriHe3AlpTPCTOF      ->IntegralAndError(1,160 ,1,260,eSigmaEntriesDeuTriHe3AlpTPCTOF); // not correct
    //Int_t nSigmaEntriesDeuTriHe3AlpTPC        = fHistTPCnSigmavRigvDeuTriHe3AlpTPC     ->IntegralAndError(1,2000,1,400,eSigmaEntriesDeuTriHe3AlpTPC);
    //Int_t nSigmaEntriesDeuTriHe3AlpTPCTOF     = fHistTPCnSigmavRigvDeuTriHe3AlpTPCTOF  ->IntegralAndError(1,2000,1,400,eSigmaEntriesDeuTriHe3AlpTPCTOF);

    // for the numerator
    Double_t ePIDEntriesDeuteronTPC[m]        = {0};
    Double_t ePIDEntriesDeuteronTPCTOF[m]     = {0};
    Double_t ePIDEntriesTritonTPC[m]          = {0};
    Double_t ePIDEntriesTritonTPCTOF[m]       = {0};
    Double_t ePIDEntriesHelium3TPC[m]         = {0};
    Double_t ePIDEntriesHelium3TPCTOF[m]      = {0};
    Double_t ePIDEntriesAlphaTPC[m]           = {0};
    Double_t ePIDEntriesAlphaTPCTOF[m]        = {0};
    Double_t ePIDEntriesHelium3AlphaTPC[m]    = {0};
    Double_t ePIDEntriesHelium3AlphaTPCTOF[m] = {0};
    Double_t ePIDEntriesDeuTriHe3AlpTPC[m]    = {0};
    Double_t ePIDEntriesDeuTriHe3AlpTPCTOF[m] = {0};

    Int_t nPIDEntriesDeuteronTPC[m]           = {0};
    Int_t nPIDEntriesDeuteronTPCTOF[m]        = {0};
    Int_t nPIDEntriesTritonTPC[m]             = {0};
    Int_t nPIDEntriesTritonTPCTOF[m]          = {0};
    Int_t nPIDEntriesHelium3TPC[m]            = {0};
    Int_t nPIDEntriesHelium3TPCTOF[m]         = {0};
    Int_t nPIDEntriesAlphaTPC[m]              = {0};
    Int_t nPIDEntriesAlphaTPCTOF[m]           = {0};
    Int_t nPIDEntriesHelium3AlphaTPC[m]       = {0};
    Int_t nPIDEntriesHelium3AlphaTPCTOF[m]    = {0};
    Int_t nPIDEntriesDeuTriHe3AlpTPC[m]       = {0};
    Int_t nPIDEntriesDeuTriHe3AlpTPCTOF[m]    = {0};

    // for the solution
    Double_t eEfficiencyDeuteronTPC[m]        = {0};
    Double_t eEfficiencyDeuteronTPCTOF[m]     = {0};
    Double_t eEfficiencyTritonTPC[m]          = {0};
    Double_t eEfficiencyTritonTPCTOF[m]       = {0};
    Double_t eEfficiencyHelium3TPC[m]         = {0};
    Double_t eEfficiencyHelium3TPCTOF[m]      = {0};
    Double_t eEfficiencyAlphaTPC[m]           = {0};
    Double_t eEfficiencyAlphaTPCTOF[m]        = {0};
    Double_t eEfficiencyHelium3AlphaTPC[m]    = {0};
    Double_t eEfficiencyHelium3AlphaTPCTOF[m] = {0};
    Double_t eEfficiencyDeuTriHe3AlpTPC[m]    = {0};
    Double_t eEfficiencyDeuTriHe3AlpTPCTOF[m] = {0};

    Double_t EfficiencyDeuteronTPC[m]         = {0};
    Double_t EfficiencyDeuteronTPCTOF[m]      = {0};
    Double_t EfficiencyTritonTPC[m]           = {0};
    Double_t EfficiencyTritonTPCTOF[m]        = {0};
    Double_t EfficiencyHelium3TPC[m]          = {0};
    Double_t EfficiencyHelium3TPCTOF[m]       = {0};
    Double_t EfficiencyAlphaTPC[m]            = {0};
    Double_t EfficiencyAlphaTPCTOF[m]         = {0};
    Double_t EfficiencyHelium3AlphaTPC[m]     = {0};
    Double_t EfficiencyHelium3AlphaTPCTOF[m]  = {0};
    Double_t EfficiencyDeuTriHe3AlpTPC[m]     = {0};
    Double_t EfficiencyDeuTriHe3AlpTPCTOF[m]  = {0};

    Double_t errx[m]                          = {0};


    for(Int_t i = 0; i < m; i++) {
      mMeanQRange[i]                   = mMinBin + i*10;

      nPIDEntriesDeuteronTPC[i]        = fHistTRDpTvPIDvDeuteronTPC_clear   ->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesDeuteronTPC[i]);
      nPIDEntriesDeuteronTPCTOF[i]     = fHistTRDpTvPIDvDeuteronTPCTOF_clear->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesDeuteronTPCTOF[i]);
      nPIDEntriesTritonTPC[i]          = fHistTRDpTvPIDvTritonTPC_clear     ->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesTritonTPC[i]);
      nPIDEntriesTritonTPCTOF[i]       = fHistTRDpTvPIDvTritonTPCTOF_clear  ->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesTritonTPCTOF[i]);
      nPIDEntriesHelium3TPC[i]         = fHistTRDpTvPIDvHelium3TPC          ->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesHelium3TPC[i]);
      nPIDEntriesHelium3TPCTOF[i]      = fHistTRDpTvPIDvHelium3TPCTOF       ->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesHelium3TPCTOF[i]);
      nPIDEntriesAlphaTPC[i]           = fHistTRDpTvPIDvAlphaTPC            ->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesAlphaTPC[i]);
      if(AllOrAnti == "All")nPIDEntriesAlphaTPCTOF[i] = fHistTRDpTvPIDvAlphaTPCTOF->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesAlphaTPCTOF[i]);
      nPIDEntriesHelium3AlphaTPC[i]    = fHistTRDpTvPIDvHelium3AlphaTPC     ->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesHelium3AlphaTPC[i]);
      nPIDEntriesHelium3AlphaTPCTOF[i] = fHistTRDpTvPIDvHelium3AlphaTPCTOF  ->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesHelium3AlphaTPCTOF[i]);
      nPIDEntriesDeuTriHe3AlpTPC[i]    = fHistTRDpTvPIDvDeuTriHe3AlpTPC     ->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesDeuTriHe3AlpTPC[i]);
      nPIDEntriesDeuTriHe3AlpTPCTOF[i] = fHistTRDpTvPIDvDeuTriHe3AlpTPCTOF  ->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesDeuTriHe3AlpTPCTOF[i]);

      EfficiencyDeuteronTPC[i]         = (Double_t)nPIDEntriesDeuteronTPC[i]        / nSigmaEntriesDeuteronTPC;
      EfficiencyDeuteronTPCTOF[i]      = (Double_t)nPIDEntriesDeuteronTPCTOF[i]     / nSigmaEntriesDeuteronTPCTOF;
      EfficiencyTritonTPC[i]           = (Double_t)nPIDEntriesTritonTPC[i]          / nSigmaEntriesTritonTPC;
      EfficiencyTritonTPCTOF[i]        = (Double_t)nPIDEntriesTritonTPCTOF[i]       / nSigmaEntriesTritonTPCTOF;
      EfficiencyHelium3TPC[i]          = (Double_t)nPIDEntriesHelium3TPC[i]         / nSigmaEntriesHelium3TPC;
      EfficiencyHelium3TPCTOF[i]       = (Double_t)nPIDEntriesHelium3TPCTOF[i]      / nSigmaEntriesHelium3TPCTOF;
      EfficiencyAlphaTPC[i]            = (Double_t)nPIDEntriesAlphaTPC[i]           / nSigmaEntriesAlphaTPC;
      if(AllOrAnti == "All")EfficiencyAlphaTPCTOF[i] = (Double_t)nPIDEntriesAlphaTPCTOF[i] / nSigmaEntriesAlphaTPCTOF;
      EfficiencyHelium3AlphaTPC[i]     = (Double_t)nPIDEntriesHelium3AlphaTPC[i]    / nSigmaEntriesHelium3AlphaTPC;
      EfficiencyHelium3AlphaTPCTOF[i]  = (Double_t)nPIDEntriesHelium3AlphaTPCTOF[i] / nSigmaEntriesHelium3AlphaTPCTOF;
      EfficiencyDeuTriHe3AlpTPC[i]     = (Double_t)nPIDEntriesDeuTriHe3AlpTPC[i]    / nSigmaEntriesDeuTriHe3AlpTPC;
      EfficiencyDeuTriHe3AlpTPCTOF[i]  = (Double_t)nPIDEntriesDeuTriHe3AlpTPCTOF[i] / nSigmaEntriesDeuTriHe3AlpTPCTOF;
                                                                                                                                                       //own function, see below
      eEfficiencyDeuteronTPC[i]        = GaussError(nPIDEntriesDeuteronTPC[i]       ,ePIDEntriesDeuteronTPC[i]       ,nSigmaEntriesDeuteronTPC       ,eSigmaEntriesDeuteronTPC);
      eEfficiencyDeuteronTPCTOF[i]     = GaussError(nPIDEntriesDeuteronTPCTOF[i]    ,ePIDEntriesDeuteronTPCTOF[i]    ,nSigmaEntriesDeuteronTPCTOF    ,eSigmaEntriesDeuteronTPCTOF);
      eEfficiencyTritonTPC[i]          = GaussError(nPIDEntriesTritonTPC[i]         ,ePIDEntriesTritonTPC[i]         ,nSigmaEntriesTritonTPC         ,eSigmaEntriesTritonTPC);
      eEfficiencyTritonTPCTOF[i]       = GaussError(nPIDEntriesTritonTPCTOF[i]      ,ePIDEntriesTritonTPCTOF[i]      ,nSigmaEntriesTritonTPCTOF      ,eSigmaEntriesTritonTPCTOF);
      eEfficiencyHelium3TPC[i]         = GaussError(nPIDEntriesHelium3TPC[i]        ,ePIDEntriesHelium3TPC[i]        ,nSigmaEntriesHelium3TPC        ,eSigmaEntriesHelium3TPC);
      eEfficiencyHelium3TPCTOF[i]      = GaussError(nPIDEntriesHelium3TPCTOF[i]     ,ePIDEntriesHelium3TPCTOF[i]     ,nSigmaEntriesHelium3TPCTOF     ,eSigmaEntriesHelium3TPCTOF);
      eEfficiencyAlphaTPC[i]           = GaussError(nPIDEntriesAlphaTPC[i]          ,ePIDEntriesAlphaTPC[i]          ,nSigmaEntriesAlphaTPC          ,eSigmaEntriesAlphaTPC);
      if(AllOrAnti == "All")eEfficiencyAlphaTPCTOF[i] = GaussError(nPIDEntriesAlphaTPCTOF[i],ePIDEntriesAlphaTPCTOF[i],nSigmaEntriesAlphaTPCTOF      ,eSigmaEntriesAlphaTPCTOF);
      eEfficiencyHelium3AlphaTPC[i]    = GaussError(nPIDEntriesHelium3AlphaTPC[i]   ,ePIDEntriesHelium3AlphaTPC[i]   ,nSigmaEntriesHelium3AlphaTPC   ,eSigmaEntriesHelium3AlphaTPC);
      eEfficiencyHelium3AlphaTPCTOF[i] = GaussError(nPIDEntriesHelium3AlphaTPCTOF[i],ePIDEntriesHelium3AlphaTPCTOF[i],nSigmaEntriesHelium3AlphaTPCTOF,eSigmaEntriesHelium3AlphaTPCTOF);
      eEfficiencyDeuTriHe3AlpTPC[i]    = GaussError(nPIDEntriesDeuTriHe3AlpTPC[i]   ,ePIDEntriesDeuTriHe3AlpTPC[i]   ,nSigmaEntriesDeuTriHe3AlpTPC   ,eSigmaEntriesDeuTriHe3AlpTPC);
      eEfficiencyDeuTriHe3AlpTPCTOF[i] = GaussError(nPIDEntriesDeuTriHe3AlpTPCTOF[i],ePIDEntriesDeuTriHe3AlpTPCTOF[i],nSigmaEntriesDeuTriHe3AlpTPCTOF,eSigmaEntriesDeuTriHe3AlpTPCTOF);

      errx[i]                          = 0.;

      //cout <<"Deuteron efficiency with TPC preselection and PID cut at "               << mMeanQRange[i] <<":\t("<< EfficiencyDeuteronTPC[i]       <<"\t+- "<<eEfficiencyDeuteronTPC[i]       <<")"<< endl;
      //cout <<"Deuteron efficiency with TPC&TOF preselection and PID cut at "           << mMeanQRange[i] <<":\t("<< EfficiencyDeuteronTPCTOF[i]    <<"\t+- "<<eEfficiencyDeuteronTPCTOF[i]    <<")"<< endl;
      //cout <<"Triton efficiency with TPC preselection and PID cut at "                 << mMeanQRange[i] <<":\t("<< EfficiencyTritonTPC[i]         <<"\t+- "<<eEfficiencyTritonTPC[i]         <<")"<< endl;
      //cout <<"Triton efficiency with TPC&TOF preselection and PID cut at "             << mMeanQRange[i] <<":\t("<< EfficiencyTritonTPCTOF[i]      <<"\t+- "<<eEfficiencyTritonTPCTOF[i]      <<")"<< endl;
      //cout <<"Helium3 efficiency with TPC preselection and PID cut at "                << mMeanQRange[i] <<":\t("<< EfficiencyHelium3TPC[i]        <<"\t+- "<<eEfficiencyHelium3TPC[i]        <<")"<< endl;
      //cout <<"Helium3 efficiency with TPC&TOF preselection and PID cut at "            << mMeanQRange[i] <<":\t("<< EfficiencyHelium3TPCTOF[i]     <<"\t+- "<<eEfficiencyHelium3TPCTOF[i]     <<")"<< endl;
      //cout <<"Alpha efficiency with TPC preselection and PID cut at "                  << mMeanQRange[i] <<":\t("<< EfficiencyAlphaTPC[i]          <<"\t+- "<<eEfficiencyAlphaTPC[i]          <<")"<< endl;
      //cout <<"Alpha efficiency with TPC&TOF preselection and PID cut at "              << mMeanQRange[i] <<":\t("<< EfficiencyAlphaTPCTOF[i]       <<"\t+- "<<eEfficiencyAlphaTPCTOF[i]       <<")"<< endl;
      //cout <<"Z=2 particle efficiency with TPC preselection and PID cut at "           << mMeanQRange[i] <<":\t("<< EfficiencyHelium3AlphaTPC[i]   <<"\t+- "<<eEfficiencyHelium3AlphaTPC[i]   <<")"<< endl;
      //cout <<"Z=2 particle efficiency with TPC&TOF preselection and PID cut at "       << mMeanQRange[i] <<":\t("<< EfficiencyHelium3AlphaTPCTOF[i]<<"\t+- "<<eEfficiencyHelium3AlphaTPCTOF[i]<<")"<< endl;
      //cout <<"d & t & He3 & Alpha efficiency with TPC preselection and PID cut at "    << mMeanQRange[i] <<":\t("<< EfficiencyDeuteronTPC[i]       <<"\t+- "<<eEfficiencyDeuTriHe3AlpTPC[i]   <<")"<< endl;
      //cout <<"d & t & He3 & Alpha efficiency with TPC&TOF preselection and PID cut at "<< mMeanQRange[i] <<":\t("<< EfficiencyDeuteronTPCTOF[i]    <<"\t+- "<<eEfficiencyDeuTriHe3AlpTPCTOF[i]<<")"<< endl;
    }


    TGraphErrors *fEfficiencyDeuteronTPC = new TGraphErrors(m,mMeanQRange,EfficiencyDeuteronTPC,errx,eEfficiencyDeuteronTPC);
    if(AllOrAnti == "All") fEfficiencyDeuteronTPC->SetTitle("Deuteron Efficiency (TPC ps);<Q_{s}> (a.u.);Efficiency");
    if(AllOrAnti == "Anti")fEfficiencyDeuteronTPC->SetTitle("Anti-Deuteron Efficiency (TPC ps);<Q_{s}> (a.u.);Efficiency");
    fEfficiencyDeuteronTPC->SetMarkerStyle(7);
    fEfficiencyDeuteronTPC->SetMarkerColor(4);
    TGraphErrors *fEfficiencyDeuteronTPCTOF = new TGraphErrors(m,mMeanQRange,EfficiencyDeuteronTPCTOF,errx,eEfficiencyDeuteronTPCTOF);
    if(AllOrAnti == "All") fEfficiencyDeuteronTPCTOF->SetTitle("Deuteron Efficiency (TPC&TOF ps);<Q_{s}> (a.u.);Efficiency");
    if(AllOrAnti == "Anti")fEfficiencyDeuteronTPCTOF->SetTitle("Anti-Deuteron Efficiency (TPC&TOF ps);<Q_{s}> (a.u.);Efficiency");
    fEfficiencyDeuteronTPCTOF->SetMarkerStyle(7);
    fEfficiencyDeuteronTPCTOF->SetMarkerColor(4);

    TGraphErrors *fEfficiencyTritonTPC = new TGraphErrors(m,mMeanQRange,EfficiencyTritonTPC,errx,eEfficiencyTritonTPC);
    if(AllOrAnti == "All") fEfficiencyTritonTPC->SetTitle("Triton Efficiency (TPC ps);<Q_{s}> (a.u.);Efficiency");
    if(AllOrAnti == "Anti")fEfficiencyTritonTPC->SetTitle("Anti-Triton Efficiency (TPC ps);<Q_{s}> (a.u.);Efficiency");
    fEfficiencyTritonTPC->SetMarkerStyle(7);
    fEfficiencyTritonTPC->SetMarkerColor(4);
    TGraphErrors *fEfficiencyTritonTPCTOF = new TGraphErrors(m,mMeanQRange,EfficiencyTritonTPCTOF,errx,eEfficiencyTritonTPCTOF);
    if(AllOrAnti == "All") fEfficiencyTritonTPCTOF->SetTitle("Triton Efficiency (TPC&TOF ps);<Q_{s}> (a.u.);Efficiency");
    if(AllOrAnti == "Anti")fEfficiencyTritonTPCTOF->SetTitle("Anti-Triton Efficiency (TPC&TOF ps);<Q_{s}> (a.u.);Efficiency");
    fEfficiencyTritonTPCTOF->SetMarkerStyle(7);
    fEfficiencyTritonTPCTOF->SetMarkerColor(4);

    TGraphErrors *fEfficiencyHelium3TPC = new TGraphErrors(m,mMeanQRange,EfficiencyHelium3TPC,errx,eEfficiencyHelium3TPC);
    if(AllOrAnti == "All") fEfficiencyHelium3TPC->SetTitle("Helium3 Efficiency (TPC ps);<Q_{s}> (a.u.);Efficiency");
    if(AllOrAnti == "Anti")fEfficiencyHelium3TPC->SetTitle("Anti-Helium3 Efficiency (TPC ps);<Q_{s}> (a.u.);Efficiency");
    fEfficiencyHelium3TPC->SetMarkerStyle(7);
    fEfficiencyHelium3TPC->SetMarkerColor(4);
    TGraphErrors *fEfficiencyHelium3TPCTOF = new TGraphErrors(m,mMeanQRange,EfficiencyHelium3TPCTOF,errx,eEfficiencyHelium3TPCTOF);
    if(AllOrAnti == "All") fEfficiencyHelium3TPCTOF->SetTitle("Helium3 Efficiency (TPC&TOF ps);<Q_{s}> (a.u.);Efficiency");
    if(AllOrAnti == "Anti")fEfficiencyHelium3TPCTOF->SetTitle("Anti-Helium3 Efficiency (TPC&TOF ps);<Q_{s}> (a.u.);Efficiency");
    fEfficiencyHelium3TPCTOF->SetMarkerStyle(7);
    fEfficiencyHelium3TPCTOF->SetMarkerColor(4);

    TGraphErrors *fEfficiencyAlphaTPC = new TGraphErrors(m,mMeanQRange,EfficiencyAlphaTPC,errx,eEfficiencyAlphaTPC);
    if(AllOrAnti == "All") fEfficiencyAlphaTPC->SetTitle("Alpha Efficiency (TPC ps);<Q_{s}> (a.u.);Efficiency");
    if(AllOrAnti == "Anti")fEfficiencyAlphaTPC->SetTitle("Anti-Alpha Efficiency (TPC ps);<Q_{s}> (a.u.);Efficiency");
    fEfficiencyAlphaTPC->SetMarkerStyle(7);
    fEfficiencyAlphaTPC->SetMarkerColor(4);
    TGraphErrors *fEfficiencyAlphaTPCTOF = new TGraphErrors(m,mMeanQRange,EfficiencyAlphaTPCTOF,errx,eEfficiencyAlphaTPCTOF);
    if(AllOrAnti == "All") fEfficiencyAlphaTPCTOF->SetTitle("Alpha Efficiency (TPC&TOF ps);<Q_{s}> (a.u.);Efficiency");
    if(AllOrAnti == "Anti")fEfficiencyAlphaTPCTOF->SetTitle("Anti-Alpha Efficiency (TPC&TOF ps);<Q_{s}> (a.u.);Efficiency");
    fEfficiencyAlphaTPCTOF->SetMarkerStyle(7);
    fEfficiencyAlphaTPCTOF->SetMarkerColor(4);

    TGraphErrors *fEfficiencyHelium3AlphaTPC = new TGraphErrors(m,mMeanQRange,EfficiencyHelium3AlphaTPC,errx,eEfficiencyHelium3AlphaTPC);
    if(AllOrAnti == "All") fEfficiencyHelium3AlphaTPC->SetTitle("Z=2 Particle Efficiency (TPC ps);<Q_{s}> (a.u.);Efficiency");
    if(AllOrAnti == "Anti")fEfficiencyHelium3AlphaTPC->SetTitle("Z=2 Anti-Particle Efficiency (TPC ps);<Q_{s}> (a.u.);Efficiency");
    fEfficiencyHelium3AlphaTPC->SetMarkerStyle(7);
    fEfficiencyHelium3AlphaTPC->SetMarkerColor(4);
    TGraphErrors *fEfficiencyHelium3AlphaTPCTOF = new TGraphErrors(m,mMeanQRange,EfficiencyHelium3AlphaTPCTOF,errx,eEfficiencyHelium3AlphaTPCTOF);
    if(AllOrAnti == "All") fEfficiencyHelium3AlphaTPCTOF->SetTitle("Z=2 Particle Efficiency (TPC&TOF ps);<Q_{s}> (a.u.);Efficiency");
    if(AllOrAnti == "Anti")fEfficiencyHelium3AlphaTPCTOF->SetTitle("Z=2 Anti-Particle Efficiency (TPC&TOF ps);<Q_{s}> (a.u.);Efficiency");
    fEfficiencyHelium3AlphaTPCTOF->SetMarkerStyle(7);
    fEfficiencyHelium3AlphaTPCTOF->SetMarkerColor(4);

    TGraphErrors *fEfficiencyDeuTriHe3AlpTPC = new TGraphErrors(m,mMeanQRange,EfficiencyDeuTriHe3AlpTPC,errx,eEfficiencyDeuTriHe3AlpTPC);
    if(AllOrAnti == "All") fEfficiencyDeuTriHe3AlpTPC->SetTitle("d & t & He3 & Alpha Efficiency (TPC ps);<Q_{s}> (a.u.);Efficiency");
    if(AllOrAnti == "Anti")fEfficiencyDeuTriHe3AlpTPC->SetTitle("Anti-d & -t & -He3 & -Alpha Efficiency (TPC ps);<Q_{s}> (a.u.);Efficiency");
    fEfficiencyDeuTriHe3AlpTPC->SetMarkerStyle(7);
    fEfficiencyDeuTriHe3AlpTPC->SetMarkerColor(4);
    TGraphErrors *fEfficiencyDeuTriHe3AlpTPCTOF = new TGraphErrors(m,mMeanQRange,EfficiencyDeuTriHe3AlpTPCTOF,errx,eEfficiencyDeuTriHe3AlpTPCTOF);
    if(AllOrAnti == "All") fEfficiencyDeuTriHe3AlpTPCTOF->SetTitle("d & t & He3 & Alpha Efficiency (TPC&TOF ps);<Q_{s}> (a.u.);Efficiency");
    if(AllOrAnti == "Anti")fEfficiencyDeuTriHe3AlpTPCTOF->SetTitle("Anti-d & -t & -He3 & -Alpha Efficiency (TPC&TOF ps);<Q_{s}> (a.u.);Efficiency");
    fEfficiencyDeuTriHe3AlpTPCTOF->SetMarkerStyle(7);
    fEfficiencyDeuTriHe3AlpTPCTOF->SetMarkerColor(4);


  // rejection
    // for the denominator
    Double_t eSigmaEntriesAllexceptDeuteronTPC         = 0.;
    Double_t eSigmaEntriesAllexceptDeuteronTPCTOF      = 0.;
    Double_t eSigmaEntriesAllexceptTritonTPC           = 0.;
    Double_t eSigmaEntriesAllexceptTritonTPCTOF        = 0.;
    Double_t eSigmaEntriesAllexceptHelium3TPC          = 0.;
    Double_t eSigmaEntriesAllexceptHelium3TPCTOF       = 0.;
    Double_t eSigmaEntriesAllexceptAlphaTPC            = 0.;
    Double_t eSigmaEntriesAllexceptAlphaTPCTOF         = 0.;
    Double_t eSigmaEntriesAllexceptHelium3AlphaTPC     = 0.;
    Double_t eSigmaEntriesAllexceptHelium3AlphaTPCTOF  = 0.;
    Double_t eSigmaEntriesAllexceptDeuTriHe3AlpTPC     = 0.;
    Double_t eSigmaEntriesAllexceptDeuTriHe3AlpTPCTOF  = 0.;

    //Int_t nSigmaEntriesAllexceptDeuteronTPC            = fHistTPCnSigmavRigvAllexceptDeuteronTPC_clear   ->IntegralAndError(1,2000,1,400,eSigmaEntriesAllexceptDeuteronTPC);
    //Int_t nSigmaEntriesAllexceptDeuteronTPCTOF         = fHistTPCnSigmavRigvAllexceptDeuteronTPCTOF_clear->IntegralAndError(1,2000,1,400,eSigmaEntriesAllexceptDeuteronTPCTOF);
    Int_t nSigmaEntriesAllexceptDeuteronTPC            = fHistTRDpTvPIDvAllexceptDeuteronTPC_clear   ->IntegralAndError(1,500,1,260,eSigmaEntriesAllexceptDeuteronTPC);    // not correct
    Int_t nSigmaEntriesAllexceptDeuteronTPCTOF         = fHistTRDpTvPIDvAllexceptDeuteronTPCTOF_clear->IntegralAndError(1,500,1,260,eSigmaEntriesAllexceptDeuteronTPCTOF); // not correct
    //Int_t nSigmaEntriesAllexceptTritonTPC              = fHistTPCnSigmavRigvAllexceptTritonTPC_clear     ->IntegralAndError(1,2000,1,400,eSigmaEntriesAllexceptTritonTPC);
    //Int_t nSigmaEntriesAllexceptTritonTPCTOF           = fHistTPCnSigmavRigvAllexceptTritonTPCTOF_clear  ->IntegralAndError(1,2000,1,400,eSigmaEntriesAllexceptTritonTPCTOF);
    Int_t nSigmaEntriesAllexceptTritonTPC              = fHistTRDpTvPIDvAllexceptTritonTPC_clear     ->IntegralAndError(1,500,1,260,eSigmaEntriesAllexceptTritonTPC);      // not correct
    Int_t nSigmaEntriesAllexceptTritonTPCTOF           = fHistTRDpTvPIDvAllexceptTritonTPCTOF_clear  ->IntegralAndError(1,500,1,260,eSigmaEntriesAllexceptTritonTPCTOF);   // not correct
    //Int_t nSigmaEntriesAllexceptHelium3TPC             = fHistTPCnSigmavRigvAllexceptHelium3TPC          ->IntegralAndError(1,2000,1,400,eSigmaEntriesAllexceptHelium3TPC);
    //Int_t nSigmaEntriesAllexceptHelium3TPCTOF          = fHistTPCnSigmavRigvAllexceptHelium3TPCTOF       ->IntegralAndError(1,2000,1,400,eSigmaEntriesAllexceptHelium3TPCTOF);
    Int_t nSigmaEntriesAllexceptHelium3TPC             = fHistTRDpTvPIDvAllexceptHelium3TPC          ->IntegralAndError(1,500,1,260,eSigmaEntriesHelium3TPC);              // not correct
    Int_t nSigmaEntriesAllexceptHelium3TPCTOF          = fHistTRDpTvPIDvAllexceptHelium3TPCTOF       ->IntegralAndError(1,500,1,260,eSigmaEntriesHelium3TPCTOF);           // not correct
    //Int_t nSigmaEntriesAllexceptAlphaTPC               = fHistTPCnSigmavRigvAllexceptAlphaTPC            ->IntegralAndError(1,2000,1,400,eSigmaEntriesAlphaTPC);
    //Int_t nSigmaEntriesAllexceptAlphaTPCTOF            = fHistTPCnSigmavRigvAllexceptAlphaTPCTOF         ->IntegralAndError(1,2000,1,400,eSigmaEntriesAlphaTPCTOF);
    Int_t nSigmaEntriesAllexceptAlphaTPC               = fHistTRDpTvPIDvAllexceptAlphaTPC            ->IntegralAndError(1,500,1,260,eSigmaEntriesAlphaTPC);                // not correct
    Int_t nSigmaEntriesAllexceptAlphaTPCTOF            = fHistTRDpTvPIDvAllexceptAlphaTPCTOF         ->IntegralAndError(1,500,1,260,eSigmaEntriesAlphaTPCTOF);             // not correct
    //Int_t nSigmaEntriesAllexceptHelium3AlphaTPC        = fHistTPCnSigmavRigvAllexceptHelium3AlphaTPC     ->IntegralAndError(1,2000,1,400,eSigmaEntriesHelium3AlphaTPC);
    //Int_t nSigmaEntriesAllexceptHelium3AlphaTPCTOF     = fHistTPCnSigmavRigvAllexceptHelium3AlphaTPCTOF  ->IntegralAndError(1,2000,1,400,eSigmaEntriesHelium3AlphaTPCTOF);
    Int_t nSigmaEntriesAllexceptHelium3AlphaTPC        = fHistTRDpTvPIDvAllexceptHelium3AlphaTPC         ->IntegralAndError(1,500 ,1,260,eSigmaEntriesHelium3AlphaTPC);    // not correct
    Int_t nSigmaEntriesAllexceptHelium3AlphaTPCTOF     = fHistTRDpTvPIDvAllexceptHelium3AlphaTPCTOF      ->IntegralAndError(1,500 ,1,260,eSigmaEntriesHelium3AlphaTPCTOF); // not correct
    //Int_t nSigmaEntriesAllexceptDeuTriHe3AlpTPC        = fHistTPCnSigmavRigvAllexceptDeuTriHe3AlpTPC     ->IntegralAndError(1,2000,1,400,eSigmaEntriesDeuTriHe3AlpTPC);
    //Int_t nSigmaEntriesAllexceptDeuTriHe3AlpTPCTOF     = fHistTPCnSigmavRigvAllexceptDeuTriHe3AlpTPCTOF  ->IntegralAndError(1,2000,1,400,eSigmaEntriesDeuTriHe3AlpTPCTOF);
    Int_t nSigmaEntriesAllexceptDeuTriHe3AlpTPC        = fHistTRDpTvPIDvAllexceptDeuTriHe3AlpTPC         ->IntegralAndError(1,500 ,1,260,eSigmaEntriesDeuTriHe3AlpTPC);    // not correct
    Int_t nSigmaEntriesAllexceptDeuTriHe3AlpTPCTOF     = fHistTRDpTvPIDvAllexceptDeuTriHe3AlpTPCTOF      ->IntegralAndError(1,500 ,1,260,eSigmaEntriesDeuTriHe3AlpTPCTOF); // not correct

    // for the numerator
    Double_t ePIDEntriesAllexceptDeuteronTPC[m]        = {0};
    Double_t ePIDEntriesAllexceptDeuteronTPCTOF[m]     = {0};
    Double_t ePIDEntriesAllexceptTritonTPC[m]          = {0};
    Double_t ePIDEntriesAllexceptTritonTPCTOF[m]       = {0};
    Double_t ePIDEntriesAllexceptHelium3TPC[m]         = {0};
    Double_t ePIDEntriesAllexceptHelium3TPCTOF[m]      = {0};
    Double_t ePIDEntriesAllexceptAlphaTPC[m]           = {0};
    Double_t ePIDEntriesAllexceptAlphaTPCTOF[m]        = {0};
    Double_t ePIDEntriesAllexceptHelium3AlphaTPC[m]    = {0};
    Double_t ePIDEntriesAllexceptHelium3AlphaTPCTOF[m] = {0};
    Double_t ePIDEntriesAllexceptDeuTriHe3AlpTPC[m]    = {0};
    Double_t ePIDEntriesAllexceptDeuTriHe3AlpTPCTOF[m] = {0};

    Int_t nPIDEntriesAllexceptDeuteronTPC[m]           = {0};
    Int_t nPIDEntriesAllexceptDeuteronTPCTOF[m]        = {0};
    Int_t nPIDEntriesAllexceptTritonTPC[m]             = {0};
    Int_t nPIDEntriesAllexceptTritonTPCTOF[m]          = {0};
    Int_t nPIDEntriesAllexceptHelium3TPC[m]            = {0};
    Int_t nPIDEntriesAllexceptHelium3TPCTOF[m]         = {0};
    Int_t nPIDEntriesAllexceptAlphaTPC[m]              = {0};
    Int_t nPIDEntriesAllexceptAlphaTPCTOF[m]           = {0};
    Int_t nPIDEntriesAllexceptHelium3AlphaTPC[m]       = {0};
    Int_t nPIDEntriesAllexceptHelium3AlphaTPCTOF[m]    = {0};
    Int_t nPIDEntriesAllexceptDeuTriHe3AlpTPC[m]       = {0};
    Int_t nPIDEntriesAllexceptDeuTriHe3AlpTPCTOF[m]    = {0};

    // for the solution
    Double_t eRejectionDeuteronTPC[m]                  = {0};
    Double_t eRejectionDeuteronTPCTOF[m]               = {0};
    Double_t eRejectionTritonTPC[m]                    = {0};
    Double_t eRejectionTritonTPCTOF[m]                 = {0};
    Double_t eRejectionHelium3TPC[m]                   = {0};
    Double_t eRejectionHelium3TPCTOF[m]                = {0};
    Double_t eRejectionAlphaTPC[m]                     = {0};
    Double_t eRejectionAlphaTPCTOF[m]                  = {0};
    Double_t eRejectionHelium3AlphaTPC[m]              = {0};
    Double_t eRejectionHelium3AlphaTPCTOF[m]           = {0};
    Double_t eRejectionDeuTriHe3AlpTPC[m]              = {0};
    Double_t eRejectionDeuTriHe3AlpTPCTOF[m]           = {0};

    Double_t RejectionDeuteronTPC[m]                   = {0};
    Double_t RejectionDeuteronTPCTOF[m]                = {0};
    Double_t RejectionTritonTPC[m]                     = {0};
    Double_t RejectionTritonTPCTOF[m]                  = {0};
    Double_t RejectionHelium3TPC[m]                    = {0};
    Double_t RejectionHelium3TPCTOF[m]                 = {0};
    Double_t RejectionAlphaTPC[m]                      = {0};
    Double_t RejectionAlphaTPCTOF[m]                   = {0};
    Double_t RejectionHelium3AlphaTPC[m]               = {0};
    Double_t RejectionHelium3AlphaTPCTOF[m]            = {0};
    Double_t RejectionDeuTriHe3AlpTPC[m]               = {0};
    Double_t RejectionDeuTriHe3AlpTPCTOF[m]            = {0};


    for(Int_t i = 0; i < m; i++) {
      mMeanQRange[i]                            = mMinBin + i*10;

      nPIDEntriesAllexceptDeuteronTPC[i]        = fHistTRDpTvPIDvAllexceptDeuteronTPC_clear   ->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesAllexceptDeuteronTPC[i]);
      nPIDEntriesAllexceptDeuteronTPCTOF[i]     = fHistTRDpTvPIDvAllexceptDeuteronTPCTOF_clear->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesAllexceptDeuteronTPCTOF[i]);
      nPIDEntriesAllexceptTritonTPC[i]          = fHistTRDpTvPIDvAllexceptTritonTPC_clear     ->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesAllexceptTritonTPC[i]);
      nPIDEntriesAllexceptTritonTPCTOF[i]       = fHistTRDpTvPIDvAllexceptTritonTPCTOF_clear  ->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesAllexceptTritonTPCTOF[i]);
      nPIDEntriesAllexceptHelium3TPC[i]         = fHistTRDpTvPIDvAllexceptHelium3TPC          ->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesAllexceptHelium3TPC[i]);
      nPIDEntriesAllexceptHelium3TPCTOF[i]      = fHistTRDpTvPIDvAllexceptHelium3TPCTOF       ->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesAllexceptHelium3TPCTOF[i]);
      nPIDEntriesAllexceptAlphaTPC[i]           = fHistTRDpTvPIDvAllexceptAlphaTPC            ->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesAllexceptAlphaTPC[i]);
      nPIDEntriesAllexceptAlphaTPCTOF[i]        = fHistTRDpTvPIDvAllexceptAlphaTPCTOF         ->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesAllexceptAlphaTPCTOF[i]);
      nPIDEntriesAllexceptHelium3AlphaTPC[i]    = fHistTRDpTvPIDvAllexceptHelium3AlphaTPC     ->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesAllexceptHelium3AlphaTPC[i]);
      nPIDEntriesAllexceptHelium3AlphaTPCTOF[i] = fHistTRDpTvPIDvAllexceptHelium3AlphaTPCTOF  ->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesAllexceptHelium3AlphaTPCTOF[i]);
      nPIDEntriesAllexceptDeuTriHe3AlpTPC[i]    = fHistTRDpTvPIDvAllexceptDeuTriHe3AlpTPC     ->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesAllexceptDeuTriHe3AlpTPC[i]);
      nPIDEntriesAllexceptDeuTriHe3AlpTPCTOF[i] = fHistTRDpTvPIDvAllexceptDeuTriHe3AlpTPCTOF  ->IntegralAndError(1,160,mMeanQRange[i],260,ePIDEntriesAllexceptDeuTriHe3AlpTPCTOF[i]);

      RejectionDeuteronTPC[i]         = (Double_t)nPIDEntriesAllexceptDeuteronTPC[i]        / nSigmaEntriesAllexceptDeuteronTPC;
      RejectionDeuteronTPCTOF[i]      = (Double_t)nPIDEntriesAllexceptDeuteronTPCTOF[i]     / nSigmaEntriesAllexceptDeuteronTPCTOF;
      RejectionTritonTPC[i]           = (Double_t)nPIDEntriesAllexceptTritonTPC[i]          / nSigmaEntriesAllexceptTritonTPC;
      RejectionTritonTPCTOF[i]        = (Double_t)nPIDEntriesAllexceptTritonTPCTOF[i]       / nSigmaEntriesAllexceptTritonTPCTOF;
      RejectionHelium3TPC[i]          = (Double_t)nPIDEntriesAllexceptHelium3TPC[i]         / nSigmaEntriesAllexceptHelium3TPC;
      RejectionHelium3TPCTOF[i]       = (Double_t)nPIDEntriesAllexceptHelium3TPCTOF[i]      / nSigmaEntriesAllexceptHelium3TPCTOF;
      RejectionAlphaTPC[i]            = (Double_t)nPIDEntriesAllexceptAlphaTPC[i]           / nSigmaEntriesAllexceptAlphaTPC;
      RejectionAlphaTPCTOF[i]         = (Double_t)nPIDEntriesAllexceptAlphaTPCTOF[i]        / nSigmaEntriesAllexceptAlphaTPCTOF;
      RejectionHelium3AlphaTPC[i]     = (Double_t)nPIDEntriesAllexceptHelium3AlphaTPC[i]    / nSigmaEntriesAllexceptHelium3AlphaTPC;
      RejectionHelium3AlphaTPCTOF[i]  = (Double_t)nPIDEntriesAllexceptHelium3AlphaTPCTOF[i] / nSigmaEntriesAllexceptHelium3AlphaTPCTOF;
      RejectionDeuTriHe3AlpTPC[i]     = (Double_t)nPIDEntriesAllexceptDeuTriHe3AlpTPC[i]    / nSigmaEntriesAllexceptDeuTriHe3AlpTPC;
      RejectionDeuTriHe3AlpTPCTOF[i]  = (Double_t)nPIDEntriesAllexceptDeuTriHe3AlpTPCTOF[i] / nSigmaEntriesAllexceptDeuTriHe3AlpTPCTOF;
                                                                                                                                                                                      //own function, see below
      eRejectionDeuteronTPC[i]        = GaussError(nPIDEntriesAllexceptDeuteronTPC[i]       ,ePIDEntriesAllexceptDeuteronTPC[i]       ,nSigmaEntriesAllexceptDeuteronTPC    ,eSigmaEntriesAllexceptDeuteronTPC);
      eRejectionDeuteronTPCTOF[i]     = GaussError(nPIDEntriesAllexceptDeuteronTPCTOF[i]    ,ePIDEntriesAllexceptDeuteronTPCTOF[i]   ,nSigmaEntriesAllexceptDeuteronTPCTOF,eSigmaEntriesAllexceptDeuteronTPCTOF);
      eRejectionTritonTPC[i]          = GaussError(nPIDEntriesAllexceptTritonTPC[i]         ,ePIDEntriesAllexceptTritonTPC[i]         ,nSigmaEntriesAllexceptTritonTPC         ,eSigmaEntriesAllexceptTritonTPC);
      eRejectionTritonTPCTOF[i]       = GaussError(nPIDEntriesAllexceptTritonTPCTOF[i]      ,ePIDEntriesAllexceptTritonTPCTOF[i]      ,nSigmaEntriesAllexceptTritonTPCTOF   ,eSigmaEntriesAllexceptTritonTPCTOF);
      eRejectionHelium3TPC[i]         = GaussError(nPIDEntriesAllexceptHelium3TPC[i]        ,ePIDEntriesAllexceptHelium3TPC[i]        ,nSigmaEntriesAllexceptHelium3TPC       ,eSigmaEntriesAllexceptHelium3TPC);
      eRejectionHelium3TPCTOF[i]      = GaussError(nPIDEntriesAllexceptHelium3TPCTOF[i]     ,ePIDEntriesAllexceptHelium3TPCTOF[i]     ,nSigmaEntriesAllexceptHelium3TPCTOF ,eSigmaEntriesAllexceptHelium3TPCTOF);
      eRejectionAlphaTPC[i]           = GaussError(nPIDEntriesAllexceptAlphaTPC[i]          ,ePIDEntriesAllexceptAlphaTPC[i]          ,nSigmaEntriesAllexceptAlphaTPC          ,eSigmaEntriesAllexceptAlphaTPC);
      eRejectionAlphaTPCTOF[i]        = GaussError(nPIDEntriesAllexceptAlphaTPCTOF[i]       ,ePIDEntriesAllexceptAlphaTPCTOF[i]       ,nSigmaEntriesAllexceptAlphaTPCTOF     ,eSigmaEntriesAllexceptAlphaTPCTOF);
      eRejectionHelium3AlphaTPC[i]    = GaussError(nPIDEntriesAllexceptHelium3AlphaTPC[i]   ,ePIDEntriesAllexceptHelium3AlphaTPC[i],nSigmaEntriesAllexceptHelium3AlphaTPC,eSigmaEntriesAllexceptHelium3AlphaTPC);
      eRejectionHelium3AlphaTPCTOF[i] = GaussError(nPIDEntriesAllexceptHelium3AlphaTPCTOF[i],ePIDEntriesAllexceptHelium3AlphaTPCTOF[i],nSigmaEntriesAllexceptHelium3AlphaTPCTOF,eSigmaEntriesAllexceptHelium3AlphaTPCTOF);
      eRejectionDeuTriHe3AlpTPC[i]    = GaussError(nPIDEntriesAllexceptDeuTriHe3AlpTPC[i]   ,ePIDEntriesAllexceptDeuTriHe3AlpTPC[i],nSigmaEntriesAllexceptDeuTriHe3AlpTPC,eSigmaEntriesAllexceptDeuTriHe3AlpTPC);
      eRejectionDeuTriHe3AlpTPCTOF[i] = GaussError(nPIDEntriesAllexceptDeuTriHe3AlpTPCTOF[i],ePIDEntriesAllexceptDeuTriHe3AlpTPCTOF[i],nSigmaEntriesAllexceptDeuTriHe3AlpTPCTOF,eSigmaEntriesAllexceptDeuTriHe3AlpTPCTOF);

      errx[i]                          = 0.;

      //cout <<"Deuteron rejection with TPC preselection and PID cut at "               << mMeanQRange[i] <<":\t("<< RejectionDeuteronTPC[i]       <<"\t+- "<<eRejectionDeuteronTPC[i]       <<")"<< endl;
      //cout <<"Deuteron rejection with TPC&TOF preselection and PID cut at "           << mMeanQRange[i] <<":\t("<< RejectionDeuteronTPCTOF[i]    <<"\t+- "<<eRejectionDeuteronTPCTOF[i]    <<")"<< endl;
      //cout <<"Triton rejection with TPC preselection and PID cut at "                 << mMeanQRange[i] <<":\t("<< RejectionTritonTPC[i]         <<"\t+- "<<eRejectionTritonTPC[i]         <<")"<< endl;
      //cout <<"Triton rejection with TPC&TOF preselection and PID cut at "             << mMeanQRange[i] <<":\t("<< RejectionTritonTPCTOF[i]      <<"\t+- "<<eRejectionTritonTPCTOF[i]      <<")"<< endl;
      //cout <<"Helium3 rejection with TPC preselection and PID cut at "                << mMeanQRange[i] <<":\t("<< RejectionHelium3TPC[i]        <<"\t+- "<<eRejectionHelium3TPC[i]        <<")"<< endl;
      //cout <<"Helium3 rejection with TPC&TOF preselection and PID cut at "            << mMeanQRange[i] <<":\t("<< RejectionHelium3TPCTOF[i]     <<"\t+- "<<eRejectionHelium3TPCTOF[i]     <<")"<< endl;
      //cout <<"Alpha rejection with TPC preselection and PID cut at "                  << mMeanQRange[i] <<":\t("<< RejectionAlphaTPC[i]          <<"\t+- "<<eRejectionAlphaTPC[i]          <<")"<< endl;
      //cout <<"Alpha rejection with TPC&TOF preselection and PID cut at "              << mMeanQRange[i] <<":\t("<< RejectionAlphaTPCTOF[i]       <<"\t+- "<<eRejectionAlphaTPCTOF[i]       <<")"<< endl;
      //cout <<"Z=2 particle rejection with TPC preselection and PID cut at "           << mMeanQRange[i] <<":\t("<< RejectionHelium3AlphaTPC[i]   <<"\t+- "<<eRejectionHelium3AlphaTPC[i]   <<")"<< endl;
      //cout <<"Z=2 particle rejection with TPC&TOF preselection and PID cut at "       << mMeanQRange[i] <<":\t("<< RejectionHelium3AlphaTPCTOF[i]<<"\t+- "<<eRejectionHelium3AlphaTPCTOF[i]<<")"<< endl;
      //cout <<"d & t & He3 & Alpha rejection with TPC preselection and PID cut at "    << mMeanQRange[i] <<":\t("<< RejectionDeuteronTPC[i]       <<"\t+- "<<eRejectionDeuTriHe3AlpTPC[i]   <<")"<< endl;
      //cout <<"d & t & He3 & Alpha rejection with TPC&TOF preselection and PID cut at "<< mMeanQRange[i] <<":\t("<< RejectionDeuteronTPCTOF[i]    <<"\t+- "<<eRejectionDeuTriHe3AlpTPCTOF[i]<<")"<< endl;
    }


    TGraphErrors *fRejectionDeuteronTPC = new TGraphErrors(m,mMeanQRange,RejectionDeuteronTPC,errx,eRejectionDeuteronTPC);
    if(AllOrAnti == "All") fRejectionDeuteronTPC->SetTitle("Deuteron Rejection (TPC ps);<Q_{s}> (a.u.);Rejection");
    if(AllOrAnti == "Anti")fRejectionDeuteronTPC->SetTitle("Anti-Deuteron Rejection (TPC ps);<Q_{s}> (a.u.);Rejection");
    fRejectionDeuteronTPC->SetMarkerStyle(7);
    fRejectionDeuteronTPC->SetMarkerColor(4);
    TGraphErrors *fRejectionDeuteronTPCTOF = new TGraphErrors(m,mMeanQRange,RejectionDeuteronTPCTOF,errx,eRejectionDeuteronTPCTOF);
    if(AllOrAnti == "All") fRejectionDeuteronTPCTOF->SetTitle("Deuteron Rejection (TPC&TOF ps);<Q_{s}> (a.u.);Rejection");
    if(AllOrAnti == "Anti")fRejectionDeuteronTPCTOF->SetTitle("Anti-Deuteron Rejection (TPC&TOF ps);<Q_{s}> (a.u.);Rejection");
    fRejectionDeuteronTPCTOF->SetMarkerStyle(7);
    fRejectionDeuteronTPCTOF->SetMarkerColor(4);

    TGraphErrors *fRejectionTritonTPC = new TGraphErrors(m,mMeanQRange,RejectionTritonTPC,errx,eRejectionTritonTPC);
    if(AllOrAnti == "All") fRejectionTritonTPC->SetTitle("Triton Rejection (TPC ps);<Q_{s}> (a.u.);Rejection");
    if(AllOrAnti == "Anti")fRejectionTritonTPC->SetTitle("Anti-Triton Rejection (TPC ps);<Q_{s}> (a.u.);Rejection");
    fRejectionTritonTPC->SetMarkerStyle(7);
    fRejectionTritonTPC->SetMarkerColor(4);
    TGraphErrors *fRejectionTritonTPCTOF = new TGraphErrors(m,mMeanQRange,RejectionTritonTPCTOF,errx,eRejectionTritonTPCTOF);
    if(AllOrAnti == "All") fRejectionTritonTPCTOF->SetTitle("Triton Rejection (TPC&TOF ps);<Q_{s}> (a.u.);Rejection");
    if(AllOrAnti == "Anti")fRejectionTritonTPCTOF->SetTitle("Anti-Triton Rejection (TPC&TOF ps);<Q_{s}> (a.u.);Rejection");
    fRejectionTritonTPCTOF->SetMarkerStyle(7);
    fRejectionTritonTPCTOF->SetMarkerColor(4);

    TGraphErrors *fRejectionHelium3TPC = new TGraphErrors(m,mMeanQRange,RejectionHelium3TPC,errx,eRejectionHelium3TPC);
    if(AllOrAnti == "All") fRejectionHelium3TPC->SetTitle("Helium3 Rejection (TPC ps);<Q_{s}> (a.u.);Rejection");
    if(AllOrAnti == "Anti")fRejectionHelium3TPC->SetTitle("Anti-Helium3 Rejection (TPC ps);<Q_{s}> (a.u.);Rejection");
    fRejectionHelium3TPC->SetMarkerStyle(7);
    fRejectionHelium3TPC->SetMarkerColor(4);
    TGraphErrors *fRejectionHelium3TPCTOF = new TGraphErrors(m,mMeanQRange,RejectionHelium3TPCTOF,errx,eRejectionHelium3TPCTOF);
    if(AllOrAnti == "All") fRejectionHelium3TPCTOF->SetTitle("Helium3 Rejection (TPC&TOF ps);<Q_{s}> (a.u.);Rejection");
    if(AllOrAnti == "Anti")fRejectionHelium3TPCTOF->SetTitle("Anti-Helium3 Rejection (TPC&TOF ps);<Q_{s}> (a.u.);Rejection");
    fRejectionHelium3TPCTOF->SetMarkerStyle(7);
    fRejectionHelium3TPCTOF->SetMarkerColor(4);

    TGraphErrors *fRejectionAlphaTPC = new TGraphErrors(m,mMeanQRange,RejectionAlphaTPC,errx,eRejectionAlphaTPC);
    if(AllOrAnti == "All") fRejectionAlphaTPC->SetTitle("Alpha Rejection (TPC ps);<Q_{s}> (a.u.);Rejection");
    if(AllOrAnti == "Anti")fRejectionAlphaTPC->SetTitle("Anti-Alpha Rejection (TPC ps);<Q_{s}> (a.u.);Rejection");
    fRejectionAlphaTPC->SetMarkerStyle(7);
    fRejectionAlphaTPC->SetMarkerColor(4);
    TGraphErrors *fRejectionAlphaTPCTOF = new TGraphErrors(m,mMeanQRange,RejectionAlphaTPCTOF,errx,eRejectionAlphaTPCTOF);
    if(AllOrAnti == "All") fRejectionAlphaTPCTOF->SetTitle("Alpha Rejection (TPC&TOF ps);<Q_{s}> (a.u.);Rejection");
    if(AllOrAnti == "Anti")fRejectionAlphaTPCTOF->SetTitle("Anti-Alpha Rejection (TPC&TOF ps);<Q_{s}> (a.u.);Rejection");
    fRejectionAlphaTPCTOF->SetMarkerStyle(7);
    fRejectionAlphaTPCTOF->SetMarkerColor(4);

    TGraphErrors *fRejectionHelium3AlphaTPC = new TGraphErrors(m,mMeanQRange,RejectionHelium3AlphaTPC,errx,eRejectionHelium3AlphaTPC);
    if(AllOrAnti == "All") fRejectionHelium3AlphaTPC->SetTitle("Z=2 Particle Rejection (TPC ps);<Q_{s}> (a.u.);Rejection");
    if(AllOrAnti == "Anti")fRejectionHelium3AlphaTPC->SetTitle("Z=2 Anti-Particle Rejection (TPC ps);<Q_{s}> (a.u.);Rejection");
    fRejectionHelium3AlphaTPC->SetMarkerStyle(7);
    fRejectionHelium3AlphaTPC->SetMarkerColor(4);
    TGraphErrors *fRejectionHelium3AlphaTPCTOF = new TGraphErrors(m,mMeanQRange,RejectionHelium3AlphaTPCTOF,errx,eRejectionHelium3AlphaTPCTOF);
    if(AllOrAnti == "All") fRejectionHelium3AlphaTPCTOF->SetTitle("Z=2 Particle Rejection (TPC&TOF ps);<Q_{s}> (a.u.);Rejection");
    if(AllOrAnti == "Anti")fRejectionHelium3AlphaTPCTOF->SetTitle("Z=2 Anti-Particle Rejection (TPC&TOF ps);<Q_{s}> (a.u.);Rejection");
    fRejectionHelium3AlphaTPCTOF->SetMarkerStyle(7);
    fRejectionHelium3AlphaTPCTOF->SetMarkerColor(4);

    TGraphErrors *fRejectionDeuTriHe3AlpTPC = new TGraphErrors(m,mMeanQRange,RejectionDeuTriHe3AlpTPC,errx,eRejectionDeuTriHe3AlpTPC);
    if(AllOrAnti == "All") fRejectionDeuTriHe3AlpTPC->SetTitle("d & t & He3 & Alpha Rejection (TPC ps);<Q_{s}> (a.u.);Rejection");
    if(AllOrAnti == "Anti")fRejectionDeuTriHe3AlpTPC->SetTitle("Anti-d & -t & -He3 & -Alpha Rejection (TPC ps);<Q_{s}> (a.u.);Rejection");
    fRejectionDeuTriHe3AlpTPC->SetMarkerStyle(7);
    fRejectionDeuTriHe3AlpTPC->SetMarkerColor(4);
    TGraphErrors *fRejectionDeuTriHe3AlpTPCTOF = new TGraphErrors(m,mMeanQRange,RejectionDeuTriHe3AlpTPCTOF,errx,eRejectionDeuTriHe3AlpTPCTOF);
    if(AllOrAnti == "All") fRejectionDeuTriHe3AlpTPCTOF->SetTitle("d & t & He3 & Alpha Rejection (TPC&TOF ps);<Q_{s}> (a.u.);Rejection");
    if(AllOrAnti == "Anti")fRejectionDeuTriHe3AlpTPCTOF->SetTitle("Anti-d & -t & -He3 & -Alpha Rejection (TPC&TOF ps);<Q_{s}> (a.u.);Rejection");
    fRejectionDeuTriHe3AlpTPCTOF->SetMarkerStyle(7);
    fRejectionDeuTriHe3AlpTPCTOF->SetMarkerColor(4);



  // efficiency vs. rejection
    TGraphErrors *fEffRejDeuteronTPC = new TGraphErrors(m,EfficiencyDeuteronTPC,RejectionDeuteronTPC,eEfficiencyDeuteronTPC,eRejectionDeuteronTPC);
    if(AllOrAnti == "All") fEffRejDeuteronTPC->SetTitle("Deuteron Efficiency vs. Rejection (TPC ps);Efficiency;Rejection");
    if(AllOrAnti == "Anti")fEffRejDeuteronTPC->SetTitle("Anti-Deuteron Efficiency vs. Rejection (TPC ps);Efficiency;Rejection");
    fEffRejDeuteronTPC->SetMarkerStyle(7);
    fEffRejDeuteronTPC->SetMarkerColor(4);
    TGraphErrors *fEffRejDeuteronTPCTOF = new TGraphErrors(m,EfficiencyDeuteronTPCTOF,RejectionDeuteronTPCTOF,eEfficiencyDeuteronTPC,eRejectionDeuteronTPCTOF);
    if(AllOrAnti == "All") fEffRejDeuteronTPCTOF->SetTitle("Deuteron Efficiency vs. Rejection (TPC&TOF ps);Efficiency;Rejection");
    if(AllOrAnti == "Anti")fEffRejDeuteronTPCTOF->SetTitle("Anti-Deuteron Efficiency vs. Rejection (TPC&TOF ps);Efficiency;Rejection");
    fEffRejDeuteronTPCTOF->SetMarkerStyle(7);
    fEffRejDeuteronTPCTOF->SetMarkerColor(4);

    TGraphErrors *fEffRejTritonTPC = new TGraphErrors(m,EfficiencyTritonTPC,RejectionTritonTPC,eEfficiencyTritonTPC,eRejectionTritonTPC);
    if(AllOrAnti == "All") fEffRejTritonTPC->SetTitle("Triton Efficiency vs. Rejection (TPC ps);Efficiency;Rejection");
    if(AllOrAnti == "Anti")fEffRejTritonTPC->SetTitle("Anti-Triton Efficiency vs. Rejection (TPC ps);Efficiency;Rejection");
    fEffRejTritonTPC->SetMarkerStyle(7);
    fEffRejTritonTPC->SetMarkerColor(4);
    TGraphErrors *fEffRejTritonTPCTOF = new TGraphErrors(m,EfficiencyTritonTPCTOF,RejectionTritonTPCTOF,eEfficiencyTritonTPC,eRejectionTritonTPCTOF);
    if(AllOrAnti == "All") fEffRejTritonTPCTOF->SetTitle("Triton Efficiency vs. Rejection (TPC&TOF ps);Efficiency;Rejection");
    if(AllOrAnti == "Anti")fEffRejTritonTPCTOF->SetTitle("Anti-Triton Efficiency vs. Rejection (TPC&TOF ps);Efficiency;Rejection");
    fEffRejTritonTPCTOF->SetMarkerStyle(7);
    fEffRejTritonTPCTOF->SetMarkerColor(4);

    TGraphErrors *fEffRejHelium3TPC = new TGraphErrors(m,EfficiencyHelium3TPC,RejectionHelium3TPC,eEfficiencyHelium3TPC,eRejectionHelium3TPC);
    if(AllOrAnti == "All") fEffRejHelium3TPC->SetTitle("Helium3 Efficiency vs. Rejection (TPC ps);Efficiency;Rejection");
    if(AllOrAnti == "Anti")fEffRejHelium3TPC->SetTitle("Anti-Helium3 Efficiency vs. Rejection (TPC ps);Efficiency;Rejection");
    fEffRejHelium3TPC->SetMarkerStyle(7);
    fEffRejHelium3TPC->SetMarkerColor(4);
    TGraphErrors *fEffRejHelium3TPCTOF = new TGraphErrors(m,EfficiencyHelium3TPCTOF,RejectionHelium3TPCTOF,eEfficiencyHelium3TPC,eRejectionHelium3TPCTOF);
    if(AllOrAnti == "All") fEffRejHelium3TPCTOF->SetTitle("Helium3 Efficiency vs. Rejection (TPC&TOF ps);Efficiency;Rejection");
    if(AllOrAnti == "Anti")fEffRejHelium3TPCTOF->SetTitle("Anti-Helium3 Efficiency vs. Rejection (TPC&TOF ps);Efficiency;Rejection");
    fEffRejHelium3TPCTOF->SetMarkerStyle(7);
    fEffRejHelium3TPCTOF->SetMarkerColor(4);

    TGraphErrors *fEffRejAlphaTPC = new TGraphErrors(m,EfficiencyAlphaTPC,RejectionAlphaTPC,eEfficiencyAlphaTPC,eRejectionAlphaTPC);
    if(AllOrAnti == "All") fEffRejAlphaTPC->SetTitle("Alpha Efficiency vs. Rejection (TPC ps);Efficiency;Rejection");
    if(AllOrAnti == "Anti")fEffRejAlphaTPC->SetTitle("Anti-Alpha Efficiency vs. Rejection (TPC ps);Efficiency;Rejection");
    fEffRejAlphaTPC->SetMarkerStyle(7);
    fEffRejAlphaTPC->SetMarkerColor(4);
    TGraphErrors *fEffRejAlphaTPCTOF = new TGraphErrors(m,EfficiencyAlphaTPCTOF,RejectionAlphaTPCTOF,eEfficiencyAlphaTPC,eRejectionAlphaTPCTOF);
    if(AllOrAnti == "All") fEffRejAlphaTPCTOF->SetTitle("Alpha Efficiency vs. Rejection (TPC&TOF ps);Efficiency;Rejection");
    if(AllOrAnti == "Anti")fEffRejAlphaTPCTOF->SetTitle("Anti-Alpha Efficiency vs. Rejection (TPC&TOF ps);Efficiency;Rejection");
    fEffRejAlphaTPCTOF->SetMarkerStyle(7);
    fEffRejAlphaTPCTOF->SetMarkerColor(4);

    TGraphErrors *fEffRejHelium3AlphaTPC = new TGraphErrors(m,EfficiencyHelium3AlphaTPC,RejectionHelium3AlphaTPC,eEfficiencyHelium3AlphaTPC,eRejectionHelium3AlphaTPC);
    if(AllOrAnti == "All") fEffRejHelium3AlphaTPC->SetTitle("Z=2 Particle Efficiency vs. Rejection (TPC ps);Efficiency;Rejection");
    if(AllOrAnti == "Anti")fEffRejHelium3AlphaTPC->SetTitle("Z=2 Anti-Particle Efficiency vs. Rejection (TPC ps);Efficiency;Rejection");
    fEffRejHelium3AlphaTPC->SetMarkerStyle(7);
    fEffRejHelium3AlphaTPC->SetMarkerColor(4);
    TGraphErrors *fEffRejHelium3AlphaTPCTOF = new TGraphErrors(m,EfficiencyHelium3AlphaTPCTOF,RejectionHelium3AlphaTPCTOF,eEfficiencyHelium3AlphaTPC,eRejectionHelium3AlphaTPCTOF);
    if(AllOrAnti == "All") fEffRejHelium3AlphaTPCTOF->SetTitle("Z=2 Particle Efficiency vs. Rejection (TPC&TOF ps);Efficiency;Rejection");
    if(AllOrAnti == "Anti")fEffRejHelium3AlphaTPCTOF->SetTitle("Z=2 Anti-Particle Efficiency vs. Rejection (TPC&TOF ps);Efficiency;Rejection");
    fEffRejHelium3AlphaTPCTOF->SetMarkerStyle(7);
    fEffRejHelium3AlphaTPCTOF->SetMarkerColor(4);

    TGraphErrors *fEffRejDeuTriHe3AlpTPC = new TGraphErrors(m,EfficiencyDeuTriHe3AlpTPC,RejectionDeuTriHe3AlpTPC,eEfficiencyDeuTriHe3AlpTPC,eRejectionDeuTriHe3AlpTPC);
    if(AllOrAnti == "All") fEffRejDeuTriHe3AlpTPC->SetTitle("d & t & He3 & Alpha Efficiency vs. Rejection (TPC ps);Efficiency;Rejection");
    if(AllOrAnti == "Anti")fEffRejDeuTriHe3AlpTPC->SetTitle("Anti-d & -t & -He3 & -Alpha Efficiency vs. Rejection (TPC ps);Efficiency;Rejection");
    fEffRejDeuTriHe3AlpTPC->SetMarkerStyle(7);
    fEffRejDeuTriHe3AlpTPC->SetMarkerColor(4);
    TGraphErrors *fEffRejDeuTriHe3AlpTPCTOF = new TGraphErrors(m,EfficiencyDeuTriHe3AlpTPCTOF,RejectionDeuTriHe3AlpTPCTOF,eEfficiencyDeuTriHe3AlpTPC,eRejectionDeuTriHe3AlpTPCTOF);
    if(AllOrAnti == "All") fEffRejDeuTriHe3AlpTPCTOF->SetTitle("d & t & He3 & Alpha Efficiency vs. Rejection (TPC&TOF ps);Efficiency;Rejection");
    if(AllOrAnti == "Anti")fEffRejDeuTriHe3AlpTPCTOF->SetTitle("Anti-d & -t & -He3 & -Alpha Efficiency vs. Rejection (TPC&TOF ps);Efficiency;Rejection");
    fEffRejDeuTriHe3AlpTPCTOF->SetMarkerStyle(7);
    fEffRejDeuTriHe3AlpTPCTOF->SetMarkerColor(4);



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Create legends
  ///_________________
  // legends for cd(1)
  TLegend *legPID = new TLegend(.02,.68,.45,.88);
  legPID->SetFillStyle(0);
  legPID->SetBorderSize(0);
  if(Collisions == "PbPb_11h" || Collisions == "PbPb_11hStd") {
    legPID->AddEntry((TObject*)0,"Pb-Pb, #sqrt{s_{NN}} = 2.76 TeV","");
    legPID->AddEntry((TObject*)0,"run LHC11h pass2","");
  }
  else if(Collisions == "pPb_13b") {
    legPID->AddEntry((TObject*)0,"p-Pb, #sqrt{s_{NN}} = 5.023 TeV","");
    legPID->AddEntry((TObject*)0,"run LHC13b pass3","");
  }
  else if(Collisions == "pPb_13c") {
    legPID->AddEntry((TObject*)0,"p-Pb, #sqrt{s_{NN}} = 5.023 TeV","");
    legPID->AddEntry((TObject*)0,"run LHC13c pass2","");
  }
  else { //(Collisions == "pp_15f")
    legPID->AddEntry((TObject*)0,"p-p, #sqrt{s}} = 13 TeV","");
    legPID->AddEntry((TObject*)0,"run LHC15f pass2","");
  }

  // legends for cd(2)
  TLegend *legPIDDeuteronTPC = new TLegend(.05,.68,.4,.85);
  legPIDDeuteronTPC->SetFillStyle(0);
  legPIDDeuteronTPC->SetBorderSize(0);
  legPIDDeuteronTPC->AddEntry((TObject*)0,"with preselection:","");
  legPIDDeuteronTPC->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legPIDDeuteronTPC->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutDeuteronTPC),"");
  TLegend *legPIDDeuteronTPCTOF = new TLegend(.05,.63,.4,.85);
  legPIDDeuteronTPCTOF->SetFillStyle(0);
  legPIDDeuteronTPCTOF->SetBorderSize(0);
  legPIDDeuteronTPCTOF->AddEntry((TObject*)0,"with preselection:","");
  legPIDDeuteronTPCTOF->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legPIDDeuteronTPCTOF->AddEntry((TObject*)0,Form("|TOFnSigma| = %i#sigma",Sigma),"");
  legPIDDeuteronTPCTOF->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutDeuteronTPCTOF),"");
  TLegend *legPIDTritonTPC = new TLegend(.05,.68,.4,.85);
  legPIDTritonTPC->SetFillStyle(0);
  legPIDTritonTPC->SetBorderSize(0);
  legPIDTritonTPC->AddEntry((TObject*)0,"with preselection:","");
  legPIDTritonTPC->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legPIDTritonTPC->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutTritonTPC),"");
  TLegend *legPIDTritonTPCTOF = new TLegend(.05,.63,.4,.85);
  legPIDTritonTPCTOF->SetFillStyle(0);
  legPIDTritonTPCTOF->SetBorderSize(0);
  legPIDTritonTPCTOF->AddEntry((TObject*)0,"with preselection:","");
  legPIDTritonTPCTOF->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legPIDTritonTPCTOF->AddEntry((TObject*)0,Form("|TOFnSigma| = %i#sigma",Sigma),"");
  legPIDTritonTPCTOF->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutTritonTPCTOF),"");
  TLegend *legPIDDeuTriHe3AlpTPC = new TLegend(.05,.63,.45,.85);
  legPIDDeuTriHe3AlpTPC->SetFillStyle(0);
  legPIDDeuTriHe3AlpTPC->SetBorderSize(0);
  legPIDDeuTriHe3AlpTPC->AddEntry((TObject*)0,"with preselection:","");
  legPIDDeuTriHe3AlpTPC->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legPIDDeuTriHe3AlpTPC->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c} for d",ClearCutDeuteronTPC),"");
  legPIDDeuTriHe3AlpTPC->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c} for t",ClearCutTritonTPC),"");
  TLegend *legPIDDeuTriHe3AlpTPCTOF = new TLegend(.05,.58,.45,.85);
  legPIDDeuTriHe3AlpTPCTOF->SetFillStyle(0);
  legPIDDeuTriHe3AlpTPCTOF->SetBorderSize(0);
  legPIDDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,"with preselection:","");
  legPIDDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legPIDDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,Form("|TOFnSigma| = %i#sigma",Sigma),"");
  legPIDDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c} for d",ClearCutDeuteronTPCTOF),"");
  legPIDDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c} for t",ClearCutTritonTPCTOF),"");
  TLegend *legPIDTPC = new TLegend(.05,.73,.4,.85);
  legPIDTPC->SetFillStyle(0);
  legPIDTPC->SetBorderSize(0);
  legPIDTPC->AddEntry((TObject*)0,"with preselection:","");
  legPIDTPC->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  TLegend *legPIDTPCTOF = new TLegend(.05,.67,.4,.85);
  legPIDTPCTOF->SetFillStyle(0);
  legPIDTPCTOF->SetBorderSize(0);
  legPIDTPCTOF->AddEntry((TObject*)0,"with preselection:","");
  legPIDTPCTOF->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legPIDTPCTOF->AddEntry((TObject*)0,Form("|TOFnSigma| = %i#sigma",Sigma),"");

  // legends for cd(3)
  TLegend *legProjDeuteronTPC = new TLegend(.45,.63,.8,.85);
  legProjDeuteronTPC->SetFillStyle(0);
  legProjDeuteronTPC->SetBorderSize(0);
  if(AllOrAnti == "All") legProjDeuteronTPC->AddEntry(fHistTRDPIDDeuteronTPC_clear_Clone,"Deuteron","l");
  if(AllOrAnti == "Anti")legProjDeuteronTPC->AddEntry(fHistTRDPIDDeuteronTPC_clear_Clone,"Anti-Deuteron","l");
  legProjDeuteronTPC->AddEntry((TObject*)0,"with preselection:","");
  legProjDeuteronTPC->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legProjDeuteronTPC->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutDeuteronTPC),"");
  TLegend *legProjDeuteronTPCTOF = new TLegend(.45,.58,.8,.85);
  legProjDeuteronTPCTOF->SetFillStyle(0);
  legProjDeuteronTPCTOF->SetBorderSize(0);
  if(AllOrAnti == "All") legProjDeuteronTPCTOF->AddEntry(fHistTRDPIDDeuteronTPCTOF_clear_Clone,"Deuteron","l");
  if(AllOrAnti == "Anti")legProjDeuteronTPCTOF->AddEntry(fHistTRDPIDDeuteronTPCTOF_clear_Clone,"Anti-Deuteron","l");
  legProjDeuteronTPCTOF->AddEntry((TObject*)0,"with preselection:","");
  legProjDeuteronTPCTOF->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legProjDeuteronTPCTOF->AddEntry((TObject*)0,Form("|TOFnSigma| = %i#sigma",Sigma),"");
  legProjDeuteronTPCTOF->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutDeuteronTPCTOF),"");
  TLegend *legProjTritonTPC = new TLegend(.55,.6,.9,.82);
  legProjTritonTPC->SetFillStyle(0);
  legProjTritonTPC->SetBorderSize(0);
  if(AllOrAnti == "All") legProjTritonTPC->AddEntry(fHistTRDPIDTritonTPC_clear_Clone,"Triton","l");
  if(AllOrAnti == "Anti")legProjTritonTPC->AddEntry(fHistTRDPIDTritonTPC_clear_Clone,"Anti-Triton","l");
  legProjTritonTPC->AddEntry((TObject*)0,"with preselection:","");
  legProjTritonTPC->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legProjTritonTPC->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutTritonTPC),"");
  TLegend *legProjTritonTPCTOF = new TLegend(.55,.55,.9,.82);
  legProjTritonTPCTOF->SetFillStyle(0);
  legProjTritonTPCTOF->SetBorderSize(0);
  if(AllOrAnti == "All") legProjTritonTPCTOF->AddEntry(fHistTRDPIDTritonTPCTOF_clear_Clone,"Triton","l");
  if(AllOrAnti == "Anti")legProjTritonTPCTOF->AddEntry(fHistTRDPIDTritonTPCTOF_clear_Clone,"Anti-Triton","l");
  legProjTritonTPCTOF->AddEntry((TObject*)0,"with preselection:","");
  legProjTritonTPCTOF->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legProjTritonTPCTOF->AddEntry((TObject*)0,Form("|TOFnSigma| = %i#sigma",Sigma),"");
  legProjTritonTPCTOF->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutTritonTPCTOF),"");
  TLegend *legProjHelium3TPC = new TLegend(.12,.67,.47,.85);
  legProjHelium3TPC->SetFillStyle(0);
  legProjHelium3TPC->SetBorderSize(0);
  if(AllOrAnti == "All") legProjHelium3TPC->AddEntry(fHistTRDPIDHelium3TPC_Clone,"Helium3","l");
  if(AllOrAnti == "Anti")legProjHelium3TPC->AddEntry(fHistTRDPIDHelium3TPC_Clone,"Anti-Helium3","l");
  legProjHelium3TPC->AddEntry((TObject*)0,"with preselection:","");
  legProjHelium3TPC->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  TLegend *legProjHelium3TPCTOF = new TLegend(.1,.61,.45,.85);
  legProjHelium3TPCTOF->SetFillStyle(0);
  legProjHelium3TPCTOF->SetBorderSize(0);
  if(AllOrAnti == "All") legProjHelium3TPCTOF->AddEntry(fHistTRDPIDHelium3TPCTOF_Clone,"Helium3","l");
  if(AllOrAnti == "Anti")legProjHelium3TPCTOF->AddEntry(fHistTRDPIDHelium3TPCTOF_Clone,"Anti-Helium3","l");
  legProjHelium3TPCTOF->AddEntry((TObject*)0,"with preselection:","");
  legProjHelium3TPCTOF->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legProjHelium3TPCTOF->AddEntry((TObject*)0,Form("|TOFnSigma| = %i#sigma",Sigma),"");
  TLegend *legProjAlphaTPC = new TLegend(.45,.67,.8,.85);
  legProjAlphaTPC->SetFillStyle(0);
  legProjAlphaTPC->SetBorderSize(0);
  if(AllOrAnti == "All") legProjAlphaTPC->AddEntry(fHistTRDPIDAlphaTPC_Clone,"Alpha","l");
  if(AllOrAnti == "Anti")legProjAlphaTPC->AddEntry(fHistTRDPIDAlphaTPC_Clone,"Anti-Alpha","l");
  legProjAlphaTPC->AddEntry((TObject*)0,"with preselection:","");
  legProjAlphaTPC->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  TLegend *legProjAlphaTPCTOF = new TLegend(.12,.61,.47,.85);
  legProjAlphaTPCTOF->SetFillStyle(0);
  legProjAlphaTPCTOF->SetBorderSize(0);
  if(AllOrAnti == "All") legProjAlphaTPCTOF->AddEntry(fHistTRDPIDAlphaTPCTOF_Clone,"Alpha","l");
  if(AllOrAnti == "Anti")legProjAlphaTPCTOF->AddEntry(fHistTRDPIDAlphaTPCTOF_Clone,"Anti-Alpha","l");
  legProjAlphaTPCTOF->AddEntry((TObject*)0,"with preselection:","");
  legProjAlphaTPCTOF->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legProjAlphaTPCTOF->AddEntry((TObject*)0,Form("|TOFnSigma| = %i#sigma",Sigma),"");
  TLegend *legProjHelium3AlphaTPC = new TLegend(.12,.67,.47,.85);
  legProjHelium3AlphaTPC->SetFillStyle(0);
  legProjHelium3AlphaTPC->SetBorderSize(0);
  if(AllOrAnti == "All") legProjHelium3AlphaTPC->AddEntry(fHistTRDPIDHelium3AlphaTPC_Clone,"Helium3 & Alpha","l");
  if(AllOrAnti == "Anti")legProjHelium3AlphaTPC->AddEntry(fHistTRDPIDHelium3AlphaTPC_Clone,"Anti-Helium3 & -Alpha","l");
  legProjHelium3AlphaTPC->AddEntry((TObject*)0,"with preselection:","");
  legProjHelium3AlphaTPC->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  TLegend *legProjHelium3AlphaTPCTOF = new TLegend(.1,.61,.45,.85);
  legProjHelium3AlphaTPCTOF->SetFillStyle(0);
  legProjHelium3AlphaTPCTOF->SetBorderSize(0);
  if(AllOrAnti == "All") legProjHelium3AlphaTPCTOF->AddEntry(fHistTRDPIDHelium3AlphaTPCTOF_Clone,"Helium3 & Alpha","l");
  if(AllOrAnti == "Anti")legProjHelium3AlphaTPCTOF->AddEntry(fHistTRDPIDHelium3AlphaTPCTOF_Clone,"Anti-Helium3 & -Alpha","l");
  legProjHelium3AlphaTPCTOF->AddEntry((TObject*)0,"with preselection:","");
  legProjHelium3AlphaTPCTOF->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legProjHelium3AlphaTPCTOF->AddEntry((TObject*)0,Form("|TOFnSigma| = %i#sigma",Sigma),"");
  TLegend *legProjDeuTriHe3AlpTPC = new TLegend(.38,.45,.78,.85);
  legProjDeuTriHe3AlpTPC->SetFillStyle(0);
  legProjDeuTriHe3AlpTPC->SetBorderSize(0);
  if(AllOrAnti == "All") legProjDeuTriHe3AlpTPC->AddEntry(fHistTRDPIDDeuTriHe3AlpTPC_Clone,"d & t & He3 & Alpha","l");
  if(AllOrAnti == "Anti")legProjDeuTriHe3AlpTPC->AddEntry(fHistTRDPIDDeuTriHe3AlpTPC_Clone,"#bar{d} & #bar{t} & {}^{3}#bar{He} & #bar{#alpha}","l");
  legProjDeuTriHe3AlpTPC->AddEntry((TObject*)0,"with preselection:","");
  legProjDeuTriHe3AlpTPC->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legProjDeuTriHe3AlpTPC->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c} for d",ClearCutDeuteronTPC),"");
  legProjDeuTriHe3AlpTPC->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c} for t",ClearCutTritonTPC),"");
  TLegend *legProjDeuTriHe3AlpTPCTOF = new TLegend(.38,.39,.78,.85);
  legProjDeuTriHe3AlpTPCTOF->SetFillStyle(0);
  legProjDeuTriHe3AlpTPCTOF->SetBorderSize(0);
  if(AllOrAnti == "All") legProjDeuTriHe3AlpTPCTOF->AddEntry(fHistTRDPIDDeuTriHe3AlpTPCTOF_Clone,"d & t & He3 & Alpha","l");
  if(AllOrAnti == "Anti")legProjDeuTriHe3AlpTPCTOF->AddEntry(fHistTRDPIDDeuTriHe3AlpTPCTOF_Clone,"#bar{d} & #bar{t} & {}^{3}#bar{He} & #bar{#alpha}","l");
  legProjDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,"with preselection:","");
  legProjDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legProjDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,Form("|TOFnSigma| = %i#sigma",Sigma),"");
  legProjDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c} for d",ClearCutDeuteronTPCTOF),"");
  legProjDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c} for t",ClearCutTritonTPCTOF),"");

  // legends for cd(4)
  TLegend *legMeanallDeuteronTPC = new TLegend(.15,.7,.45,.9);
  legMeanallDeuteronTPC->SetFillStyle(0);
  legMeanallDeuteronTPC->SetBorderSize(0);
  legMeanallDeuteronTPC->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDAllexceptDeuteronTPC_clear->GetMean()),"");
  if(AllOrAnti == "All") TLegend *legMeanDeuteronTPC = new TLegend(.17,.36,.47,.56);
  if(AllOrAnti == "Anti")TLegend *legMeanDeuteronTPC = new TLegend(.26,.36,.56,.56);
  legMeanDeuteronTPC->SetFillStyle(0);
  legMeanDeuteronTPC->SetBorderSize(0);
  legMeanDeuteronTPC->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDDeuteronTPC_clear_Clone2->GetMean()),"");
  TLegend *legMeanallDeuteronTPCTOF = new TLegend(.15,.7,.45,.9);
  legMeanallDeuteronTPCTOF->SetFillStyle(0);
  legMeanallDeuteronTPCTOF->SetBorderSize(0);
  legMeanallDeuteronTPCTOF->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDAllexceptDeuteronTPCTOF_clear->GetMean()),"");
  if(AllOrAnti == "All") TLegend *legMeanDeuteronTPCTOF = new TLegend(.15,.42,.45,.62);
  if(AllOrAnti == "Anti")TLegend *legMeanDeuteronTPCTOF = new TLegend(.26,.42,.56,.62);
  legMeanDeuteronTPCTOF->SetFillStyle(0);
  legMeanDeuteronTPCTOF->SetBorderSize(0);
  legMeanDeuteronTPCTOF->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDDeuteronTPCTOF_clear_Clone2->GetMean()),"");
  TLegend *legParticleDeuteron = new TLegend(.5,.6,.8,.77);
  legParticleDeuteron->SetFillStyle(0);
  legParticleDeuteron->SetBorderSize(0);
  if(AllOrAnti == "All") legParticleDeuteron->AddEntry(fHistTRDPIDDeuteronTPC_clear_Clone2,"d, #bar{d}","l");
  if(AllOrAnti == "Anti")legParticleDeuteron->AddEntry(fHistTRDPIDDeuteronTPC_clear_Clone2,"#bar{d}","l");
  legParticleDeuteron->AddEntry(fHistTRDPIDAllexceptDeuteronTPC_clear,"else","l");

  TLegend *legMeanallTritonTPC = new TLegend(.15,.7,.45,.9);
  legMeanallTritonTPC->SetFillStyle(0);
  legMeanallTritonTPC->SetBorderSize(0);
  legMeanallTritonTPC->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDAllexceptTritonTPC_clear->GetMean()),"");
  if(AllOrAnti == "All") TLegend *legMeanTritonTPC = new TLegend(.3,.24,.62,.44);
  if(AllOrAnti == "Anti")TLegend *legMeanTritonTPC = new TLegend(.34,.34,.66,.54);
  legMeanTritonTPC->SetFillStyle(0);
  legMeanTritonTPC->SetBorderSize(0);
  legMeanTritonTPC->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDTritonTPC_clear_Clone2->GetMean()),"");
  TLegend *legMeanallTritonTPCTOF = new TLegend(.15,.7,.45,.9);
  legMeanallTritonTPCTOF->SetFillStyle(0);
  legMeanallTritonTPCTOF->SetBorderSize(0);
  legMeanallTritonTPCTOF->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDAllexceptTritonTPCTOF_clear->GetMean()),"");
  TLegend *legMeanTritonTPCTOF = new TLegend(.27,.36,.59,.56);
  legMeanTritonTPCTOF->SetFillStyle(0);
  legMeanTritonTPCTOF->SetBorderSize(0);
  legMeanTritonTPCTOF->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDTritonTPCTOF_clear_Clone2->GetMean()),"");
  TLegend *legParticleTriton = new TLegend(.5,.6,.8,.77);
  legParticleTriton->SetFillStyle(0);
  legParticleTriton->SetBorderSize(0);
  if(AllOrAnti == "All") legParticleTriton->AddEntry(fHistTRDPIDTritonTPC_clear_Clone2,"t, #bar{t}","l");
  if(AllOrAnti == "Anti")legParticleTriton->AddEntry(fHistTRDPIDTritonTPC_clear_Clone2,"#bar{t}","l");
  legParticleTriton->AddEntry(fHistTRDPIDAllexceptTritonTPC_clear,"else","l");

  TLegend *legMeanallHelium3TPC = new TLegend(.15,.7,.45,.9);
  legMeanallHelium3TPC->SetFillStyle(0);
  legMeanallHelium3TPC->SetBorderSize(0);
  legMeanallHelium3TPC->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDAllexceptHelium3TPC->GetMean()),"");
  if(AllOrAnti == "All") TLegend *legMeanHelium3TPC = new TLegend(.43,.31,.75,.52);
  if(AllOrAnti == "Anti")TLegend *legMeanHelium3TPC = new TLegend(.42,.35,.74,.56);
  legMeanHelium3TPC->SetFillStyle(0);
  legMeanHelium3TPC->SetBorderSize(0);
  legMeanHelium3TPC->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDHelium3TPC_Clone2->GetMean()),"");
  TLegend *legMeanallHelium3TPCTOF = new TLegend(.15,.7,.45,.9);
  legMeanallHelium3TPCTOF->SetFillStyle(0);
  legMeanallHelium3TPCTOF->SetBorderSize(0);
  legMeanallHelium3TPCTOF->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDAllexceptHelium3TPCTOF->GetMean()),"");
  TLegend *legMeanHelium3TPCTOF = new TLegend(.43,.34,.75,.55);
  legMeanHelium3TPCTOF->SetFillStyle(0);
  legMeanHelium3TPCTOF->SetBorderSize(0);
  legMeanHelium3TPCTOF->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDHelium3TPCTOF_Clone2->GetMean()),"");
  TLegend *legParticleHelium3 = new TLegend(.5,.6,.8,.77);
  legParticleHelium3->SetFillStyle(0);
  legParticleHelium3->SetBorderSize(0);
  if(AllOrAnti == "All") legParticleHelium3->AddEntry(fHistTRDPIDHelium3TPC_Clone2,"{}^{3}He,{}^{3}#bar{He}","l");
  if(AllOrAnti == "Anti")legParticleHelium3->AddEntry(fHistTRDPIDHelium3TPC_Clone2,"{}^{3}#bar{He}","l");
  legParticleHelium3->AddEntry(fHistTRDPIDAllexceptHelium3TPC,"else","l");

  TLegend *legMeanallAlphaTPC = new TLegend(.15,.7,.45,.9);
  legMeanallAlphaTPC->SetFillStyle(0);
  legMeanallAlphaTPC->SetBorderSize(0);
  legMeanallAlphaTPC->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDAllexceptAlphaTPC->GetMean()),"");
  TLegend *legMeanAlphaTPC = new TLegend(.32,.3,.64,.5);
  legMeanAlphaTPC->SetFillStyle(0);
  legMeanAlphaTPC->SetBorderSize(0);
  legMeanAlphaTPC->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDAlphaTPC_Clone2->GetMean()),"");
  TLegend *legMeanallAlphaTPCTOF = new TLegend(.15,.7,.45,.9);
  legMeanallAlphaTPCTOF->SetFillStyle(0);
  legMeanallAlphaTPCTOF->SetBorderSize(0);
  legMeanallAlphaTPCTOF->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDAllexceptAlphaTPCTOF->GetMean()),"");
  TLegend *legMeanAlphaTPCTOF = new TLegend(.54,.27,.86,.47);
  legMeanAlphaTPCTOF->SetFillStyle(0);
  legMeanAlphaTPCTOF->SetBorderSize(0);
  legMeanAlphaTPCTOF->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDAlphaTPCTOF_Clone2->GetMean()),"");
  TLegend *legParticleAlpha = new TLegend(.5,.6,.8,.77);
  legParticleAlpha->SetFillStyle(0);
  legParticleAlpha->SetBorderSize(0);
  if(AllOrAnti == "All") legParticleAlpha->AddEntry(fHistTRDPIDAlphaTPC_Clone2,"#alpha, #bar{#alpha}","l");
  if(AllOrAnti == "Anti")legParticleAlpha->AddEntry(fHistTRDPIDAlphaTPC_Clone2,"#bar{#alpha}","l");
  legParticleAlpha->AddEntry(fHistTRDPIDAllexceptAlphaTPC,"else","l");

  TLegend *legMeanallHelium3AlphaTPC = new TLegend(.15,.7,.45,.9);
  legMeanallHelium3AlphaTPC->SetFillStyle(0);
  legMeanallHelium3AlphaTPC->SetBorderSize(0);
  legMeanallHelium3AlphaTPC->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDAllexceptHelium3AlphaTPC->GetMean()),"");
  TLegend *legMeanHelium3AlphaTPC = new TLegend(.42,.32,.74,.52);
  legMeanHelium3AlphaTPC->SetFillStyle(0);
  legMeanHelium3AlphaTPC->SetBorderSize(0);
  legMeanHelium3AlphaTPC->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDHelium3AlphaTPC->GetMean()),"");
  TLegend *legMeanallHelium3AlphaTPCTOF = new TLegend(.15,.7,.45,.9);
  legMeanallHelium3AlphaTPCTOF->SetFillStyle(0);
  legMeanallHelium3AlphaTPCTOF->SetBorderSize(0);
  legMeanallHelium3AlphaTPCTOF->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDAllexceptHelium3AlphaTPCTOF->GetMean()),"");
  TLegend *legMeanHelium3AlphaTPCTOF = new TLegend(.44,.34,.76,.54);
  legMeanHelium3AlphaTPCTOF->SetFillStyle(0);
  legMeanHelium3AlphaTPCTOF->SetBorderSize(0);
  legMeanHelium3AlphaTPCTOF->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDHelium3AlphaTPCTOF->GetMean()),"");
  TLegend *legParticleHelium3Alpha = new TLegend(.5,.6,.85,.77);
  legParticleHelium3Alpha->SetFillStyle(0);
  legParticleHelium3Alpha->SetBorderSize(0);
  if(AllOrAnti == "All") legParticleHelium3Alpha->AddEntry(fHistTRDPIDHelium3AlphaTPC,"{}^{3}He,{}^{3}#bar{He}, #alpha, #bar{#alpha}","l");
  if(AllOrAnti == "Anti")legParticleHelium3Alpha->AddEntry(fHistTRDPIDHelium3AlphaTPC,"{}^{3}#bar{He}, #bar{#alpha}","l");
  legParticleHelium3Alpha->AddEntry(fHistTRDPIDAllexceptHelium3AlphaTPC,"else","l");

  TLegend *legMeanallDeuTriHe3AlpTPC = new TLegend(.15,.7,.45,.9);
  legMeanallDeuTriHe3AlpTPC->SetFillStyle(0);
  legMeanallDeuTriHe3AlpTPC->SetBorderSize(0);
  legMeanallDeuTriHe3AlpTPC->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDAllexceptDeuTriHe3AlpTPC->GetMean()),"");
  if(AllOrAnti == "All") TLegend *legMeanDeuTriHe3AlpTPC = new TLegend(.17,.36,.47,.56);
  if(AllOrAnti == "Anti")TLegend *legMeanDeuTriHe3AlpTPC = new TLegend(.32,.21,.62,.41);
  legMeanDeuTriHe3AlpTPC->SetFillStyle(0);
  legMeanDeuTriHe3AlpTPC->SetBorderSize(0);
  legMeanDeuTriHe3AlpTPC->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDDeuTriHe3AlpTPC->GetMean()),"");
  TLegend *legMeanallDeuTriHe3AlpTPCTOF = new TLegend(.15,.7,.45,.9);
  legMeanallDeuTriHe3AlpTPCTOF->SetFillStyle(0);
  legMeanallDeuTriHe3AlpTPCTOF->SetBorderSize(0);
  legMeanallDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDAllexceptDeuTriHe3AlpTPCTOF->GetMean()),"");
  if(AllOrAnti == "All") TLegend *legMeanDeuTriHe3AlpTPCTOF = new TLegend(.15,.44,.45,.59);
  if(AllOrAnti == "Anti")TLegend *legMeanDeuTriHe3AlpTPCTOF = new TLegend(.32,.21,.62,.41);
  legMeanDeuTriHe3AlpTPCTOF->SetFillStyle(0);
  legMeanDeuTriHe3AlpTPCTOF->SetBorderSize(0);
  legMeanDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,Form("Mean = %.2f",fHistTRDPIDDeuTriHe3AlpTPCTOF->GetMean()),"");
  TLegend *legParticleDeuTriHe3Alp = new TLegend(.5,.45,.85,.82);
  legParticleDeuTriHe3Alp->SetFillStyle(0);
  legParticleDeuTriHe3Alp->SetBorderSize(0);
  if(AllOrAnti == "All") legParticleDeuTriHe3Alp->AddEntry(fHistTRDPIDDeuTriHe3AlpTPC,"#splitline{d, #bar{d}, t, #bar{t},}{{}^{3}He,{}^{3}#bar{He}, #alpha, #bar{#alpha}}","l");
  if(AllOrAnti == "Anti")legParticleDeuTriHe3Alp->AddEntry(fHistTRDPIDDeuTriHe3AlpTPC,"#bar{d}, #bar{t},{}^{3}#bar{He}, #bar{#alpha}","l");
  legParticleDeuTriHe3Alp->AddEntry(fHistTRDPIDAllexceptDeuTriHe3AlpTPC,"else","l");

  // legends for cd(5)
  TLegend *legPIDCutPurDeuteronTPC = new TLegend(.01,y1,.54,.88);
  legPIDCutPurDeuteronTPC->SetFillStyle(0);
  legPIDCutPurDeuteronTPC->SetBorderSize(0);
  legPIDCutPurDeuteronTPC->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutPurDeuteronTPC->AddEntry((TObject*)0,"Deuteron purity for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutPurDeuteronTPC->AddEntry((TObject*)0,Form("#bullet %i:  %.5f #pm %.5f",Deu + i*10,PurityDeuteronTPC[nDeu + i],ePurityDeuteronTPC[nDeu + i]),"");
  }
  TLegend *legPIDCutPurDeuteronTPCTOF = new TLegend(.35,y1,.88,.88);
  legPIDCutPurDeuteronTPCTOF->SetFillStyle(0);
  legPIDCutPurDeuteronTPCTOF->SetBorderSize(0);
  legPIDCutPurDeuteronTPCTOF->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutPurDeuteronTPCTOF->AddEntry((TObject*)0,"Deuteron purity for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutPurDeuteronTPCTOF->AddEntry((TObject*)0,Form("#bullet %i:  %.5f #pm %.5f",Deu + i*10,PurityDeuteronTPCTOF[nDeu + i],ePurityDeuteronTPCTOF[nDeu + i]),"");
  }
  TLegend *legPIDCutPurTritonTPC = new TLegend(.01,y1,.54,.88);
  legPIDCutPurTritonTPC->SetFillStyle(0);
  legPIDCutPurTritonTPC->SetBorderSize(0);
  legPIDCutPurTritonTPC->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutPurTritonTPC->AddEntry((TObject*)0,"Triton purity for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutPurTritonTPC->AddEntry((TObject*)0,Form("#bullet %i:  %.5f #pm %.5f",Tri + i*10,PurityTritonTPC[nTri + i],ePurityTritonTPC[nTri + i]),"");
  }
  TLegend *legPIDCutPurTritonTPCTOF = new TLegend(.01,y1,.54,.88);
  legPIDCutPurTritonTPCTOF->SetFillStyle(0);
  legPIDCutPurTritonTPCTOF->SetBorderSize(0);
  legPIDCutPurTritonTPCTOF->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutPurTritonTPCTOF->AddEntry((TObject*)0,"Triton purity for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutPurTritonTPCTOF->AddEntry((TObject*)0,Form("#bullet %i:  %.5f #pm %.5f",Tri + i*10,PurityTritonTPCTOF[nTri + i],ePurityTritonTPCTOF[nTri + i]),"");
  }
  TLegend *legPIDCutPurHelium3TPC = new TLegend(.01,y1,.54,.88);
  legPIDCutPurHelium3TPC->SetFillStyle(0);
  legPIDCutPurHelium3TPC->SetBorderSize(0);
  legPIDCutPurHelium3TPC->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutPurHelium3TPC->AddEntry((TObject*)0,"Helium3 purity for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutPurHelium3TPC->AddEntry((TObject*)0,Form("#bullet %i:  %.6f #pm %.6f",He3 + i*10,PurityHelium3TPC[nHe3 + i],ePurityHelium3TPC[nHe3 + i]),"");
  }
  TLegend *legPIDCutPurHelium3TPCTOF = new TLegend(.01,y1,.54,.88);
  legPIDCutPurHelium3TPCTOF->SetFillStyle(0);
  legPIDCutPurHelium3TPCTOF->SetBorderSize(0);
  legPIDCutPurHelium3TPCTOF->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutPurHelium3TPCTOF->AddEntry((TObject*)0,"Helium3 purity for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutPurHelium3TPCTOF->AddEntry((TObject*)0,Form("#bullet %i:  %.6f #pm %.6f",He3 + i*10,PurityHelium3TPCTOF[nHe3 + i],ePurityHelium3TPCTOF[nHe3 + i]),"");
  }
  TLegend *legPIDCutPurAlphaTPC = new TLegend(.01,y1,.54,.88);
  legPIDCutPurAlphaTPC->SetFillStyle(0);
  legPIDCutPurAlphaTPC->SetBorderSize(0);
  legPIDCutPurAlphaTPC->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutPurAlphaTPC->AddEntry((TObject*)0,"Alpha purity for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutPurAlphaTPC->AddEntry((TObject*)0,Form("#bullet %i:  %.6f #pm %.6f",Alp + i*10,PurityAlphaTPC[nAlp + i],ePurityAlphaTPC[nAlp + i]),"");
  }
  TLegend *legPIDCutPurAlphaTPCTOF = new TLegend(.01,y1,.54,.88);
  legPIDCutPurAlphaTPCTOF->SetFillStyle(0);
  legPIDCutPurAlphaTPCTOF->SetBorderSize(0);
  legPIDCutPurAlphaTPCTOF->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutPurAlphaTPCTOF->AddEntry((TObject*)0,"Alpha purity for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutPurAlphaTPCTOF->AddEntry((TObject*)0,Form("#bullet %i:  %.6f #pm %.6f",Alp + i*10,PurityAlphaTPCTOF[nAlp + i],ePurityAlphaTPCTOF[nAlp + i]),"");
  }
  TLegend *legPIDCutPurHelium3AlphaTPC = new TLegend(.01,y1,.54,.88);
  legPIDCutPurHelium3AlphaTPC->SetFillStyle(0);
  legPIDCutPurHelium3AlphaTPC->SetBorderSize(0);
  legPIDCutPurHelium3AlphaTPC->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutPurHelium3AlphaTPC->AddEntry((TObject*)0,"Z=2 particle purity for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutPurHelium3AlphaTPC->AddEntry((TObject*)0,Form("#bullet %i:  %.6f #pm %.6f",He3Alp + i*10,PurityHelium3AlphaTPC[nHe3Alp + i],ePurityHelium3AlphaTPC[nHe3Alp + i]),"");
  }
  TLegend *legPIDCutPurHelium3AlphaTPCTOF = new TLegend(.01,y1,.54,.88);
  legPIDCutPurHelium3AlphaTPCTOF->SetFillStyle(0);
  legPIDCutPurHelium3AlphaTPCTOF->SetBorderSize(0);
  legPIDCutPurHelium3AlphaTPCTOF->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutPurHelium3AlphaTPCTOF->AddEntry((TObject*)0,"Z=2 particle purity for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutPurHelium3AlphaTPCTOF->AddEntry((TObject*)0,Form("#bullet %i:  %.6f #pm %.6f",He3Alp + i*10,PurityHelium3AlphaTPCTOF[nHe3Alp + i],ePurityHelium3AlphaTPCTOF[nHe3Alp + i]),"");
  }
  TLegend *legPIDCutPurDeuTriHe3AlpTPC = new TLegend(.01,y1,.54,.88);
  legPIDCutPurDeuTriHe3AlpTPC->SetFillStyle(0);
  legPIDCutPurDeuTriHe3AlpTPC->SetBorderSize(0);
  legPIDCutPurDeuTriHe3AlpTPC->AddEntry((TObject*)0,"Possible candidates:","");
  if(AllOrAnti == "All") legPIDCutPurDeuTriHe3AlpTPC->AddEntry((TObject*)0,"d & t & He3 & #alpha purity for <Q_{s}> cuts at:","");
  if(AllOrAnti == "Anti")legPIDCutPurDeuTriHe3AlpTPC->AddEntry((TObject*)0,"#bar{d} & #bar{t} & {}^{3}#bar{He} & #bar{#alpha} purity for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutPurDeuTriHe3AlpTPC->AddEntry((TObject*)0,Form("#bullet %i:  %.6f #pm %.6f",DeuTriHe3Alp + i*10,PurityDeuTriHe3AlpTPC[nDeuTriHe3Alp + i],ePurityDeuTriHe3AlpTPC[nDeuTriHe3Alp + i]),"");
  }
  TLegend *legPIDCutPurDeuTriHe3AlpTPCTOF = new TLegend(.01,y1,.54,.88);
  legPIDCutPurDeuTriHe3AlpTPCTOF->SetFillStyle(0);
  legPIDCutPurDeuTriHe3AlpTPCTOF->SetBorderSize(0);
  legPIDCutPurDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,"Possible candidates:","");
  if(AllOrAnti == "All") legPIDCutPurDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,"d & t & He3 & #alpha purity for <Q_{s}> cuts at:","");
  if(AllOrAnti == "Anti")legPIDCutPurDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,"#bar{d} & #bar{t} & {}^{3}#bar{He} & #bar{#alpha} purity for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutPurDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,Form("#bullet %i:  %.6f #pm %.6f",He3Alp + i*10,PurityDeuTriHe3AlpTPCTOF[nDeuTriHe3Alp + i],ePurityDeuTriHe3AlpTPCTOF[nDeuTriHe3Alp + i]),"");
  }

  TLegend *legnSigmaCutDeuteronTPC = new TLegend(.05,.13,.4,.31);  // same for cd(5,7,8,9)
  legnSigmaCutDeuteronTPC->SetFillStyle(0);
  legnSigmaCutDeuteronTPC->SetBorderSize(0);
  legnSigmaCutDeuteronTPC->AddEntry((TObject*)0,"with preselection:","");
  legnSigmaCutDeuteronTPC->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legnSigmaCutDeuteronTPC->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutDeuteronTPC),"");
  TLegend *legnSigmaCutDeuteronTPCTOF = new TLegend(.05,.13,.4,.36);  // same for cd(5,7,8,9)
  legnSigmaCutDeuteronTPCTOF->SetFillStyle(0);
  legnSigmaCutDeuteronTPCTOF->SetBorderSize(0);
  legnSigmaCutDeuteronTPCTOF->AddEntry((TObject*)0,"with preselection:","");
  legnSigmaCutDeuteronTPCTOF->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legnSigmaCutDeuteronTPCTOF->AddEntry((TObject*)0,Form("|TOFnSigma| = %i#sigma",Sigma),"");
  legnSigmaCutDeuteronTPCTOF->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutDeuteronTPCTOF),"");
  TLegend *legnSigmaCutTritonTPC = new TLegend(.05,.13,.4,.31);  // same for cd(5,7,8,9)
  legnSigmaCutTritonTPC->SetFillStyle(0);
  legnSigmaCutTritonTPC->SetBorderSize(0);
  legnSigmaCutTritonTPC->AddEntry((TObject*)0,"with preselection:","");
  legnSigmaCutTritonTPC->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legnSigmaCutTritonTPC->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutTritonTPC),"");
  TLegend *legnSigmaCutTritonTPCTOF = new TLegend(.05,.13,.4,.36);  // same for cd(5,7,8,9)
  legnSigmaCutTritonTPCTOF->SetFillStyle(0);
  legnSigmaCutTritonTPCTOF->SetBorderSize(0);
  legnSigmaCutTritonTPCTOF->AddEntry((TObject*)0,"with preselection:","");
  legnSigmaCutTritonTPCTOF->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legnSigmaCutTritonTPCTOF->AddEntry((TObject*)0,Form("|TOFnSigma| = %i#sigma",Sigma),"");
  legnSigmaCutTritonTPCTOF->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutTritonTPCTOF),"");
  TLegend *legnSigmaCutHelium3AlphaTPC = new TLegend(.05,.13,.4,.25);  // same for cd(5,7,8,9)
  legnSigmaCutHelium3AlphaTPC->SetFillStyle(0);
  legnSigmaCutHelium3AlphaTPC->SetBorderSize(0);
  legnSigmaCutHelium3AlphaTPC->AddEntry((TObject*)0,"with preselection:","");
  legnSigmaCutHelium3AlphaTPC->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  TLegend *legnSigmaCutHelium3AlphaTPCTOF = new TLegend(.05,.13,.4,.31);  // same for cd(5,7,8,9)
  legnSigmaCutHelium3AlphaTPCTOF->SetFillStyle(0);
  legnSigmaCutHelium3AlphaTPCTOF->SetBorderSize(0);
  legnSigmaCutHelium3AlphaTPCTOF->AddEntry((TObject*)0,"with preselection:","");
  legnSigmaCutHelium3AlphaTPCTOF->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legnSigmaCutHelium3AlphaTPCTOF->AddEntry((TObject*)0,Form("|TOFnSigma| = %i#sigma",Sigma),"");
  TLegend *legnSigmaCutDeuTriHe3AlpTPC = new TLegend(.05,.13,.4,.36);  // same for cd(5,7,8,9)
  legnSigmaCutDeuTriHe3AlpTPC->SetFillStyle(0);
  legnSigmaCutDeuTriHe3AlpTPC->SetBorderSize(0);
  legnSigmaCutDeuTriHe3AlpTPC->AddEntry((TObject*)0,"with preselection:","");
  legnSigmaCutDeuTriHe3AlpTPC->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legnSigmaCutDeuTriHe3AlpTPC->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c} for d",ClearCutDeuteronTPC),"");
  legnSigmaCutDeuTriHe3AlpTPC->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c} for t",ClearCutTritonTPC),"");
  TLegend *legnSigmaCutDeuTriHe3AlpTPC9 = new TLegend(.53,.13,.88,.36);  // same for cd(5,7,8,9)
  legnSigmaCutDeuTriHe3AlpTPC9->SetFillStyle(0);
  legnSigmaCutDeuTriHe3AlpTPC9->SetBorderSize(0);
  legnSigmaCutDeuTriHe3AlpTPC9->AddEntry((TObject*)0,"with preselection:","");
  legnSigmaCutDeuTriHe3AlpTPC9->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legnSigmaCutDeuTriHe3AlpTPC9->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c} for d",ClearCutDeuteronTPC),"");
  legnSigmaCutDeuTriHe3AlpTPC9->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c} for t",ClearCutTritonTPC),"");
  TLegend *legnSigmaCutDeuTriHe3AlpTPCTOF = new TLegend(.05,.13,.4,.42);  // same for cd(5,7,8,9)
  legnSigmaCutDeuTriHe3AlpTPCTOF->SetFillStyle(0);
  legnSigmaCutDeuTriHe3AlpTPCTOF->SetBorderSize(0);
  legnSigmaCutDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,"with preselection:","");
  legnSigmaCutDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legnSigmaCutDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,Form("|TOFnSigma| = %i#sigma",Sigma),"");
  legnSigmaCutDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c} for d",ClearCutDeuteronTPCTOF),"");
  legnSigmaCutDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c} for t",ClearCutTritonTPCTOF),"");
  TLegend *legnSigmaCutDeuTriHe3AlpTPCTOF9 = new TLegend(.53,.13,.88,.42);  // same for cd(5,7,8,9)
  legnSigmaCutDeuTriHe3AlpTPCTOF9->SetFillStyle(0);
  legnSigmaCutDeuTriHe3AlpTPCTOF9->SetBorderSize(0);
  legnSigmaCutDeuTriHe3AlpTPCTOF9->AddEntry((TObject*)0,"with preselection:","");
  legnSigmaCutDeuTriHe3AlpTPCTOF9->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legnSigmaCutDeuTriHe3AlpTPCTOF9->AddEntry((TObject*)0,Form("|TOFnSigma| = %i#sigma",Sigma),"");
  legnSigmaCutDeuTriHe3AlpTPCTOF9->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c} for d",ClearCutDeuteronTPCTOF),"");
  legnSigmaCutDeuTriHe3AlpTPCTOF9->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c} for t",ClearCutTritonTPCTOF),"");

  // legends for cd(6)

  // legends for cd(7)
  TLegend *legPIDCutEffDeuteronTPC = new TLegend(.36,y1,.89,.88);
  legPIDCutEffDeuteronTPC->SetFillStyle(0);
  legPIDCutEffDeuteronTPC->SetBorderSize(0);
  legPIDCutEffDeuteronTPC->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutEffDeuteronTPC->AddEntry((TObject*)0,"Deuteron efficiency for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutEffDeuteronTPC->AddEntry((TObject*)0,Form("#bullet %i:  %.4f #pm %.4f",Deu + i*10,EfficiencyDeuteronTPC[mDeu + i],eEfficiencyDeuteronTPC[mDeu + i]),"");
  }
  TLegend *legPIDCutEffDeuteronTPCTOF = new TLegend(.36,y1,.89,.88);
  legPIDCutEffDeuteronTPCTOF->SetFillStyle(0);
  legPIDCutEffDeuteronTPCTOF->SetBorderSize(0);
  legPIDCutEffDeuteronTPCTOF->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutEffDeuteronTPCTOF->AddEntry((TObject*)0,"Deuteron efficiency for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutEffDeuteronTPCTOF->AddEntry((TObject*)0,Form("#bullet %i:  %.4f #pm %.4f",Deu + i*10,EfficiencyDeuteronTPCTOF[mDeu + i],eEfficiencyDeuteronTPCTOF[mDeu + i]),"");
  }
  TLegend *legPIDCutEffTritonTPC = new TLegend(.36,y1,.89,.88);
  legPIDCutEffTritonTPC->SetFillStyle(0);
  legPIDCutEffTritonTPC->SetBorderSize(0);
  legPIDCutEffTritonTPC->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutEffTritonTPC->AddEntry((TObject*)0,"Triton efficiency for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutEffTritonTPC->AddEntry((TObject*)0,Form("#bullet %i:  %.3f #pm %.3f",Tri + i*10,EfficiencyTritonTPC[mTri + i],eEfficiencyTritonTPC[mTri + i]),"");
  }
  TLegend *legPIDCutEffTritonTPCTOF = new TLegend(.36,y1,.89,.88);
  legPIDCutEffTritonTPCTOF->SetFillStyle(0);
  legPIDCutEffTritonTPCTOF->SetBorderSize(0);
  legPIDCutEffTritonTPCTOF->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutEffTritonTPCTOF->AddEntry((TObject*)0,"Triton efficiency for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutEffTritonTPCTOF->AddEntry((TObject*)0,Form("#bullet %i:  %.3f #pm %.3f",Tri + i*10,EfficiencyTritonTPCTOF[mTri + i],eEfficiencyTritonTPCTOF[mTri + i]),"");
  }
  TLegend *legPIDCutEffHelium3TPC = new TLegend(.36,y1,.89,.88);
  legPIDCutEffHelium3TPC->SetFillStyle(0);
  legPIDCutEffHelium3TPC->SetBorderSize(0);
  legPIDCutEffHelium3TPC->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutEffHelium3TPC->AddEntry((TObject*)0,"Helium3 efficiency for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutEffHelium3TPC->AddEntry((TObject*)0,Form("#bullet %i:  %.3f #pm %.3f",He3 + i*10,EfficiencyHelium3TPC[mHe3 + i],eEfficiencyHelium3TPC[mHe3 + i]),"");
  }
  TLegend *legPIDCutEffHelium3TPCTOF = new TLegend(.36,y1,.89,.88);
  legPIDCutEffHelium3TPCTOF->SetFillStyle(0);
  legPIDCutEffHelium3TPCTOF->SetBorderSize(0);
  legPIDCutEffHelium3TPCTOF->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutEffHelium3TPCTOF->AddEntry((TObject*)0,"Helium3 efficiency for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutEffHelium3TPCTOF->AddEntry((TObject*)0,Form("#bullet %i:  %.3f #pm %.3f",He3 + i*10,EfficiencyHelium3TPCTOF[mHe3 + i],eEfficiencyHelium3TPCTOF[mHe3 + i]),"");
  }
  TLegend *legPIDCutEffAlphaTPC = new TLegend(.36,y1,.89,.88);
  legPIDCutEffAlphaTPC->SetFillStyle(0);
  legPIDCutEffAlphaTPC->SetBorderSize(0);
  legPIDCutEffAlphaTPC->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutEffAlphaTPC->AddEntry((TObject*)0,"Alpha efficiency for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutEffAlphaTPC->AddEntry((TObject*)0,Form("#bullet %i:  %.3f #pm %.3f",Alp + i*10,EfficiencyAlphaTPC[mAlp + i],eEfficiencyAlphaTPC[mAlp + i]),"");
  }
  TLegend *legPIDCutEffAlphaTPCTOF = new TLegend(.36,y1,.89,.88);
  legPIDCutEffAlphaTPCTOF->SetFillStyle(0);
  legPIDCutEffAlphaTPCTOF->SetBorderSize(0);
  legPIDCutEffAlphaTPCTOF->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutEffAlphaTPCTOF->AddEntry((TObject*)0,"Alpha efficiency for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutEffAlphaTPCTOF->AddEntry((TObject*)0,Form("#bullet %i:  %.3f #pm %.3f",Alp + i*10,EfficiencyAlphaTPCTOF[mAlp + i],eEfficiencyAlphaTPCTOF[mAlp + i]),"");
  }
  TLegend *legPIDCutEffHelium3AlphaTPC = new TLegend(.36,y1,.89,.88);
  legPIDCutEffHelium3AlphaTPC->SetFillStyle(0);
  legPIDCutEffHelium3AlphaTPC->SetBorderSize(0);
  legPIDCutEffHelium3AlphaTPC->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutEffHelium3AlphaTPC->AddEntry((TObject*)0,"Z=2 particle efficiency for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutEffHelium3AlphaTPC->AddEntry((TObject*)0,Form("#bullet %i:  %.3f #pm %.3f",He3Alp + i*10,EfficiencyHelium3AlphaTPC[mHe3Alp + i],eEfficiencyHelium3AlphaTPC[mHe3Alp + i]),"");
  }
  TLegend *legPIDCutEffHelium3AlphaTPCTOF = new TLegend(.36,y1,.89,.88);
  legPIDCutEffHelium3AlphaTPCTOF->SetFillStyle(0);
  legPIDCutEffHelium3AlphaTPCTOF->SetBorderSize(0);
  legPIDCutEffHelium3AlphaTPCTOF->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutEffHelium3AlphaTPCTOF->AddEntry((TObject*)0,"Z=2 particle efficiency for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutEffHelium3AlphaTPCTOF->AddEntry((TObject*)0,Form("#bullet %i:  %.3f #pm %.3f",He3Alp + i*10,EfficiencyHelium3AlphaTPCTOF[mHe3Alp + i],eEfficiencyHelium3AlphaTPCTOF[mHe3Alp + i]),"");
  }
  TLegend *legPIDCutEffDeuTriHe3AlpTPC = new TLegend(.36,y1,.89,.88);
  legPIDCutEffDeuTriHe3AlpTPC->SetFillStyle(0);
  legPIDCutEffDeuTriHe3AlpTPC->SetBorderSize(0);
  legPIDCutEffDeuTriHe3AlpTPC->AddEntry((TObject*)0,"Possible candidates:","");
  if(AllOrAnti == "All") legPIDCutEffDeuTriHe3AlpTPC->AddEntry((TObject*)0,"d & t & He3 & #alpha efficiency for <Q_{s}> cuts at:","");
  if(AllOrAnti == "Anti")legPIDCutEffDeuTriHe3AlpTPC->AddEntry((TObject*)0,"#bar{d} & #bar{t} & {}^{3}#bar{He} & #bar{#alpha} efficiency for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutEffDeuTriHe3AlpTPC->AddEntry((TObject*)0,Form("#bullet %i:  %.4f #pm %.4f",DeuTriHe3Alp + i*10,EfficiencyDeuTriHe3AlpTPC[mDeuTriHe3Alp + i],eEfficiencyDeuTriHe3AlpTPC[mDeuTriHe3Alp + i]),"");
  }
  TLegend *legPIDCutEffDeuTriHe3AlpTPCTOF = new TLegend(.36,y1,.89,.88);
  legPIDCutEffDeuTriHe3AlpTPCTOF->SetFillStyle(0);
  legPIDCutEffDeuTriHe3AlpTPCTOF->SetBorderSize(0);
  legPIDCutEffDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,"Possible candidates:","");
  if(AllOrAnti == "All") legPIDCutEffDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,"d & t & He3 & #alpha efficiency for <Q_{s}> cuts at:","");
  if(AllOrAnti == "Anti")legPIDCutEffDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,"#bar{d} & #bar{t} & {}^{3}#bar{He} & #bar{#alpha} efficiency for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutEffDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,Form("#bullet %i:  %.4f #pm %.4f",He3Alp + i*10,EfficiencyDeuTriHe3AlpTPCTOF[mDeuTriHe3Alp + i],eEfficiencyDeuTriHe3AlpTPCTOF[mDeuTriHe3Alp + i]),"");
  }

  // legends for cd(8)
  TLegend *legPIDCutRejDeuteronTPC = new TLegend(.35,y1,.88,.88);
  legPIDCutRejDeuteronTPC->SetFillStyle(0);
  legPIDCutRejDeuteronTPC->SetBorderSize(0);
  legPIDCutRejDeuteronTPC->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutRejDeuteronTPC->AddEntry((TObject*)0,"Deuteron rejection for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutRejDeuteronTPC->AddEntry((TObject*)0,Form("#bullet %i:  %.8f #pm %.8f",Deu + i*10,RejectionDeuteronTPC[mDeu + i],eRejectionDeuteronTPC[mDeu + i]),"");
  }
  TLegend *legPIDCutRejDeuteronTPCTOF = new TLegend(.35,y1,.88,.88);
  legPIDCutRejDeuteronTPCTOF->SetFillStyle(0);
  legPIDCutRejDeuteronTPCTOF->SetBorderSize(0);
  legPIDCutRejDeuteronTPCTOF->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutRejDeuteronTPCTOF->AddEntry((TObject*)0,"Deuteron rejection for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutRejDeuteronTPCTOF->AddEntry((TObject*)0,Form("#bullet %i:  %.8f #pm %.8f",Deu + i*10,RejectionDeuteronTPCTOF[mDeu + i],eRejectionDeuteronTPCTOF[mDeu + i]),"");
  }
  TLegend *legPIDCutRejTritonTPC = new TLegend(.35,y1,.88,.88);
  legPIDCutRejTritonTPC->SetFillStyle(0);
  legPIDCutRejTritonTPC->SetBorderSize(0);
  legPIDCutRejTritonTPC->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutRejTritonTPC->AddEntry((TObject*)0,"Triton rejection for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutRejTritonTPC->AddEntry((TObject*)0,Form("#bullet %i:  %.8f #pm %.8f",Tri + i*10,RejectionTritonTPC[mTri + i],eRejectionTritonTPC[mTri + i]),"");
  }
  TLegend *legPIDCutRejTritonTPCTOF = new TLegend(.35,y1,.88,.88);
  legPIDCutRejTritonTPCTOF->SetFillStyle(0);
  legPIDCutRejTritonTPCTOF->SetBorderSize(0);
  legPIDCutRejTritonTPCTOF->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutRejTritonTPCTOF->AddEntry((TObject*)0,"Triton rejection for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutRejTritonTPCTOF->AddEntry((TObject*)0,Form("#bullet %i:  %.8f #pm %.8f",Tri + i*10,RejectionTritonTPCTOF[mTri + i],eRejectionTritonTPCTOF[mTri + i]),"");
  }
  TLegend *legPIDCutRejHelium3TPC = new TLegend(.35,y1,.88,.88);
  legPIDCutRejHelium3TPC->SetFillStyle(0);
  legPIDCutRejHelium3TPC->SetBorderSize(0);
  legPIDCutRejHelium3TPC->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutRejHelium3TPC->AddEntry((TObject*)0,"Helium3 rejection for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutRejHelium3TPC->AddEntry((TObject*)0,Form("#bullet %i:  %.8f #pm %.8f",He3 + i*10,RejectionHelium3TPC[mHe3 + i],eRejectionHelium3TPC[mHe3 + i]),"");
  }
  TLegend *legPIDCutRejHelium3TPCTOF = new TLegend(.35,y1,.88,.88);
  legPIDCutRejHelium3TPCTOF->SetFillStyle(0);
  legPIDCutRejHelium3TPCTOF->SetBorderSize(0);
  legPIDCutRejHelium3TPCTOF->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutRejHelium3TPCTOF->AddEntry((TObject*)0,"Helium3 rejection for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutRejHelium3TPCTOF->AddEntry((TObject*)0,Form("#bullet %i:  %.8f #pm %.8f",He3 + i*10,RejectionHelium3TPCTOF[mHe3 + i],eRejectionHelium3TPCTOF[mHe3 + i]),"");
  }
  TLegend *legPIDCutRejAlphaTPC = new TLegend(.35,y1,.88,.88);
  legPIDCutRejAlphaTPC->SetFillStyle(0);
  legPIDCutRejAlphaTPC->SetBorderSize(0);
  legPIDCutRejAlphaTPC->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutRejAlphaTPC->AddEntry((TObject*)0,"Alpha rejection for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutRejAlphaTPC->AddEntry((TObject*)0,Form("#bullet %i:  %.8f #pm %.8f",Alp + i*10,RejectionAlphaTPC[mAlp + i],eRejectionAlphaTPC[mAlp + i]),"");
  }
  TLegend *legPIDCutRejAlphaTPCTOF = new TLegend(.35,y1,.88,.88);
  legPIDCutRejAlphaTPCTOF->SetFillStyle(0);
  legPIDCutRejAlphaTPCTOF->SetBorderSize(0);
  legPIDCutRejAlphaTPCTOF->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutRejAlphaTPCTOF->AddEntry((TObject*)0,"Alpha rejection for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutRejAlphaTPCTOF->AddEntry((TObject*)0,Form("#bullet %i:  %.8f #pm %.8f",Alp + i*10,RejectionAlphaTPCTOF[mAlp + i],eRejectionAlphaTPCTOF[mAlp + i]),"");
  }
  TLegend *legPIDCutRejHelium3AlphaTPC = new TLegend(.35,y1,.88,.88);
  legPIDCutRejHelium3AlphaTPC->SetFillStyle(0);
  legPIDCutRejHelium3AlphaTPC->SetBorderSize(0);
  legPIDCutRejHelium3AlphaTPC->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutRejHelium3AlphaTPC->AddEntry((TObject*)0,"Z=2 particle rejection for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutRejHelium3AlphaTPC->AddEntry((TObject*)0,Form("#bullet %i:  %.8f #pm %.8f",He3Alp + i*10,RejectionHelium3AlphaTPC[mHe3Alp + i],eRejectionHelium3AlphaTPC[mHe3Alp + i]),"");
  }
  TLegend *legPIDCutRejHelium3AlphaTPCTOF = new TLegend(.35,y1,.88,.88);
  legPIDCutRejHelium3AlphaTPCTOF->SetFillStyle(0);
  legPIDCutRejHelium3AlphaTPCTOF->SetBorderSize(0);
  legPIDCutRejHelium3AlphaTPCTOF->AddEntry((TObject*)0,"Possible candidates:","");
  legPIDCutRejHelium3AlphaTPCTOF->AddEntry((TObject*)0,"Z=2 particle rejection for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutRejHelium3AlphaTPCTOF->AddEntry((TObject*)0,Form("#bullet %i:  %.8f #pm %.8f",He3Alp + i*10,RejectionHelium3AlphaTPCTOF[mHe3Alp + i],eRejectionHelium3AlphaTPCTOF[mHe3Alp + i]),"");
  }
  TLegend *legPIDCutRejDeuTriHe3AlpTPC = new TLegend(.35,y1,.88,.88);
  legPIDCutRejDeuTriHe3AlpTPC->SetFillStyle(0);
  legPIDCutRejDeuTriHe3AlpTPC->SetBorderSize(0);
  legPIDCutRejDeuTriHe3AlpTPC->AddEntry((TObject*)0,"Possible candidates:","");
  if(AllOrAnti == "All") legPIDCutRejDeuTriHe3AlpTPC->AddEntry((TObject*)0,"d & t & He3 & #alpha rejection for <Q_{s}> cuts at:","");
  if(AllOrAnti == "Anti")legPIDCutRejDeuTriHe3AlpTPC->AddEntry((TObject*)0,"#bar{d} & #bar{t} & {}^{3}#bar{He} & #bar{#alpha} rejection for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutRejDeuTriHe3AlpTPC->AddEntry((TObject*)0,Form("#bullet %i:  %.8f #pm %.8f",DeuTriHe3Alp + i*10,RejectionDeuTriHe3AlpTPC[mDeuTriHe3Alp + i],eRejectionDeuTriHe3AlpTPC[mDeuTriHe3Alp + i]),"");
  }
  TLegend *legPIDCutRejDeuTriHe3AlpTPCTOF = new TLegend(.35,y1,.88,.88);
  legPIDCutRejDeuTriHe3AlpTPCTOF->SetFillStyle(0);
  legPIDCutRejDeuTriHe3AlpTPCTOF->SetBorderSize(0);
  legPIDCutRejDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,"Possible candidates:","");
  if(AllOrAnti == "All") legPIDCutRejDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,"d & t & He3 & #alpha rejection for <Q_{s}> cuts at:","");
  if(AllOrAnti == "Anti")legPIDCutRejDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,"#bar{d} & #bar{t} & {}^{3}#bar{He} & #bar{#alpha} rejection for <Q_{s}> cuts at:","");
  for(Int_t i = 0; i < nopc; i++) {
    legPIDCutRejDeuTriHe3AlpTPCTOF->AddEntry((TObject*)0,Form("#bullet %i:  %.8f #pm %.8f",He3Alp + i*10,RejectionDeuTriHe3AlpTPCTOF[mDeuTriHe3Alp + i],eRejectionDeuTriHe3AlpTPCTOF[mDeuTriHe3Alp + i]),"");
  }

  // legends for cd(9)
  TLegend *legPIDCut180Helium3AlphaTPC = new TLegend(.27,.47,.5,.54);
  legPIDCut180Helium3AlphaTPC->SetFillStyle(0);
  legPIDCut180Helium3AlphaTPC->SetBorderSize(0);
  legPIDCut180Helium3AlphaTPC->AddEntry((TObject*)0,"<Q_{s}> Cut: 180","");
  TLegend *legPIDCut170Helium3AlphaTPC = new TLegend(.5,.57,.62,.65);
  legPIDCut170Helium3AlphaTPC->SetFillStyle(0);
  legPIDCut170Helium3AlphaTPC->SetBorderSize(0);
  legPIDCut170Helium3AlphaTPC->AddEntry((TObject*)0,"170","");
  TLegend *legPIDCut160Helium3AlphaTPC = new TLegend(.56,.6,.68,.68);
  legPIDCut160Helium3AlphaTPC->SetFillStyle(0);
  legPIDCut160Helium3AlphaTPC->SetBorderSize(0);
  legPIDCut160Helium3AlphaTPC->AddEntry((TObject*)0,"160","");
  TLegend *legPIDCut180Helium3AlphaTPCTOF = new TLegend(.27,.48,.5,.55);
  legPIDCut180Helium3AlphaTPCTOF->SetFillStyle(0);
  legPIDCut180Helium3AlphaTPCTOF->SetBorderSize(0);
  legPIDCut180Helium3AlphaTPCTOF->AddEntry((TObject*)0,"<Q_{s}> Cut: 180","");
  TLegend *legPIDCut170Helium3AlphaTPCTOF = new TLegend(.48,.57,.6,.65);
  legPIDCut170Helium3AlphaTPCTOF->SetFillStyle(0);
  legPIDCut170Helium3AlphaTPCTOF->SetBorderSize(0);
  legPIDCut170Helium3AlphaTPCTOF->AddEntry((TObject*)0,"170","");
  TLegend *legPIDCut160Helium3AlphaTPCTOF = new TLegend(.53,.61,.65,.69);
  legPIDCut160Helium3AlphaTPCTOF->SetFillStyle(0);
  legPIDCut160Helium3AlphaTPCTOF->SetBorderSize(0);
  legPIDCut160Helium3AlphaTPCTOF->AddEntry((TObject*)0,"160","");

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Create and fill canvases
  ///_________________________________________________________________________________________________
  TCanvas *c1 = new TCanvas("c1",Form("Deuteron Efficiency, Rejection, Purity for several cuts (TPC%iSigma ps)",Sigma));
  c1->Divide(3,3);
  c1->cd(1)->SetLogz();
    fHistTRDpTvPID->Draw("colz");
    legPID->Draw("same");
  c1->cd(2)->SetLogz();
    fHistTRDpTvPIDvDeuteronTPC_clear->Draw("colz");
    legPIDDeuteronTPC->Draw("same");
  c1->cd(3);
    fHistTRDPIDDeuteronTPC_clear_Clone->Draw();
    legProjDeuteronTPC->Draw("same");
  c1->cd(4);
    fHistTRDPIDallDeuteronTPC->Draw();
    fHistTRDPIDAllexceptDeuteronTPC_clear->Draw("same");
    fHistTRDPIDDeuteronTPC_clear_Clone2->Draw("same");
    legMeanDeuteronTPC->Draw("same");
    legMeanallDeuteronTPC->Draw("same");
    legParticleDeuteron->Draw("same");
  c1->cd(5);
    fPurityDeuteronTPC->Draw("AP");
    legPIDCutPurDeuteronTPC->Draw("same");
    legnSigmaCutDeuteronTPC->Draw("same");
  c1->cd(6);

  c1->cd(7);
    fEfficiencyDeuteronTPC->Draw("AP");
    legPIDCutEffDeuteronTPC->Draw("same");
    legnSigmaCutDeuteronTPC->Draw("same");
  c1->cd(8)->SetLogy();
    fRejectionDeuteronTPC->Draw("AP");
    legPIDCutRejDeuteronTPC->Draw("same");
    legnSigmaCutDeuteronTPC->Draw("same");
  c1->cd(9)->SetLogy();
    fEffRejDeuteronTPC->Draw("AP");
    legnSigmaCutDeuteronTPC->Draw("same");

  TCanvas *c2 = new TCanvas("c2",Form("Deuteron Efficiency, Rejection, Purity for several cuts ((TPC&TOF)%iSigma ps)",Sigma));
  c2->Divide(3,3);
  c2->cd(1)->SetLogz();
    fHistTRDpTvPID->Draw("colz");
    legPID->Draw("same");
  c2->cd(2)->SetLogz();
    fHistTRDpTvPIDvDeuteronTPCTOF_clear->Draw("colz");
    legPIDDeuteronTPCTOF->Draw("same");
  c2->cd(3);
    fHistTRDPIDDeuteronTPCTOF_clear_Clone->Draw();
    legProjDeuteronTPCTOF->Draw("same");
  c2->cd(4);
    fHistTRDPIDallDeuteronTPCTOF->Draw();
    fHistTRDPIDAllexceptDeuteronTPCTOF_clear->Draw("same");
    fHistTRDPIDDeuteronTPCTOF_clear_Clone2->Draw("same");
    legMeanDeuteronTPCTOF->Draw("same");
    legMeanallDeuteronTPCTOF->Draw("same");
    legParticleDeuteron->Draw("same");
  c2->cd(5);
    fPurityDeuteronTPCTOF->Draw("AP");
    legPIDCutPurDeuteronTPCTOF->Draw("same");
    legnSigmaCutDeuteronTPCTOF->Draw("same");
  c2->cd(6);

  c2->cd(7);
    fEfficiencyDeuteronTPCTOF->Draw("AP");
    legPIDCutEffDeuteronTPCTOF->Draw("same");
    legnSigmaCutDeuteronTPCTOF->Draw("same");
  c2->cd(8)->SetLogy();
    fRejectionDeuteronTPCTOF->Draw("AP");
    legPIDCutRejDeuteronTPCTOF->Draw("same");
    legnSigmaCutDeuteronTPCTOF->Draw("same");
  c2->cd(9)->SetLogy();
    fEffRejDeuteronTPCTOF->Draw("AP");
    legnSigmaCutDeuteronTPCTOF->Draw("same");

  TCanvas *c3 = new TCanvas("c3",Form("Triton Efficiency, Rejection, Purity for several cuts (TPC%iSigma ps)",Sigma));
  c3->Divide(3,3);
  c3->cd(1)->SetLogz();
    fHistTRDpTvPID->Draw("colz");
    legPID->Draw("same");
  c3->cd(2)->SetLogz();
    fHistTRDpTvPIDvTritonTPC_clear->Draw("colz");
    legPIDTritonTPC->Draw("same");
  c3->cd(3);
    fHistTRDPIDTritonTPC_clear_Clone->Draw();
    legProjTritonTPC->Draw("same");
  c3->cd(4);
    fHistTRDPIDallTritonTPC->Draw();
    fHistTRDPIDAllexceptTritonTPC_clear->Draw("same");
    fHistTRDPIDTritonTPC_clear_Clone2->Draw("same");
    legMeanTritonTPC->Draw("same");
    legMeanallTritonTPC->Draw("same");
    legParticleTriton->Draw("same");
  c3->cd(5);
    fPurityTritonTPC->Draw("AP");
    legPIDCutPurTritonTPC->Draw("same");
    legnSigmaCutTritonTPC->Draw("same");
  c3->cd(6);

  c3->cd(7);
    fEfficiencyTritonTPC->Draw("AP");
    legPIDCutEffTritonTPC->Draw("same");
    legnSigmaCutTritonTPC->Draw("same");
  c3->cd(8)->SetLogy();
    fRejectionTritonTPC->Draw("AP");
    legPIDCutRejTritonTPC->Draw("same");
    legnSigmaCutTritonTPC->Draw("same");
  c3->cd(9)->SetLogy();
    fEffRejTritonTPC->Draw("AP");
    legnSigmaCutTritonTPC->Draw("same");

  TCanvas *c4 = new TCanvas("c4",Form("Triton Efficiency, Rejection, Purity for several cuts ((TPC&TOF)%iSigma ps)",Sigma));
  c4->Divide(3,3);
  c4->cd(1)->SetLogz();
    fHistTRDpTvPID->Draw("colz");
    legPID->Draw("same");
  c4->cd(2)->SetLogz();
    fHistTRDpTvPIDvTritonTPCTOF_clear->Draw("colz");
    legPIDTritonTPCTOF->Draw("same");
  c4->cd(3);
    fHistTRDPIDTritonTPCTOF_clear_Clone->Draw();
    legProjTritonTPCTOF->Draw("same");
  c4->cd(4);
    fHistTRDPIDallTritonTPCTOF->Draw();
    fHistTRDPIDAllexceptTritonTPCTOF_clear->Draw("same");
    fHistTRDPIDTritonTPCTOF_clear_Clone2->Draw("same");
    legMeanTritonTPCTOF->Draw("same");
    legMeanallTritonTPCTOF->Draw("same");
    legParticleTriton->Draw("same");
  c4->cd(5);
    fPurityTritonTPCTOF->Draw("AP");
    legPIDCutPurTritonTPCTOF->Draw("same");
    legnSigmaCutTritonTPCTOF->Draw("same");
  c4->cd(6);

  c4->cd(7);
    fEfficiencyTritonTPCTOF->Draw("AP");
    legPIDCutEffTritonTPCTOF->Draw("same");
    legnSigmaCutTritonTPCTOF->Draw("same");
  c4->cd(8)->SetLogy();
    fRejectionTritonTPCTOF->Draw("AP");
    legPIDCutRejTritonTPCTOF->Draw("same");
    legnSigmaCutTritonTPCTOF->Draw("same");
  c4->cd(9)->SetLogy();
    fEffRejTritonTPCTOF->Draw("AP");
    legnSigmaCutTritonTPCTOF->Draw("same");

  TCanvas *c5 = new TCanvas("c5",Form("Helium3 Efficiency, Rejection, Purity for several cuts (TPC%iSigma ps)",Sigma));
  c5->Divide(3,3);
  c5->cd(1)->SetLogz();
    fHistTRDpTvPID->Draw("colz");
    legPID->Draw("same");
  c5->cd(2)->SetLogz();
    fHistTRDpTvPIDvHelium3TPC->Draw("colz");
    legPIDTPC->Draw("same");
  c5->cd(3);
    fHistTRDPIDHelium3TPC_Clone->Draw();
    legProjHelium3TPC->Draw("same");
  c5->cd(4);
    fHistTRDPIDallHelium3TPC->Draw();
    fHistTRDPIDAllexceptHelium3TPC->Draw("same");
    fHistTRDPIDHelium3TPC_Clone2->Draw("same");
    legMeanHelium3TPC->Draw("same");
    legMeanallHelium3TPC->Draw("same");
    legParticleHelium3->Draw("same");
  c5->cd(5);
    fPurityHelium3TPC->Draw("AP");
    legPIDCutPurHelium3TPC->Draw("same");
    legnSigmaCutHelium3AlphaTPC->Draw("same");
  c5->cd(6)->SetLogy();
    fPurityHelium3TPC->Draw("AP");
  c5->cd(7);
    fEfficiencyHelium3TPC->Draw("AP");
    legPIDCutEffHelium3TPC->Draw("same");
    legnSigmaCutHelium3AlphaTPC->Draw("same");
  c5->cd(8)->SetLogy();
    fRejectionHelium3TPC->Draw("AP");
    legPIDCutRejHelium3TPC->Draw("same");
    legnSigmaCutHelium3AlphaTPC->Draw("same");
  c5->cd(9)->SetLogy();
    fEffRejHelium3TPC->Draw("AP");
    legnSigmaCutHelium3AlphaTPC->Draw("same");

  TCanvas *c6 = new TCanvas("c6",Form("Helium3 Efficiency, Rejection, Purity for several cuts ((TPC&TOF)%iSigma ps)",Sigma));
  c6->Divide(3,3);
  c6->cd(1)->SetLogz();
    fHistTRDpTvPID->Draw("colz");
    legPID->Draw("same");
  c6->cd(2)->SetLogz();
    fHistTRDpTvPIDvHelium3TPCTOF->Draw("colz");
    legPIDTPCTOF->Draw("same");
  c6->cd(3);
    fHistTRDPIDHelium3TPCTOF_Clone->Draw();
    legProjHelium3TPCTOF->Draw("same");
  c6->cd(4);
    fHistTRDPIDallHelium3TPCTOF->Draw();
    fHistTRDPIDAllexceptHelium3TPCTOF->Draw("same");
    fHistTRDPIDHelium3TPCTOF_Clone2->Draw("same");
    legMeanHelium3TPCTOF->Draw("same");
    legMeanallHelium3TPCTOF->Draw("same");
    legParticleHelium3->Draw("same");
  c6->cd(5);
    fPurityHelium3TPCTOF->Draw("AP");
    legPIDCutPurHelium3TPCTOF->Draw("same");
    legnSigmaCutHelium3AlphaTPCTOF->Draw("same");
  c6->cd(6)->SetLogy();
    fPurityHelium3TPCTOF->Draw("AP");
  c6->cd(7);
    fEfficiencyHelium3TPCTOF->Draw("AP");
    legPIDCutEffHelium3TPCTOF->Draw("same");
    legnSigmaCutHelium3AlphaTPCTOF->Draw("same");
  c6->cd(8)->SetLogy();
    fRejectionHelium3TPCTOF->Draw("AP");
    legPIDCutRejHelium3TPCTOF->Draw("same");
    legnSigmaCutHelium3AlphaTPCTOF->Draw("same");
  c6->cd(9)->SetLogy();
    fEffRejHelium3TPCTOF->Draw("AP");
    legnSigmaCutHelium3AlphaTPCTOF->Draw("same");

  TCanvas *c7 = new TCanvas("c7",Form("Alpha Efficiency, Rejection, Purity for several cuts (TPC%iSigma ps)",Sigma));
  c7->Divide(3,3);
  c7->cd(1)->SetLogz();
    fHistTRDpTvPID->Draw("colz");
    legPID->Draw("same");
  c7->cd(2)->SetLogz();
    fHistTRDpTvPIDvAlphaTPC->Draw("colz");
    legPIDTPC->Draw("same");
  c7->cd(3);
    fHistTRDPIDAlphaTPC_Clone->Draw();
    legProjAlphaTPC->Draw("same");
  c7->cd(4);
    fHistTRDPIDallAlphaTPC->Draw();
    fHistTRDPIDAllexceptAlphaTPC->Draw("same");
    fHistTRDPIDAlphaTPC_Clone2->Draw("same");
    legMeanAlphaTPC->Draw("same");
    legMeanallAlphaTPC->Draw("same");
    legParticleAlpha->Draw("same");
  c7->cd(5);
    fPurityAlphaTPC->Draw("AP");
    legPIDCutPurAlphaTPC->Draw("same");
    legnSigmaCutHelium3AlphaTPC->Draw("same");
  c7->cd(6)->SetLogy();
    fPurityAlphaTPC->Draw("AP");
  c7->cd(7);
    fEfficiencyAlphaTPC->Draw("AP");
    legPIDCutEffAlphaTPC->Draw("same");
    legnSigmaCutHelium3AlphaTPC->Draw("same");
  c7->cd(8)->SetLogy();
    fRejectionAlphaTPC->Draw("AP");
    legPIDCutRejAlphaTPC->Draw("same");
    legnSigmaCutHelium3AlphaTPC->Draw("same");
  c7->cd(9)->SetLogy();
    fEffRejAlphaTPC->Draw("AP");
    legnSigmaCutHelium3AlphaTPC->Draw("same");

  TCanvas *c8 = new TCanvas("c8",Form("Alpha Efficiency, Rejection, Purity for several cuts ((TPC&TOF)%iSigma ps)",Sigma));
  c8->Divide(3,3);
  c8->cd(1)->SetLogz();
    fHistTRDpTvPID->Draw("colz");
    legPID->Draw("same");
  c8->cd(2)->SetLogz();
    fHistTRDpTvPIDvAlphaTPCTOF->Draw("colz");
    legPIDTPCTOF->Draw("same");
  c8->cd(3);
    fHistTRDPIDAlphaTPCTOF_Clone->Draw();
    legProjAlphaTPCTOF->Draw("same");
  c8->cd(4);
    fHistTRDPIDallAlphaTPCTOF->Draw();
    fHistTRDPIDAllexceptAlphaTPCTOF->Draw("same");
    fHistTRDPIDAlphaTPCTOF_Clone2->Draw("same");
    legMeanAlphaTPCTOF->Draw("same");
    legMeanallAlphaTPCTOF->Draw("same");
    legParticleAlpha->Draw("same");
  c8->cd(5);
    fPurityAlphaTPCTOF->Draw("AP");
    legPIDCutPurAlphaTPCTOF->Draw("same");
    legnSigmaCutHelium3AlphaTPCTOF->Draw("same");
  c8->cd(6)->SetLogy();
    fPurityAlphaTPCTOF->Draw("AP");
  c8->cd(7);
    fEfficiencyAlphaTPCTOF->Draw("AP");
    legPIDCutEffAlphaTPCTOF->Draw("same");
    legnSigmaCutHelium3AlphaTPCTOF->Draw("same");
  c8->cd(8)->SetLogy();
    fRejectionAlphaTPCTOF->Draw("AP");
    legPIDCutRejAlphaTPCTOF->Draw("same");
    legnSigmaCutHelium3AlphaTPCTOF->Draw("same");
  c8->cd(9)->SetLogy();
    fEffRejAlphaTPCTOF->Draw("AP");
    legnSigmaCutHelium3AlphaTPCTOF->Draw("same");

  TCanvas *c9 = new TCanvas("c9",Form("Z=2 particle Efficiency, Rejection, Purity for several cuts (TPC%iSigma ps)",Sigma));
  c9->Divide(3,3);
  c9->cd(1)->SetLogz();
    fHistTRDpTvPID->Draw("colz");
    legPID->Draw("same");
  c9->cd(2)->SetLogz();
    fHistTRDpTvPIDvHelium3AlphaTPC->Draw("colz");
    legPIDTPC->Draw("same");
  c9->cd(3);
    fHistTRDPIDHelium3AlphaTPC_Clone->Draw();
    legProjHelium3AlphaTPC->Draw("same");
  c9->cd(4);
    fHistTRDPIDallHelium3AlphaTPC->Draw();
    fHistTRDPIDAllexceptHelium3AlphaTPC->Draw("same");
    fHistTRDPIDHelium3AlphaTPC->Draw("same");
    legMeanHelium3AlphaTPC->Draw("same");
    legMeanallHelium3AlphaTPC->Draw("same");
    legParticleHelium3Alpha->Draw("same");
  c9->cd(5);
    fPurityHelium3AlphaTPC->Draw("AP");
    legPIDCutPurHelium3AlphaTPC->Draw("same");
    legnSigmaCutHelium3AlphaTPC->Draw("same");
  c9->cd(6)->SetLogy();
    fPurityHelium3AlphaTPC->Draw("AP");
  c9->cd(7);
    fEfficiencyHelium3AlphaTPC->Draw("AP");
    legPIDCutEffHelium3AlphaTPC->Draw("same");
    legnSigmaCutHelium3AlphaTPC->Draw("same");
  c9->cd(8)->SetLogy();
    fRejectionHelium3AlphaTPC->Draw("AP");
    legPIDCutRejHelium3AlphaTPC->Draw("same");
    legnSigmaCutHelium3AlphaTPC->Draw("same");
  c9->cd(9)->SetLogy();
    fEffRejHelium3AlphaTPC->Draw("AP");
    //legPIDCut180Helium3AlphaTPC->Draw("same");
    //legPIDCut170Helium3AlphaTPC->Draw("same");
    //legPIDCut160Helium3AlphaTPC->Draw("same");
    legnSigmaCutHelium3AlphaTPC->Draw("same");

  TCanvas *c10 = new TCanvas("c10",Form("Z=2 particle Efficiency, Rejection, Purity for several cuts ((TPC&TOF)%iSigma ps)",Sigma));
  c10->Divide(3,3);
  c10->cd(1)->SetLogz();
    fHistTRDpTvPID->Draw("colz");
    legPID->Draw("same");
  c10->cd(2)->SetLogz();
    fHistTRDpTvPIDvHelium3AlphaTPCTOF->Draw("colz");
    legPIDTPCTOF->Draw("same");
  c10->cd(3);
    fHistTRDPIDHelium3AlphaTPCTOF_Clone->Draw();
    legProjHelium3AlphaTPCTOF->Draw("same");
  c10->cd(4);
    fHistTRDPIDallHelium3AlphaTPCTOF->Draw();
    fHistTRDPIDAllexceptHelium3AlphaTPCTOF->Draw("same");
    fHistTRDPIDHelium3AlphaTPCTOF->Draw("same");
    legMeanHelium3AlphaTPCTOF->Draw("same");
    legMeanallHelium3AlphaTPCTOF->Draw("same");
    legParticleHelium3Alpha->Draw("same");
  c10->cd(5);
    fPurityHelium3AlphaTPCTOF->Draw("AP");
    legPIDCutPurHelium3AlphaTPCTOF->Draw("same");
    legnSigmaCutHelium3AlphaTPCTOF->Draw("same");
  c10->cd(6)->SetLogy();
    fPurityHelium3AlphaTPCTOF->Draw("AP");
  c10->cd(7);
    fEfficiencyHelium3AlphaTPCTOF->Draw("AP");
    legPIDCutEffHelium3AlphaTPCTOF->Draw("same");
    legnSigmaCutHelium3AlphaTPCTOF->Draw("same");
  c10->cd(8)->SetLogy();
    fRejectionHelium3AlphaTPCTOF->Draw("AP");
    legPIDCutRejHelium3AlphaTPCTOF->Draw("same");
    legnSigmaCutHelium3AlphaTPCTOF->Draw("same");
  c10->cd(9)->SetLogy();
    fEffRejHelium3AlphaTPCTOF->Draw("AP");
    //legPIDCut180Helium3AlphaTPCTOF->Draw("same");
    //legPIDCut170Helium3AlphaTPCTOF->Draw("same");
    //legPIDCut160Helium3AlphaTPCTOF->Draw("same");
    legnSigmaCutHelium3AlphaTPCTOF->Draw("same");

  TCanvas *c11 = new TCanvas("c11",Form("d & t & He3 & Alpha Efficiency, Rejection, Purity for several cuts (TPC%iSigma ps)",Sigma));
  c11->Divide(3,3);
  c11->cd(1)->SetLogz();
    fHistTRDpTvPID->Draw("colz");
    legPID->Draw("same");
  c11->cd(2)->SetLogz();
    fHistTRDpTvPIDvDeuTriHe3AlpTPC->Draw("colz");
    legPIDDeuTriHe3AlpTPC->Draw("same");
  c11->cd(3);
    fHistTRDPIDDeuTriHe3AlpTPC_Clone->Draw();
    legProjDeuTriHe3AlpTPC->Draw("same");
  c11->cd(4);
    fHistTRDPIDallDeuTriHe3AlpTPC->Draw();
    fHistTRDPIDAllexceptDeuTriHe3AlpTPC->Draw("same");
    fHistTRDPIDDeuTriHe3AlpTPC->Draw("same");
    legMeanDeuTriHe3AlpTPC->Draw("same");
    legMeanallDeuTriHe3AlpTPC->Draw("same");
    legParticleDeuTriHe3Alp->Draw("same");
  c11->cd(5);
    fPurityDeuTriHe3AlpTPC->Draw("AP");
    legPIDCutPurDeuTriHe3AlpTPC->Draw("same");
    legnSigmaCutDeuTriHe3AlpTPC->Draw("same");
  c11->cd(6)->SetLogy();
    fPurityDeuTriHe3AlpTPC->Draw("AP");
  c11->cd(7);
    fEfficiencyDeuTriHe3AlpTPC->Draw("AP");
    legPIDCutEffDeuTriHe3AlpTPC->Draw("same");
    legnSigmaCutDeuTriHe3AlpTPC->Draw("same");
  c11->cd(8)->SetLogy();
    fRejectionDeuTriHe3AlpTPC->Draw("AP");
    legPIDCutRejDeuTriHe3AlpTPC->Draw("same");
    legnSigmaCutDeuTriHe3AlpTPC->Draw("same");
  c11->cd(9)->SetLogy();
    fEffRejDeuTriHe3AlpTPC->Draw("AP");
    legnSigmaCutDeuTriHe3AlpTPC9->Draw("same");

  TCanvas *c12 = new TCanvas("c12",Form("d & t & He3 & Alpha Efficiency, Rejection, Purity for several cuts ((TPC&TOF)%iSigma ps)",Sigma));
  c12->Divide(3,3);
  c12->cd(1)->SetLogz();
    fHistTRDpTvPID->Draw("colz");
    legPID->Draw("same");
  c12->cd(2)->SetLogz();
    fHistTRDpTvPIDvDeuTriHe3AlpTPCTOF->Draw("colz");
    legPIDDeuTriHe3AlpTPCTOF->Draw("same");
  c12->cd(3);
    fHistTRDPIDDeuTriHe3AlpTPCTOF_Clone->Draw();
    legProjDeuTriHe3AlpTPCTOF->Draw("same");
  c12->cd(4);
    fHistTRDPIDallDeuTriHe3AlpTPCTOF->Draw();
    fHistTRDPIDAllexceptDeuTriHe3AlpTPCTOF->Draw("same");
    fHistTRDPIDDeuTriHe3AlpTPCTOF->Draw("same");
    legMeanDeuTriHe3AlpTPCTOF->Draw("same");
    legMeanallDeuTriHe3AlpTPCTOF->Draw("same");
    legParticleDeuTriHe3Alp->Draw("same");
  c12->cd(5);
    fPurityDeuTriHe3AlpTPCTOF->Draw("AP");
    legPIDCutPurDeuTriHe3AlpTPCTOF->Draw("same");
    legnSigmaCutDeuTriHe3AlpTPCTOF->Draw("same");
  c12->cd(6)->SetLogy();
    fPurityDeuTriHe3AlpTPCTOF->Draw("AP");
  c12->cd(7);
    fEfficiencyDeuTriHe3AlpTPCTOF->Draw("AP");
    legPIDCutEffDeuTriHe3AlpTPCTOF->Draw("same");
    legnSigmaCutDeuTriHe3AlpTPCTOF->Draw("same");
  c12->cd(8)->SetLogy();
    fRejectionDeuTriHe3AlpTPCTOF->Draw("AP");
    legPIDCutRejDeuTriHe3AlpTPCTOF->Draw("same");
    legnSigmaCutDeuTriHe3AlpTPCTOF->Draw("same");
  c12->cd(9)->SetLogy();
    fEffRejDeuTriHe3AlpTPCTOF->Draw("AP");
    legnSigmaCutDeuTriHe3AlpTPCTOF9->Draw("same");


  //TCanvas *c13 = new TCanvas ("c13","d & t & He3 & Alpha");
  //c13->Divide(1,1);
  //c13->cd(1);
  //  fHistTRDPIDallDeuTriHe3AlpTPC->Draw();
  //  fHistTRDPIDAllexceptDeuTriHe3AlpTPC->Draw("same");
  //  fHistTRDPIDDeuTriHe3AlpTPC->Draw("same");
  //  legMeanDeuTriHe3AlpTPC->Draw("same");
  //  legMeanallDeuTriHe3AlpTPC->Draw("same");
  //  legParticleDeuTriHe3Alp->Draw("same");

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Save canvases
  ///____________________________________________________________________________________________
  if(AllOrAnti == "All") {
    c1 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_DeuteronTPC%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c1 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_DeuteronTPC%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
    c2 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_DeuteronTPCTOF%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c2 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_DeuteronTPCTOF%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
    c3 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_TritonTPC%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c3 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_TritonTPC%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
    c4 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_TritonTPCTOF%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c4 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_TritonTPCTOF%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
    c5 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_Helium3TPC%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c5 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_Helium3TPC%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
    c6 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_Helium3TPCTOF%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c6 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_Helium3TPCTOF%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
    c7 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AlphaTPC%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c7 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AlphaTPC%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
    c8 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AlphaTPCTOF%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c8 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AlphaTPCTOF%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
    c9 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_Z2particlesTPC%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c9 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_Z2particlesTPC%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
    c10->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_Z2particlesTPCTOF%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c10->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_Z2particlesTPCTOF%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
    c11->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_DeuTriHe3AlpTPC%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c11->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_DeuTriHe3AlpTPC%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
    c12->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_DeuTriHe3AlpTPCTOF%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c12->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_DeuTriHe3AlpTPCTOF%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
  }
  if(AllOrAnti == "Anti") {
    c1 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiDeuteronTPC%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c1 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiDeuteronTPC%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
    c2 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiDeuteronTPCTOF%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c2 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiDeuteronTPCTOF%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
    c3 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiTritonTPC%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c3 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiTritonTPC%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
    c4 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiTritonTPCTOF%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c4 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiTritonTPCTOF%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
    c5 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiHelium3TPC%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c5 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiHelium3TPC%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
    c6 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiHelium3TPCTOF%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c6 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiHelium3TPCTOF%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
    c7 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiAlphaTPC%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c7 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiAlphaTPC%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
    c8 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiAlphaTPCTOF%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c8 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiAlphaTPCTOF%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
    c9 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiZ2particlesTPC%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c9 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiZ2particlesTPC%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
    c10->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiZ2particlesTPCTOF%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c10->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiZ2particlesTPCTOF%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
    c11->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiDeuTriHe3AlpTPC%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c11->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiDeuTriHe3AlpTPC%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
    c12->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiDeuTriHe3AlpTPCTOF%iSigma%s.pdf",Collisions,Sigma,Sigma,Collisions));
    c12->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/bbrudnyj_ReadTreeTestTRDTriggerERP/plots%iSigma/bbrudnyj_AntiDeuTriHe3AlpTPCTOF%iSigma%s.root",Collisions,Sigma,Sigma,Collisions));
  }

  printf("\n\ndone!\n");


}


Double_t GaussError(Int_t a, Double_t da, Int_t b, Double_t db) // for function f = a/b
{

  return TMath::Sqrt( (da/b)*(da/b) + (((a*db)/b)/b)*(((a*db)/b)/b) );

}
