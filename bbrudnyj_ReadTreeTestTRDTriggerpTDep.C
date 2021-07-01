/*--------------------------------------------------------------------------------------------------*\
|                                                                                                    |
|    Macro that makes projection histograms with branches of "nmartin_TestTriggerTRDTree.root"       |
|                                                                                                    |
|    - pT Dependence for  Efficiency, Rejection and Purity for several cuts in <Q_s> for             |
|      Deuteron, Triton, Helium3, Alpha, Helium3 & Alpha, d & t & He3 & Alpha and their Antis        |
|                                                                                                    |
| Author: Benjamin Brudnyj (2016)                                                                    | 
\*--------------------------------------------------------------------------------------------------*/
void bbrudnyj_ReadTreeTestTRDTriggerpTDep() {

  static const Int_t n              = 26;  // for purity
  static const Int_t m              = 26;  // for efficiency and rejection
  static const Int_t nParticles     = 6;
  static const Int_t nPreselections = 2;
  static const Int_t nAnalysis      = 3;
  static const Int_t nFiles         = 2;
  static const Int_t nPads          = 9;
  static const Int_t npTArrays      = 8;  //0-1, 1-1.5, 1.5-2, 2-2.5, 2.5-3, 3-3.5, 3.5-4, 4-7 GeV/c

  gStyle->SetTitleW(.8f);
  gStyle->SetTitleH(.08f);
  gStyle->SetTitleSize(.045,"xy");
  //gStyle->SetTitleOffset(.8,"x");
  gROOT->ForceStyle();

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Define Variables
  ///____________________________________
  Int_t Sigma                         = 3;                     // 2 or 3 is testing in this analysis

  char *AllOrAnti                     = "All";                 // which particles shall be analyse?
  //char *AllOrAnti                     = "Anti";

  //char *Collisions = "PbPb_11h_tracklets";
  char *Collisions                    = "PbPb_11h"; char *Tracklets = "";
  //char *Collisions                    = "PbPb_11hStd";
  //char *Collisions                    = "pPb_13b";
  //char *Collisions                    = "pPb_13c";
  //char *Collisions                    = "pp_15f";

  //char *Tracklets                    = "4trkl"; Double_t fnTrkl = 4; // exists only in PbPb_11h_tracklets
  //char *Tracklets                    = "5trkl"; Double_t fnTrkl = 5;
  //char *Tracklets                    = "6trkl"; Double_t fnTrkl = 6;

  char *Particle[nParticles]          = {"Deuteron","Triton","Helium3","Alpha","Helium3Alpha","DeuTriHe3Alp"};

  char *Preselection[nPreselections]  = {"TPC","TPCTOF"};

  char *Analyse[nAnalysis]            = {"Pur","Eff","Rej"};
  char *Analyselong[nAnalysis]        = {"Purity","Trigger efficiency","Rejection"};

  char *File[nFiles]                  = {"pdf","root"};

  Int_t pTArray500Anti_i[npTArrays] = {241,236,231,226,221,216,211,181}; // for hist 1 to 500 bins
  Int_t pTArray500Anti_f[npTArrays] = {250,240,235,230,225,220,215,210};
  Int_t pTArray500All_i[npTArrays]  = {251,261,266,271,276,281,286,291};
  Int_t pTArray500All_f[npTArrays]  = {260,265,270,275,280,285,290,320};
  Int_t pTArray160Anti_i[npTArrays] = { 71, 66, 61, 56, 51, 46, 41, 11}; // for hist 1 to 160 bins
  Int_t pTArray160Anti_f[npTArrays] = { 80, 70, 65, 60, 55, 50, 45, 40};
  Int_t pTArray160All_i[npTArrays]  = { 81, 91, 96,101,106,111,116,121};
  Int_t pTArray160All_f[npTArrays]  = { 90, 95,100,105,110,115,120,150};

  Double_t pTRange[npTArrays]         = {0.5,1.25,1.75,2.25,2.75,3.25,3.75,5.5};
  Double_t DrawY[npTArrays]           = {0};
  Double_t DraweX[npTArrays]          = {0.5,0.25,0.25,0.25,0.25,0.25,0.25,1.5};
  Double_t DraweY[npTArrays]          = {0};

  if(Sigma == 2) {
    Double_t ClearCutDeuteronTPC     = 1.5; // Rigidity cuts on TPCnSigma histograms
    Double_t ClearCutDeuteronTPCTOF  = 1.7;
    Double_t ClearCutTritonTPC       = 1.8;
    Double_t ClearCutTritonTPCTOF    = 2.;
  }
  else if(Sigma == 3) {
    Double_t ClearCutDeuteronTPC     = 1.4;
    Double_t ClearCutDeuteronTPCTOF  = 1.7;
    Double_t ClearCutTritonTPC       = 1.7;
    Double_t ClearCutTritonTPCTOF    = 2.;
  }
  else {
    Double_t ClearCutDeuteronTPC     = 20.;
    Double_t ClearCutDeuteronTPCTOF  = 20.;
    Double_t ClearCutTritonTPC       = 20.;
    Double_t ClearCutTritonTPCTOF    = 20.;
  }


  // PID histograms
  TH2D         *fHistTRDpTvPID[nParticles][nPreselections];
  TH2D         *fHistTRDpTvPIDAllexcept[nParticles][nPreselections];
  // Projection histograms
  TH1D         *fHistTRDPID[nParticles][nPreselections];
  TH1D         *fHistTRDPIDAllexcept[nParticles][nPreselections];
  TH1D         *fHistTRDPIDall[nParticles][nPreselections];

  // Analysis graphs
  TGraphErrors *fpTDepGraph[nParticles][nAnalysis][nPreselections][nPads];

  TCanvas      *c[nParticles][nAnalysis][nPreselections];
  TCanvas      *cBA[nAnalysis];

  TLegend      *legCollisions;
  TLegend      *legPIDCut[nPads];
  TLegend      *legRigCut[nParticles][nPreselections];
  TLegend      *legPreselection[nParticles][nPreselections];
  TLegend      *legTracklets = new TLegend(.54,.15,.92,.25);
  legTracklets->SetFillStyle(0);
  legTracklets->SetBorderSize(0);
  TLegend      *legTracklets2 = new TLegend(.05,.15,.45,.25);
  legTracklets2->SetFillStyle(0);
  legTracklets2->SetBorderSize(0);

  printf(Form("\nCreating and saving 36 canvases with Sigma = %i for %s particles from %s...\n",Sigma,AllOrAnti,Collisions));
  if(Tracklets != "") printf(Form("with %.f tracklet tracks\n\n",fnTrkl));

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Get data from file /lustre/nyx/alice/users/bbrudnyj/trunk/nmartin_nuclei/AliAnalysisTaskTestTriggerTRD.cxx
  /// Get histograms from /lustre/nyx/alice/users/bbrudnyj/codes/bbrudnyj_ReadTreeTestTRDTrigger2PID.C
  ///                 and                                    .../bbrudnyj_ReadTreeTestTRDTrigger3.C
  ///_________________________________________________________________________________________________
  TFile *cinput1     = TFile::Open(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger2/bbrudnyj_TestTRDTrigger2%s%iSigma%sPID.root",Collisions,Tracklets,Collisions,Sigma,Tracklets));
  TList *inlistPID   = (TList*)cinput1->Get(Form("bbrudnyj_TestTRDTrigger2%s%iSigma%sPID;1",Collisions,Sigma,Tracklets));
  TList *inlistProj  = (TList*)cinput1->Get(Form("bbrudnyj_TestTRDTrigger2%s%iSigma%sProj;1",Collisions,Sigma,Tracklets));

  TFile *cinput2     = TFile::Open(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger3/bbrudnyj_TestTRDTrigger3%s%iSigma%s.root",Collisions,Tracklets,Collisions,Sigma,Tracklets));
  TList *inlistSigma = (TList*)cinput2->Get(Form("bbrudnyj_TestTRDTrigger3%s%iSigma%s;1",Collisions,Sigma,Tracklets));

  // PID histograms
    TH2D *fHistTRDpTvPIDall         = (TH2D*)inlistPID->FindObject("fHistTRDpTvPID");

    if(AllOrAnti == "All") {
      // PID histograms
      fHistTRDpTvPID[0][0]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvDeuteronTPC_clear");
      fHistTRDpTvPID[1][0]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvTritonTPC_clear");
      fHistTRDpTvPID[2][0]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvHelium3TPC");
      fHistTRDpTvPID[3][0]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAlphaTPC");
      fHistTRDpTvPID[4][0]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvHelium3AlphaTPC");
      fHistTRDpTvPID[5][0]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvDeuTriHe3AlpTPC");

      fHistTRDpTvPIDAllexcept[0][0] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptDeuteronTPC_clear");
      fHistTRDpTvPIDAllexcept[1][0] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptTritonTPC_clear");
      fHistTRDpTvPIDAllexcept[2][0] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptHelium3TPC");
      fHistTRDpTvPIDAllexcept[3][0] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAlphaTPC");
      fHistTRDpTvPIDAllexcept[4][0] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptHelium3AlphaTPC");
      fHistTRDpTvPIDAllexcept[5][0] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptDeuTriHe3AlpTPC");

      fHistTRDpTvPID[0][1]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvDeuteronTPCTOF_clear");
      fHistTRDpTvPID[1][1]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvTritonTPCTOF_clear");
      fHistTRDpTvPID[2][1]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvHelium3TPCTOF");
      fHistTRDpTvPID[3][1]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAlphaTPCTOF");
      fHistTRDpTvPID[4][1]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvHelium3AlphaTPCTOF");
      fHistTRDpTvPID[5][1]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvDeuTriHe3AlpTPCTOF");

      fHistTRDpTvPIDAllexcept[0][1] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptDeuteronTPCTOF_clear");
      fHistTRDpTvPIDAllexcept[1][1] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptTritonTPCTOF_clear");
      fHistTRDpTvPIDAllexcept[2][1] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptHelium3TPCTOF");
      fHistTRDpTvPIDAllexcept[3][1] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAlphaTPCTOF");
      fHistTRDpTvPIDAllexcept[4][1] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptHelium3AlphaTPCTOF");
      fHistTRDpTvPIDAllexcept[5][1] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptDeuTriHe3AlpTPCTOF");

    // Projection histograms
      // TPC preselection
      fHistTRDPID[0][0]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPC_clear");
      fHistTRDPID[1][0]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPC_clear");
      fHistTRDPID[2][0]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3TPC");
      fHistTRDPID[3][0]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDAlphaTPC");
      fHistTRDPID[4][0]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3AlphaTPC");
      fHistTRDPID[5][0]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuTriHe3AlpTPC");

      fHistTRDPIDAllexcept[0][0]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptDeuteronTPC_clear");
      fHistTRDPIDAllexcept[1][0]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptTritonTPC_clear");
      fHistTRDPIDAllexcept[2][0]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptHelium3TPC");
      fHistTRDPIDAllexcept[3][0]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAlphaTPC");
      fHistTRDPIDAllexcept[4][0]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptHelium3AlphaTPC");
      fHistTRDPIDAllexcept[5][0]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptDeuTriHe3AlpTPC");

      TH1D *fHistTRDPIDallTPC       = (TH1D*)inlistProj->FindObject("fHistTRDPIDallTPC");
      fHistTRDPIDall[0][0]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallDeuteronTPC");
      fHistTRDPIDall[1][0]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallTritonTPC");
      fHistTRDPIDall[2][0]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallHelium3TPC");
      fHistTRDPIDall[3][0]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAlphaTPC");
      fHistTRDPIDall[4][0]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallHelium3AlphaTPC");
      fHistTRDPIDall[5][0]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallDeuTriHe3AlpTPC");

      TH1D *fHistTRDPIDDeuteronTPC_clear_Clone     = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPC_clear_Clone");
      TH1D *fHistTRDPIDDeuteronTPC_clear_Clone2    = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPC_clear_Clone2");
      TH1D *fHistTRDPIDTritonTPC_clear_Clone       = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPC_clear_Clone");
      TH1D *fHistTRDPIDTritonTPC_clear_Clone2      = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPC_clear_Clone2");
      TH1D *fHistTRDPIDHelium3TPC_Clone            = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3TPC_Clone");
      TH1D *fHistTRDPIDHelium3TPC_Clone2           = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3TPC_Clone2");
      TH1D *fHistTRDPIDAlphaTPC_Clone              = (TH1D*)inlistProj->FindObject("fHistTRDPIDAlphaTPC_Clone");
      TH1D *fHistTRDPIDAlphaTPC_Clone2             = (TH1D*)inlistProj->FindObject("fHistTRDPIDAlphaTPC_Clone2");
      TH1D *fHistTRDPIDHelium3AlphaTPC_Clone       = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3AlphaTPC_Clone");
      TH1D *fHistTRDPIDDeuTriHe3AlpTPC_Clone       = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuTriHe3AlpTPC_Clone");

      // TPC & TOF preselection
      fHistTRDPID[0][1]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPCTOF_clear");
      fHistTRDPID[1][1]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPCTOF_clear");
      fHistTRDPID[2][1]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3TPCTOF");
      fHistTRDPID[3][1]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDAlphaTPCTOF");
      fHistTRDPID[4][1]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3AlphaTPCTOF");
      fHistTRDPID[5][1]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuTriHe3AlpTPCTOF");

      fHistTRDPIDAllexcept[0][1]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptDeuteronTPCTOF_clear");
      fHistTRDPIDAllexcept[1][1]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptTritonTPCTOF_clear");
      fHistTRDPIDAllexcept[2][1]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptHelium3TPCTOF");
      fHistTRDPIDAllexcept[3][1]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAlphaTPCTOF");
      fHistTRDPIDAllexcept[4][1]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptHelium3AlphaTPCTOF");
      fHistTRDPIDAllexcept[5][1]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptDeuTriHe3AlpTPCTOF");

      TH1D *fHistTRDPIDallTPCTOF    = (TH1D*)inlistProj->FindObject("fHistTRDPIDallTPCTOF");
      fHistTRDPIDall[0][1]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallDeuteronTPCTOF");
      fHistTRDPIDall[1][1]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallTritonTPCTOF");
      fHistTRDPIDall[2][1]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallHelium3TPCTOF");
      fHistTRDPIDall[3][1]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAlphaTPCTOF");
      fHistTRDPIDall[4][1]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallHelium3AlphaTPCTOF");
      fHistTRDPIDall[5][1]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallDeuTriHe3AlpTPCTOF");

      TH1D *fHistTRDPIDDeuteronTPCTOF_clear_Clone  = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPCTOF_clear_Clone");
      TH1D *fHistTRDPIDDeuteronTPCTOF_clear_Clone2 = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPCTOF_clear_Clone2");
      TH1D *fHistTRDPIDTritonTPCTOF_clear_Clone    = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPCTOF_clear_Clone");
      TH1D *fHistTRDPIDTritonTPCTOF_clear_Clone2   = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPCTOF_clear_Clone2");
      TH1D *fHistTRDPIDHelium3TPCTOF_Clone         = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3TPCTOF_Clone");
      TH1D *fHistTRDPIDHelium3TPCTOF_Clone2        = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3TPCTOF_Clone2");
      TH1D *fHistTRDPIDAlphaTPCTOF_Clone           = (TH1D*)inlistProj->FindObject("fHistTRDPIDAlphaTPCTOF_Clone");
      TH1D *fHistTRDPIDAlphaTPCTOF_Clone2          = (TH1D*)inlistProj->FindObject("fHistTRDPIDAlphaTPCTOF_Clone2");
      TH1D *fHistTRDPIDHelium3AlphaTPCTOF_Clone    = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3AlphaTPCTOF_Clone");
      TH1D *fHistTRDPIDDeuTriHe3AlpTPCTOF_Clone    = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuTriHe3AlpTPCTOF_Clone");

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
      // PID histograms
      fHistTRDpTvPID[0][0]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiDeuteronTPC_clear");
      fHistTRDpTvPID[1][0]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiTritonTPC_clear");
      fHistTRDpTvPID[2][0]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiHelium3TPC");
      fHistTRDpTvPID[3][0]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiAlphaTPC");
      fHistTRDpTvPID[4][0]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiHelium3AlphaTPC");
      fHistTRDpTvPID[5][0]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiDeuTriHe3AlpTPC");

      fHistTRDpTvPIDAllexcept[0][0] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiDeuteronTPC_clear");
      fHistTRDpTvPIDAllexcept[1][0] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiTritonTPC_clear");
      fHistTRDpTvPIDAllexcept[2][0] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiHelium3TPC");
      fHistTRDpTvPIDAllexcept[3][0] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiAlphaTPC");
      fHistTRDpTvPIDAllexcept[4][0] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiHelium3AlphaTPC");
      fHistTRDpTvPIDAllexcept[5][0] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiDeuTriHe3AlpTPC");

      fHistTRDpTvPID[0][1]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiDeuteronTPCTOF_clear");
      fHistTRDpTvPID[1][1]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiTritonTPCTOF_clear");
      fHistTRDpTvPID[2][1]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiHelium3TPCTOF");
      fHistTRDpTvPID[3][1]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiAlphaTPCTOF");
      fHistTRDpTvPID[4][1]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiHelium3AlphaTPCTOF");
      fHistTRDpTvPID[5][1]          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiDeuTriHe3AlpTPCTOF");

      fHistTRDpTvPIDAllexcept[0][1] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiDeuteronTPCTOF_clear");
      fHistTRDpTvPIDAllexcept[1][1] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiTritonTPCTOF_clear");
      fHistTRDpTvPIDAllexcept[2][1] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiHelium3TPCTOF");
      fHistTRDpTvPIDAllexcept[3][1] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiAlphaTPCTOF");
      fHistTRDpTvPIDAllexcept[4][1] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiHelium3AlphaTPCTOF");
      fHistTRDpTvPIDAllexcept[5][1] = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAllexceptAntiDeuTriHe3AlpTPCTOF");

    // Projection histograms
      // TPC preselection
      fHistTRDPID[0][0]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuteronTPC_clear");
      fHistTRDPID[1][0]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiTritonTPC_clear");
      fHistTRDPID[2][0]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3TPC");
      fHistTRDPID[3][0]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiAlphaTPC");
      fHistTRDPID[4][0]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3AlphaTPC");
      fHistTRDPID[5][0]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuTriHe3AlpTPC");

      fHistTRDPIDAllexcept[0][0]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiDeuteronTPC_clear");
      fHistTRDPIDAllexcept[1][0]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiTritonTPC_clear");
      fHistTRDPIDAllexcept[2][0]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiHelium3TPC");
      fHistTRDPIDAllexcept[3][0]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiAlphaTPC");
      fHistTRDPIDAllexcept[4][0]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiHelium3AlphaTPC");
      fHistTRDPIDAllexcept[5][0]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiDeuTriHe3AlpTPC");

      fHistTRDPIDall[0][0]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiDeuteronTPC");
      fHistTRDPIDall[1][0]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiTritonTPC");
      fHistTRDPIDall[2][0]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiHelium3TPC");
      fHistTRDPIDall[3][0]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiAlphaTPC");
      fHistTRDPIDall[4][0]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiHelium3AlphaTPC");
      fHistTRDPIDall[5][0]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiDeuTriHe3AlpTPC");

      TH1D *fHistTRDPIDDeuteronTPC_clear_Clone     = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuteronTPC_clear_Clone");
      TH1D *fHistTRDPIDDeuteronTPC_clear_Clone2    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuteronTPC_clear_Clone2");
      TH1D *fHistTRDPIDTritonTPC_clear_Clone       = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiTritonTPC_clear_Clone");
      TH1D *fHistTRDPIDTritonTPC_clear_Clone2      = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiTritonTPC_clear_Clone2");
      TH1D *fHistTRDPIDHelium3TPC_Clone            = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3TPC_Clone");
      TH1D *fHistTRDPIDHelium3TPC_Clone2           = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3TPC_Clone2");
      TH1D *fHistTRDPIDAlphaTPC_Clone              = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiAlphaTPC_Clone");
      TH1D *fHistTRDPIDAlphaTPC_Clone2             = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiAlphaTPC_Clone2");
      TH1D *fHistTRDPIDHelium3AlphaTPC_Clone       = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3AlphaTPC_Clone");
      TH1D *fHistTRDPIDDeuTriHe3AlpTPC_Clone       = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuTriHe3AlpTPC_Clone");

      // TPC & TOF preselection
      fHistTRDPID[0][1]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuteronTPCTOF_clear");
      fHistTRDPID[1][1]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiTritonTPCTOF_clear");
      fHistTRDPID[2][1]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3TPCTOF");
      fHistTRDPID[3][1]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiAlphaTPCTOF");
      fHistTRDPID[4][1]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3AlphaTPCTOF");
      fHistTRDPID[5][1]             = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuTriHe3AlpTPCTOF");

      fHistTRDPIDAllexcept[0][1]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiDeuteronTPCTOF_clear");
      fHistTRDPIDAllexcept[1][1]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiTritonTPCTOF_clear");
      fHistTRDPIDAllexcept[2][1]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiHelium3TPCTOF");
      fHistTRDPIDAllexcept[3][1]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiAlphaTPCTOF");
      fHistTRDPIDAllexcept[4][1]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiHelium3AlphaTPCTOF");
      fHistTRDPIDAllexcept[5][1]    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAllexceptAntiDeuTriHe3AlpTPCTOF");

      fHistTRDPIDall[0][1]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiDeuteronTPCTOF");
      fHistTRDPIDall[1][1]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiTritonTPCTOF");
      fHistTRDPIDall[2][1]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiHelium3TPCTOF");
      fHistTRDPIDall[3][1]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiAlphaTPCTOF");
      fHistTRDPIDall[4][1]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiHelium3AlphaTPCTOF");
      fHistTRDPIDall[5][1]          = (TH1D*)inlistProj->FindObject("fHistTRDPIDallAntiDeuTriHe3AlpTPCTOF");

      TH1D *fHistTRDPIDDeuteronTPCTOF_clear_Clone  = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuteronTPCTOF_clear_Clone");
      TH1D *fHistTRDPIDDeuteronTPCTOF_clear_Clone2 = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuteronTPCTOF_clear_Clone2");
      TH1D *fHistTRDPIDTritonTPCTOF_clear_Clone    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiTritonTPCTOF_clear_Clone");
      TH1D *fHistTRDPIDTritonTPCTOF_clear_Clone2   = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiTritonTPCTOF_clear_Clone2");
      TH1D *fHistTRDPIDHelium3TPCTOF_Clone         = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3TPCTOF_Clone");
      TH1D *fHistTRDPIDHelium3TPCTOF_Clone2        = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3TPCTOF_Clone2");
      TH1D *fHistTRDPIDAlphaTPCTOF_Clone           = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiAlphaTPCTOF_Clone");
      TH1D *fHistTRDPIDAlphaTPCTOF_Clone2          = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiAlphaTPCTOF_Clone2");
      TH1D *fHistTRDPIDHelium3AlphaTPCTOF_Clone    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3AlphaTPCTOF_Clone");
      TH1D *fHistTRDPIDDeuTriHe3AlpTPCTOF_Clone    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuTriHe3AlpTPCTOF_Clone");

    // nSigma histograms
      //TH2D *fHistTPCnSigmavRigvDeuteronTPC_clear             = (TH2D*)inlistSigma->FindObject("fHistTPCnSigmavRigvAntiDeuteronTPC_clear");
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
    Double_t nMinBin                                           = 260. - n*10;
    Double_t nMeanQRange[n]                                    = {0};

    Int_t    nPID[npTArrays][nParticles][nPreselections][n]    = {{{{0}}}};
    Int_t    nPIDall[npTArrays][n]                             = {{0}};
    Double_t Purity[npTArrays][nParticles][nPreselections][n]  = {{{{0}}}};

    Double_t ePID[npTArrays][nParticles][nPreselections][n]    = {{{{0}}}};
    Double_t ePIDall[npTArrays][n]                             = {{0}};
    Double_t ePurity[npTArrays][nParticles][nPreselections][n] = {{{{0}}}};

    for(Int_t p = 0; p < npTArrays; p++) {
      for(Int_t i = 0; i < n; i++) {
        nMeanQRange[i] = nMinBin + i*10;
        nPIDall[p][i]                                        = fHistTRDpTvPIDall->IntegralAndError(pTArray500Anti_i[p],pTArray500Anti_f[p],nMeanQRange[i],260,ePIDall[p][i]);
	if(AllOrAnti == "All") nPIDall[p][i] = nPIDall[p][i] + fHistTRDpTvPIDall->IntegralAndError(pTArray500All_i[p] ,pTArray500All_f[p] ,nMeanQRange[i],260,ePIDall[p][i]);


        if(nPIDall[p][i] == 0) continue;
	if(nPIDall[p][i] == 1) ePIDall[p][i] = 1;

        for(Int_t j = 0; j < nParticles; j++) {
	  for(Int_t k = 0; k < nPreselections; k++) {

            nPID[p][j][k][i]                                           = fHistTRDpTvPID[j][k]->IntegralAndError(pTArray160Anti_i[p],pTArray160Anti_f[p],nMeanQRange[i],260,ePID[p][j][k][i]);
	    if(AllOrAnti == "All") nPID[p][j][k][i] = nPID[p][j][k][i] + fHistTRDpTvPID[j][k]->IntegralAndError(pTArray160All_i[p] ,pTArray160All_f[p] ,nMeanQRange[i],260,ePID[p][j][k][i]);

            if(nPID[p][j][k][i] == 0) continue;
	    if(nPID[p][j][k][i] == 1) ePID[p][j][k][i] = 1;

	    Purity[p][j][k][i]  = (Double_t)nPID[p][j][k][i] / (Double_t)nPIDall[p][i];
            ePurity[p][j][k][i] = GaussError((Double_t)nPID[p][j][k][i],ePID[p][j][k][i],(Double_t)nPIDall[p][i],ePIDall[p][i]); //own function, see below


	    //if(p==7 && i==10 && j==4 && k==0){
            //  cout << Form("pTArray: %.2f, PID > %.f, Particle: %s, Preselection: %s\n",pTRange[p],nMeanQRange[i],Particle[j],Preselection[k]) << endl;

            //  cout << "nPID = " << nPID[p][j][k][i] << endl;
            //  cout << "ePID = " << ePID[p][j][k][i] << endl;
            //  cout << "nPIDall = " << nPIDall[p][i] << endl;
            //  cout << "ePIDall = " << ePIDall[p][i] << "\n" << endl;

	    //  cout <<  Purity[p][j][k][i] << " err " << ePurity[p][j][k][i] << "\n" << endl;
	    //}
          }
        }
      }
    }



    // fill the graphs
    for(Int_t j = 0; j < nParticles; j++) {
      for(Int_t k = 0; k < nPreselections; k++) {
	for(Int_t l = 2; l < nPads; l++) {
          for(Int_t p = 0; p < npTArrays; p++) {
	    if(l==2) {
	      DrawY[p]  = 0. + Purity[p][j][k][l-2];
	      DraweY[p] = 0. + ePurity[p][j][k][l-2];
	    }
	    else {
	      DrawY[p]  = 0. + Purity[p][j][k][l-2+8];
	      DraweY[p] = 0. + ePurity[p][j][k][l-2+8];
	    }
	  }
          fpTDepGraph[j][0][k][l] = new TGraphErrors(npTArrays,pTRange,DrawY,DraweX,DraweY);
          if(AllOrAnti == "All") fpTDepGraph[j][0][k][l]->SetTitle(Form("%s purity (%i#sigma %s preselection);#it{p}_{T} (GeV/#it{c});Purity",Particle[j],Sigma,Preselection[k]));
	  if(AllOrAnti == "Anti")fpTDepGraph[j][0][k][l]->SetTitle(Form("Anti-%s purity (%i#sigma %s preselection);#it{p}_{T} (GeV/#it{c});Purity",Particle[j],Sigma,Preselection[k]));
          fpTDepGraph[j][0][k][l]->SetMarkerStyle(7);
          fpTDepGraph[j][0][k][l]->SetMarkerColor(4);
	  fpTDepGraph[j][0][k][l]->SetMinimum(0);
        }
      }
    }

  //___________
  // efficiency
    Double_t mMinBin                                              = 260. - m*10;
    Double_t mMeanQRange[m]                                       = {0};

    Int_t nPIDEntries[npTArrays][nParticles][nPreselections][m]    = {{{{0}}}};
    Int_t nSigmaEntries[npTArrays][nParticles][nPreselections]     = {{{0}}};
    Double_t Efficiency[npTArrays][nParticles][nPreselections][m]  = {{{{0}}}};

    Double_t ePIDEntries[npTArrays][nParticles][nPreselections][m] = {{{{0}}}};
    Double_t eSigmaEntries[npTArrays][nParticles][nPreselections]  = {{{0}}};
    Double_t eEfficiency[npTArrays][nParticles][nPreselections][m] = {{{{0}}}};

    for(Int_t p = 0; p < npTArrays; p++) {
      for(Int_t j = 0; j < nParticles; j++) {
        for(Int_t k = 0; k < nPreselections; k++) {
          nSigmaEntries[p][j][k]                                                 = fHistTRDpTvPID[j][k]->IntegralAndError(pTArray160Anti_i[p],pTArray160Anti_f[p],1,260,eSigmaEntries[p][j][k]); // not correct
          if(AllOrAnti == "All") nSigmaEntries[p][j][k] = nSigmaEntries[p][j][k] + fHistTRDpTvPID[j][k]->IntegralAndError(pTArray160All_i[p] ,pTArray160All_f[p] ,1,260,eSigmaEntries[p][j][k]);

          if(nSigmaEntries[p][j][k] == 0) continue;
          if(nSigmaEntries[p][j][k] == 1) eSigmaEntries[p][j][k] = 1;

          for(Int_t i = 0; i < m; i++) {
            mMeanQRange[i] = mMinBin + i*10;

            nPIDEntries[p][j][k][i]                                                = fHistTRDpTvPID[j][k]->IntegralAndError(pTArray160Anti_i[p],pTArray160Anti_f[p],mMeanQRange[i],260,ePIDEntries[p][j][k][i]);
            if(AllOrAnti == "All") nPIDEntries[p][j][k][i] = nPIDEntries[p][j][k][i] + fHistTRDpTvPID[j][k]->IntegralAndError(pTArray160All_i[p],pTArray160All_f[p],mMeanQRange[i],260,ePIDEntries[p][j][k][i]);

            if(nPIDEntries[p][j][k][i] == 0) continue;
            if(nPIDEntries[p][j][k][i] == 1) ePIDEntries[p][j][k][i] = 1;

            Efficiency[p][j][k][i]  = (Double_t)nPIDEntries[p][j][k][i] / (Double_t)nSigmaEntries[p][j][k];
            eEfficiency[p][j][k][i] = GaussError((Double_t)nPIDEntries[p][j][k][i],ePIDEntries[p][j][k][i],(Double_t)nSigmaEntries[p][j][k],eSigmaEntries[p][j][k]);

	    //if(p==7 && j==2 && k==0 && i==10){
            //  cout << Form("\npTArray: %.2f, PID > %.f, Particle: %s, Preselection: %s\n",pTRange[p],mMeanQRange[i],Particle[j],Preselection[k]) << endl;

            //  cout << "nPIDEntries = " << nPIDEntries[p][j][k][i] << endl;
            //  cout << "ePIDEntries = " << ePIDEntries[p][j][k][i] << endl;
            //  cout << "nSigmaEntries = " << nSigmaEntries[p][j][k] << endl;
            //  cout << "eSigmaEntries = " << eSigmaEntries[p][j][k] << "\n" << endl;

	    //  cout <<  Efficiency[p][j][k][i] << " err " << eEfficiency[p][j][k][i] << "\n" << endl;
	    //}
	  }
	}
      }
    }

    // fill the graphs
    for(Int_t j = 0; j < nParticles; j++) {
      for(Int_t k = 0; k < nPreselections; k++) {
	for(Int_t l = 2; l < nPads; l++) {
          for(Int_t p = 0; p < npTArrays; p++) {
	    if(l==2) {
	      DrawY[p]  = 0. + Efficiency[p][j][k][l-2];
	      DraweY[p] = 0. + eEfficiency[p][j][k][l-2];
	    }
	    else {
	      DrawY[p]  = 0. + Efficiency[p][j][k][l-2+8];
	      DraweY[p] = 0. + eEfficiency[p][j][k][l-2+8];
	    }
	  }
          fpTDepGraph[j][1][k][l] = new TGraphErrors(npTArrays,pTRange,DrawY,DraweX,DraweY);
          if(AllOrAnti == "All") fpTDepGraph[j][1][k][l]->SetTitle(Form("%s efficiency (%i#sigma %s preselection);#it{p}_{T} (GeV/#it{c});Efficiency",Particle[j],Sigma,Preselection[k]));
	  if(AllOrAnti == "Anti")fpTDepGraph[j][1][k][l]->SetTitle(Form("Anti-%s efficiency (%i#sigma %s preselection);#it{p}_{T} (GeV/#it{c});Efficiency",Particle[j],Sigma,Preselection[k]));
          fpTDepGraph[j][1][k][l]->SetMarkerStyle(7);
          fpTDepGraph[j][1][k][l]->SetMarkerColor(2);
          if(l>2) fpTDepGraph[j][1][k][l]->SetMinimum(0);
        }
      }
    }


  //__________
  // rejection
    Int_t nPIDEntriesAllexcept[npTArrays][nParticles][nPreselections][m]    = {{{{0}}}};
    Int_t nSigmaEntriesAllexcept[npTArrays][nParticles][nPreselections]     = {{{0}}};
    Double_t Rejection[npTArrays][nParticles][nPreselections][m]            = {{{{0}}}};

    Double_t ePIDEntriesAllexcept[npTArrays][nParticles][nPreselections][m] = {{{{0}}}};
    Double_t eSigmaEntriesAllexcept[npTArrays][nParticles][nPreselections]  = {{{0}}};
    Double_t eRejection[npTArrays][nParticles][nPreselections][m]           = {{{{0}}}};

    for(Int_t p = 0; p < npTArrays; p++) {
      for(Int_t j = 0; j < nParticles; j++) {
        for(Int_t k = 0; k < nPreselections; k++) {
          nSigmaEntriesAllexcept[p][j][k]                        = fHistTRDpTvPIDAllexcept[j][k]->IntegralAndError(pTArray500Anti_i[p],pTArray500Anti_f[p],1,260,eSigmaEntriesAllexcept[p][j][k]); // not correct
          if(AllOrAnti == "All") nSigmaEntriesAllexcept[p][j][k] = nSigmaEntriesAllexcept[p][j][k] + fHistTRDpTvPIDAllexcept[j][k]->IntegralAndError(pTArray500All_i[p],pTArray500All_f[p],1,260,eSigmaEntriesAllexcept[p][j][k]);

          if(nSigmaEntriesAllexcept[p][j][k] == 0) continue;
          if(nSigmaEntriesAllexcept[p][j][k] == 1) eSigmaEntriesAllexcept[p][j][k] = 1;

          for(Int_t i = 0; i < m; i++) {
            mMeanQRange[i] = mMinBin + i*10;

            nPIDEntriesAllexcept[p][j][k][i]                     = fHistTRDpTvPIDAllexcept[j][k]->IntegralAndError(pTArray500Anti_i[p],pTArray500Anti_f[p],mMeanQRange[i],260,ePIDEntriesAllexcept[p][j][k][i]);
            if(AllOrAnti == "All") nPIDEntriesAllexcept[p][j][k][i] = nPIDEntriesAllexcept[p][j][k][i] + fHistTRDpTvPIDAllexcept[j][k]->IntegralAndError(pTArray500All_i[p],pTArray500All_f[p],mMeanQRange[i],260,ePIDEntriesAllexcept[p][j][k][i]);

            if(nPIDEntriesAllexcept[p][j][k][i] == 0) continue;
            if(nPIDEntriesAllexcept[p][j][k][i] == 1) ePIDEntriesAllexcept[p][j][k][i] = 1;

            Rejection[p][j][k][i]  = (Double_t)nPIDEntriesAllexcept[p][j][k][i] / (Double_t)nSigmaEntriesAllexcept[p][j][k];
            eRejection[p][j][k][i] = GaussError((Double_t)nPIDEntriesAllexcept[p][j][k][i],ePIDEntriesAllexcept[p][j][k][i],(Double_t)nSigmaEntriesAllexcept[p][j][k],eSigmaEntriesAllexcept[p][j][k]);

	  }
	}
      }
    }

    // fill the graphs
    for(Int_t j = 0; j < nParticles; j++) {
      for(Int_t k = 0; k < nPreselections; k++) {
	for(Int_t l = 2; l < nPads; l++) {
          for(Int_t p = 0; p < npTArrays; p++) {
	    if(l==2) {
	      DrawY[p]  = 0. + Rejection[p][j][k][l-2];
	      DraweY[p] = 0. + eRejection[p][j][k][l-2];
	    }
	    else {
	      DrawY[p]  = 0. + Rejection[p][j][k][l-2+8];
	      DraweY[p] = 0. + eRejection[p][j][k][l-2+8];
	    }
	  }
          fpTDepGraph[j][2][k][l] = new TGraphErrors(npTArrays,pTRange,DrawY,DraweX,DraweY);
          if(AllOrAnti == "All") fpTDepGraph[j][2][k][l]->SetTitle(Form("%s rejection (%i#sigma %s preselection);#it{p}_{T} (GeV/#it{c});Rejection",Particle[j],Sigma,Preselection[k]));
	  if(AllOrAnti == "Anti")fpTDepGraph[j][2][k][l]->SetTitle(Form("Anti-%s rejection (%i#sigma %s preselection);#it{p}_{T} (GeV/#it{c});Rejection",Particle[j],Sigma,Preselection[k]));
          fpTDepGraph[j][2][k][l]->SetMarkerStyle(7);
          fpTDepGraph[j][2][k][l]->SetMarkerColor(3);
        }
      }
    }



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Create Legends
  ///_________________
  // legends for cd(1)
  if(Tracklets != "") legCollisions = new TLegend(.02,.6,.45,.88);
  else                legCollisions = new TLegend(.02,.68,.45,.88);
  legCollisions->SetFillStyle(0);
  legCollisions->SetBorderSize(0);
  if(Collisions == "PbPb_11h" || Collisions == "PbPb_11hStd" || Collisions == "PbPb_11h_tracklets") {
    legCollisions->AddEntry((TObject*)0,"Pb-Pb, #sqrt{s_{NN}} = 2.76 TeV","");
    legCollisions->AddEntry((TObject*)0,"run LHC11h pass2","");
  }
  else if(Collisions == "pPb_13b") {
    legCollisions->AddEntry((TObject*)0,"p-Pb, #sqrt{s_{NN}} = 5.023 TeV","");
    legCollisions->AddEntry((TObject*)0,"run LHC13b pass3","");
  }
  else if(Collisions == "pPb_13c") {
    legCollisions->AddEntry((TObject*)0,"p-Pb, #sqrt{s_{NN}} = 5.023 TeV","");
    legCollisions->AddEntry((TObject*)0,"run LHC13c pass2","");
  }
  else { //else if(Collisions == "pp_15f") {
    legCollisions->AddEntry((TObject*)0,"p-p, #sqrt{s} = 13 TeV","");
    legCollisions->AddEntry((TObject*)0,"run LHC15f pass2","");
  }
  if(Tracklets != "") {
    legCollisions->AddEntry((TObject*)0,Form("%.f tracklet tracks",fnTrkl),"");
    legTracklets->AddEntry((TObject*)0,Form("%.f tracklet tracks",fnTrkl),"");
    legTracklets2->AddEntry((TObject*)0,Form("%.f tracklet tracks",fnTrkl),"");
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Fill 36 canvases
  ///____________________________________
  for(Int_t j = 0; j < nParticles; j++) {
    for(Int_t ana = 0; ana < nAnalysis; ana++) {
      for(Int_t k = 0; k < nPreselections; k++) {

//        c[j][ana][k] = new TCanvas(Form("c%i%i%i",j,ana,k),Form("%s %s %s pT dependent with %iSigma %s preselection",AllOrAnti,Particle[j],Analyse[ana],Sigma,Preselection[k]));
//	c[j][ana][k]->Divide(3,3);
	// cd(1)
//	c[j][ana][k]->cd(1)->SetLogz();
//          fHistTRDpTvPIDall->Draw("colz");
//	  legCollisions->Draw("same");

	// cd(2)
	if(j<2  && k==0) legPreselection[j][k] = new TLegend(.05,.59,.45,.85);
	else if(j<2  && k==1) legPreselection[j][k] = new TLegend(.05,.55,.45,.85);
	else if(j>=2 && j<5 && k==0) legPreselection[j][k] = new TLegend(.05,.66,.45,.85);
	else if(j>=2 && j<5 && k==1) legPreselection[j][k] = new TLegend(.05,.6,.45,.85);
	else if(j==5 && k==0) legPreselection[j][k] = new TLegend(.05,.54,.45,.85);
	//else if(j==5 && k==1) legPreselection[j][k] = new TLegend(.05,.5,.45,.85);
	else if(j==5 && k==1) legPreselection[j][k] = new TLegend(.42,.32,.82,.74);
	legPreselection[j][k]->SetFillStyle(0);
        legPreselection[j][k]->SetBorderSize(0);
	if(j==0  && k==0) {
	  legPreselection[j][k]->AddEntry((TObject*)0,"with preselection:","");
	  legPreselection[j][k]->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
	  legPreselection[j][k]->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutDeuteronTPC),"");
	}
	else if(j==0  && k==1) {
	  legPreselection[j][k]->AddEntry((TObject*)0,"with preselection:","");
	  legPreselection[j][k]->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
	  legPreselection[j][k]->AddEntry((TObject*)0,Form("|TOFnSigma| = %i#sigma",Sigma),"");
	  legPreselection[j][k]->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutDeuteronTPCTOF),"");
	}
	else if(j==1  && k==0) {
	  legPreselection[j][k]->AddEntry((TObject*)0,"with preselection:","");
	  legPreselection[j][k]->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
	  legPreselection[j][k]->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutTritonTPC),"");
	}
	else if(j==1  && k==1) {
	  legPreselection[j][k]->AddEntry((TObject*)0,"with preselection:","");
	  legPreselection[j][k]->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
	  legPreselection[j][k]->AddEntry((TObject*)0,Form("|TOFnSigma| = %i#sigma",Sigma),"");
	  legPreselection[j][k]->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutTritonTPCTOF),"");
	}
	else if(j>=2 && j<5 && k==0) {
	  legPreselection[j][k]->AddEntry((TObject*)0,"with preselection:","");
	  legPreselection[j][k]->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
	}
	else if(j>=2 && j<5 && k==1) {
	  legPreselection[j][k]->AddEntry((TObject*)0,"with preselection:","");
	  legPreselection[j][k]->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
	  legPreselection[j][k]->AddEntry((TObject*)0,Form("|TOFnSigma| = %i#sigma",Sigma),"");
	}
	else if(j==5  && k==0) {
	  legPreselection[j][k]->AddEntry((TObject*)0,"with preselection:","");
	  legPreselection[j][k]->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
	  legPreselection[j][k]->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c} for d",ClearCutDeuteronTPC),"");
	  legPreselection[j][k]->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c} for t",ClearCutTritonTPC),"");
	}
	else if(j==5  && k==1) {
	  legPreselection[j][k]->AddEntry((TObject*)0,"with preselection:","");
	  legPreselection[j][k]->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
	  legPreselection[j][k]->AddEntry((TObject*)0,Form("|TOFnSigma| = %i#sigma",Sigma),"");
	  legPreselection[j][k]->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c} for d",ClearCutDeuteronTPCTOF),"");
	  legPreselection[j][k]->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c} for t",ClearCutTritonTPCTOF),"");
	}

//	c[j][ana][k]->cd(2)->SetLogz();
//	  fHistTRDpTvPID[j][k]->Draw("colz");
//	  legPreselection[j][k]->Draw("same");
//	  legTracklets2->Draw("same");

        // cd(3-9)
        if(j<5)  legRigCut[j][k] = new TLegend(.38,.6,.82,.72);
	else if(j==5) legRigCut[j][k] = new TLegend(.35,.52,.89,.72);
	legRigCut[j][k]->SetFillStyle(0);
        legRigCut[j][k]->SetBorderSize(0);
	if(j==0 && k==0) legRigCut[j][k]->AddEntry((TObject*)0,Form("Rigidity cut for d: |#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutDeuteronTPC),"");
	else if(j==0 && k==1) legRigCut[j][k]->AddEntry((TObject*)0,Form("Rigidity cut for d: |#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutDeuteronTPCTOF),"");
	else if(j==1 && k==0) legRigCut[j][k]->AddEntry((TObject*)0,Form("Rigidity cut for t: |#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutTritonTPC),"");
	else if(j==1 && k==1) legRigCut[j][k]->AddEntry((TObject*)0,Form("Rigidity cut for t: |#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutTritonTPCTOF),"");
	else if(j==5 && k==0) {
	  legRigCut[j][k]->AddEntry((TObject*)0,Form("Rigidity cut for d: |#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutDeuteronTPC),"");
          legRigCut[j][k]->AddEntry((TObject*)0,Form("Rigidity cut for t: |#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutTritonTPC),"");
	}
	else if(j==5 && k==1) {
	  legRigCut[j][k]->AddEntry((TObject*)0,Form("Rigidity cut for d: |#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutDeuteronTPCTOF),"");
          legRigCut[j][k]->AddEntry((TObject*)0,Form("Rigidity cut for t: |#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutTritonTPCTOF),"");
	}

	for(Int_t l = 2; l < nPads; l++) { // for pad 3-9 for every canvas

	  legPIDCut[l] = new TLegend(.3,.72,.9,.87);
          legPIDCut[l]->SetFillStyle(0);
          legPIDCut[l]->SetBorderSize(0);
	  if(l==2) legPIDCut[l]->AddEntry((TObject*)0,"Threshold at #bf{<#it{Q}_{s}> = 0 a.u.}","");
	  else     legPIDCut[l]->AddEntry((TObject*)0,Form("Threshold at #bf{<#it{Q}_{s}> = %i a.u.}",60+l*10),"");


//	  c[j][ana][k]->cd(l+1)->SetTicks();
//	  fpTDepGraph[j][ana][k][l]->Draw("AP");
//	  legPIDCut[l]->Draw("same");
//	  legRigCut[j][k]->Draw("same");
//	  legTracklets->Draw("same");
	}

        TLegend *legBIG = new TLegend(.36,.55,.96,.85);
	legBIG->SetFillStyle(0);
	legBIG->SetBorderSize(0);
	legBIG->AddEntry((TObject*)0,"Threshold at","");
	legBIG->AddEntry((TObject*)0,"#bf{#LT#it{Q}_{s}#GT} > #bf{100 a.u.}","");


	/*  //////////////////////////////////////////////////////////////////////////////////////////////
        /// Save 36 canvases
        ///________________________________
	for(Int_t f = 0; f < nFiles; f++) {

	c[j][ana][k]->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTriggerpTDep/plots%iSigma/%s/bbrudnyj_pTDep%s%s%s%s%iSigma%s.%s",Collisions,Tracklets,Sigma,AllOrAnti,AllOrAnti,Particle[j],Analyse[ana],Preselection[k],Sigma,Tracklets,File[f]));

        }*/
      }
    }
  }

  printf("\n\ndone!\n");


  for(Int_t ana = 0; ana < nAnalysis; ana++){

    fpTDepGraph[5][ana][1][4]->SetTitle("");
    fpTDepGraph[5][ana][1][4]->GetYaxis()->SetTitle(Form("%s",Analyselong[ana]));
    fpTDepGraph[5][ana][1][4]->GetYaxis()->SetTitleSize(.07);
    fpTDepGraph[5][ana][1][4]->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fpTDepGraph[5][ana][1][4]->GetXaxis()->SetTitleSize(.07);
    fpTDepGraph[5][ana][1][4]->SetMarkerStyle(34);
    fpTDepGraph[5][ana][1][4]->SetMarkerSize(2);
    fpTDepGraph[5][ana][1][4]->SetMarkerColor(4);

    cBA[ana] = new TCanvas(Form("cBA%i",ana),Form("TRDpTDep%s",Analyse[ana]));
    cBA[ana]->SetGrid();
    cBA[ana]->SetBottomMargin(0.16);
    cBA[ana]->SetRightMargin(0.03);
    cBA[ana]->SetLeftMargin(0.15);
    fpTDepGraph[5][ana][1][4]->Draw("AP");
    //legPIDCut[4]->Draw("same");
    //legRigCut[5][1]->Draw("same");
    //legPreselection[5][1]->Draw("same");
    legBIG->Draw("same");

    cBA[ana]->SaveAs(Form("~/TRDpTDep%s.%s",Analyse[ana],File[0]));
  }


}


Double_t GaussError(Double_t a, Double_t da, Double_t b, Double_t db) // for function f = a/b
{

  return (a/b) * TMath::Sqrt( (da/a)*(da/a) + (db/b)*(db/b) );

}
