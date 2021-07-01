/*--------------------------------------------------------------------------------------------------*\
|                                                                                                    |
|    Macro that makes projection histograms with branches of "nmartin_TestTriggerTRDTree.root"       |
|    - TRD PID vs. pT for all particles                                                              |
|      - with preselection of TPC, TOF, TPC && TOF for:<seperated particles> = Electrons,Pions,Kaons,|
|                                                      Protons, Deuterons, Tritons, Helium3, Alphas  |
|    - TRD PID Y-projection histograms for: all particles                                            |
|      - with preselection of TPC, TOF, TPC && TOF for:<seperated particles>                         |
|                                                                                                    |
|                                                                                                    |
| Author: Benjamin Brudnyj (2016)                                                                    | 
\*--------------------------------------------------------------------------------------------------*/
void bbrudnyj_ReadTreeTestTRDTriggerPID() {

  static const Int_t nFiles = 2;

  gStyle->SetTitleW(.8f);
  gStyle->SetTitleH(.08f);
  gStyle->SetTitleSize(.07,"xy");
  //gStyle->SetTitleOffset(.8,"x");
  gROOT->ForceStyle();

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Define Variables
  ///__________________
  Int_t Sigma      = 3;                   // 2 or 3 is testing in this analysis

  char *AllOrAnti  = "All";               // which particles shall be analyse?
  //char *AllOrAnti  = "Anti";

  //char *Collisions = "PbPb_11h_tracklets";
  char *Collisions = "PbPb_11h"; char *Tracklets = "";
  //char *Collisions = "PbPb_11hStd";
  //char *Collisions = "pPb_13b";
  //char *Collisions = "pPb_13c";
  //char *Collisions = "pp_15f";

  //char *Tracklets                   = "4trkl"; Double_t fnTrkl = 4; // exists only in PbPb_11h_tracklets
  //char *Tracklets                   = "5trkl"; Double_t fnTrkl = 5;
  //char *Tracklets                   = "6trkl"; Double_t fnTrkl = 6;

  //Double_t pTCut                    = 2.; char *pT2 = "pT2";
  Double_t pTCut                    = 1000.; char *pT2 = "";  // equivalent to no pT Cut

  char *File[nFiles]                = {"pdf","root"};

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

  printf(Form("\nCreating 6 canvases with Sigma = %i from %s\n",Sigma,Collisions));
  if(Tracklets != "" && pTCut != 2.) printf(Form("with %.f tracklet tracks\n\n",fnTrkl));
  else if(Tracklets != "" && pTCut == 2.) printf(Form("with %.f tracklet tracks and pT < %.f\n\n",fnTrkl,pTCut));
  else if(Tracklets == "" && pTCut == 2.) printf(Form("and pT < %.f\n\n",pTCut));

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Get data from file /lustre/nyx/alice/users/bbrudnyj/trunk/nmartin_nuclei/AliAnalysisTaskTestTriggerTRD.cxx
  /// Get histograms from /lustre/nyx/alice/users/bbrudnyj/codes/bbrudnyj_ReadTreeTestTRDTrigger2PID.C
  ///_________________________________________________________________________________________________
  TFile *cinput     = TFile::Open(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger2/bbrudnyj_TestTRDTrigger2%s%iSigma%sPID%s.root",Collisions,Tracklets,Collisions,Sigma,Tracklets,pT2));
  TList *inlistPID  = (TList*)cinput->Get(Form("bbrudnyj_TestTRDTrigger2%s%iSigma%sPID%s;1",Collisions,Sigma,Tracklets,pT2));
  TList *inlistProj = (TList*)cinput->Get(Form("bbrudnyj_TestTRDTrigger2%s%iSigma%sProj%s;1",Collisions,Sigma,Tracklets,pT2));


  if(AllOrAnti == "All") {
  // PID histograms
    TH2D *fHistTRDpTvPID                        = (TH2D*)inlistPID->FindObject("fHistTRDpTvPID");
    // TPC preselection
    TH2D *fHistTRDpTvPIDvElectronTPC            = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvElectronTPC");
    TH2D *fHistTRDpTvPIDvPionTPC                = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvPionTPC");
    TH2D *fHistTRDpTvPIDvKaonTPC                = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvKaonTPC");
    TH2D *fHistTRDpTvPIDvProtonTPC              = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvProtonTPC");
    TH2D *fHistTRDpTvPIDvDeuteronTPC            = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvDeuteronTPC");
    TH2D *fHistTRDpTvPIDvDeuteronTPC_clear      = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvDeuteronTPC_clear");
    TH2D *fHistTRDpTvPIDvTritonTPC              = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvTritonTPC");
    TH2D *fHistTRDpTvPIDvTritonTPC_clear        = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvTritonTPC_clear");
    TH2D *fHistTRDpTvPIDvHelium3TPC             = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvHelium3TPC");
    TH2D *fHistTRDpTvPIDvAlphaTPC               = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAlphaTPC");
    // TOF preselection
    TH2D *fHistTRDpTvPIDvElectronTOF            = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvElectronTOF");
    TH2D *fHistTRDpTvPIDvPionTOF                = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvPionTOF");
    TH2D *fHistTRDpTvPIDvKaonTOF                = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvKaonTOF");
    TH2D *fHistTRDpTvPIDvProtonTOF              = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvProtonTOF");
    TH2D *fHistTRDpTvPIDvDeuteronTOF            = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvDeuteronTOF");
    TH2D *fHistTRDpTvPIDvTritonTOF              = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvTritonTOF");
    TH2D *fHistTRDpTvPIDvHelium3TOF             = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvHelium3TOF");
    TH2D *fHistTRDpTvPIDvAlphaTOF               = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAlphaTOF");
    // TPC & TOF preselection
    TH2D *fHistTRDpTvPIDvElectronTPCTOF         = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvElectronTPCTOF");
    TH2D *fHistTRDpTvPIDvPionTPCTOF             = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvPionTPCTOF");
    TH2D *fHistTRDpTvPIDvKaonTPCTOF             = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvKaonTPCTOF");
    TH2D *fHistTRDpTvPIDvProtonTPCTOF           = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvProtonTPCTOF");
    TH2D *fHistTRDpTvPIDvDeuteronTPCTOF         = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvDeuteronTPCTOF");
    TH2D *fHistTRDpTvPIDvDeuteronTPCTOF_clear   = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvDeuteronTPCTOF_clear");
    TH2D *fHistTRDpTvPIDvTritonTPCTOF           = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvTritonTPCTOF");
    TH2D *fHistTRDpTvPIDvTritonTPCTOF_clear     = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvTritonTPCTOF_clear");
    TH2D *fHistTRDpTvPIDvHelium3TPCTOF          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvHelium3TPCTOF");
    TH2D *fHistTRDpTvPIDvAlphaTPCTOF            = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAlphaTPCTOF");

  // Projection histograms
    // TPC preselection
    TH1D *fHistTRDPIDallTPC                     = (TH1D*)inlistProj->FindObject("fHistTRDPIDallTPC");

    TH1D *fHistTRDPIDElectronTPC                = (TH1D*)inlistProj->FindObject("fHistTRDPIDElectronTPC");
    TH1D *fHistTRDPIDPionTPC                    = (TH1D*)inlistProj->FindObject("fHistTRDPIDPionTPC");
    TH1D *fHistTRDPIDKaonTPC                    = (TH1D*)inlistProj->FindObject("fHistTRDPIDKaonTPC");
    TH1D *fHistTRDPIDProtonTPC                  = (TH1D*)inlistProj->FindObject("fHistTRDPIDProtonTPC");
    TH1D *fHistTRDPIDDeuteronTPC                = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPC");
    TH1D *fHistTRDPIDDeuteronTPC_clear          = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPC_clear");
    TH1D *fHistTRDPIDTritonTPC                  = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPC");
    TH1D *fHistTRDPIDTritonTPC_clear            = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPC_clear");
    TH1D *fHistTRDPIDHelium3TPC                 = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3TPC");
    TH1D *fHistTRDPIDAlphaTPC                   = (TH1D*)inlistProj->FindObject("fHistTRDPIDAlphaTPC");

    TH1D *fHistTRDPIDElectronTPC_Clone          = (TH1D*)inlistProj->FindObject("fHistTRDPIDElectronTPC_Clone");
    TH1D *fHistTRDPIDPionTPC_Clone              = (TH1D*)inlistProj->FindObject("fHistTRDPIDPionTPC_Clone");
    TH1D *fHistTRDPIDKaonTPC_Clone              = (TH1D*)inlistProj->FindObject("fHistTRDPIDKaonTPC_Clone");
    TH1D *fHistTRDPIDProtonTPC_Clone            = (TH1D*)inlistProj->FindObject("fHistTRDPIDProtonTPC_Clone");
    TH1D *fHistTRDPIDDeuteronTPC_Clone          = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPC_Clone");
    TH1D *fHistTRDPIDDeuteronTPC_clear_Clone    = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPC_clear_Clone");
    TH1D *fHistTRDPIDTritonTPC_Clone            = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPC_Clone");
    TH1D *fHistTRDPIDTritonTPC_clear_Clone      = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPC_clear_Clone");
    TH1D *fHistTRDPIDHelium3TPC_Clone           = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3TPC_Clone");
    TH1D *fHistTRDPIDAlphaTPC_Clone             = (TH1D*)inlistProj->FindObject("fHistTRDPIDAlphaTPC_Clone");
    // TOF preselection
    TH1D *fHistTRDPIDallTOF                     = (TH1D*)inlistProj->FindObject("fHistTRDPIDallTOF");

    TH1D *fHistTRDPIDElectronTOF                = (TH1D*)inlistProj->FindObject("fHistTRDPIDElectronTOF");
    TH1D *fHistTRDPIDPionTOF                    = (TH1D*)inlistProj->FindObject("fHistTRDPIDPionTOF");
    TH1D *fHistTRDPIDKaonTOF                    = (TH1D*)inlistProj->FindObject("fHistTRDPIDKaonTOF");
    TH1D *fHistTRDPIDProtonTOF                  = (TH1D*)inlistProj->FindObject("fHistTRDPIDProtonTOF");
    TH1D *fHistTRDPIDDeuteronTOF                = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTOF");
    TH1D *fHistTRDPIDTritonTOF                  = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTOF");
    TH1D *fHistTRDPIDHelium3TOF                 = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3TOF");
    TH1D *fHistTRDPIDAlphaTOF                   = (TH1D*)inlistProj->FindObject("fHistTRDPIDAlphaTOF");

    TH1D *fHistTRDPIDElectronTOF_Clone          = (TH1D*)inlistProj->FindObject("fHistTRDPIDElectronTOF_Clone");
    TH1D *fHistTRDPIDPionTOF_Clone              = (TH1D*)inlistProj->FindObject("fHistTRDPIDPionTOF_Clone");
    TH1D *fHistTRDPIDKaonTOF_Clone              = (TH1D*)inlistProj->FindObject("fHistTRDPIDKaonTOF_Clone");
    TH1D *fHistTRDPIDProtonTOF_Clone            = (TH1D*)inlistProj->FindObject("fHistTRDPIDProtonTOF_Clone");
    TH1D *fHistTRDPIDDeuteronTOF_Clone          = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTOF_Clone");
    TH1D *fHistTRDPIDTritonTOF_Clone            = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTOF_Clone");
    TH1D *fHistTRDPIDHelium3TOF_Clone           = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3TOF_Clone");
    TH1D *fHistTRDPIDAlphaTOF_Clone             = (TH1D*)inlistProj->FindObject("fHistTRDPIDAlphaTOF_Clone");
    // TPC & TOF preselection
    TH1D *fHistTRDPIDallTPCTOF                  = (TH1D*)inlistProj->FindObject("fHistTRDPIDallTPCTOF");

    TH1D *fHistTRDPIDElectronTPCTOF             = (TH1D*)inlistProj->FindObject("fHistTRDPIDElectronTPCTOF");
    TH1D *fHistTRDPIDPionTPCTOF                 = (TH1D*)inlistProj->FindObject("fHistTRDPIDPionTPCTOF");
    TH1D *fHistTRDPIDKaonTPCTOF                 = (TH1D*)inlistProj->FindObject("fHistTRDPIDKaonTPCTOF");
    TH1D *fHistTRDPIDProtonTPCTOF               = (TH1D*)inlistProj->FindObject("fHistTRDPIDProtonTPCTOF");
    TH1D *fHistTRDPIDDeuteronTPCTOF             = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPCTOF");
    TH1D *fHistTRDPIDDeuteronTPCTOF_clear       = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPCTOF_clear");
    TH1D *fHistTRDPIDTritonTPCTOF               = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPCTOF");
    TH1D *fHistTRDPIDTritonTPCTOF_clear         = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPCTOF_clear");
    TH1D *fHistTRDPIDHelium3TPCTOF              = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3TPCTOF");
    TH1D *fHistTRDPIDAlphaTPCTOF                = (TH1D*)inlistProj->FindObject("fHistTRDPIDAlphaTPCTOF");

    TH1D *fHistTRDPIDElectronTPCTOF_Clone       = (TH1D*)inlistProj->FindObject("fHistTRDPIDElectronTPCTOF_Clone");
    TH1D *fHistTRDPIDPionTPCTOF_Clone           = (TH1D*)inlistProj->FindObject("fHistTRDPIDPionTPCTOF_Clone");
    TH1D *fHistTRDPIDKaonTPCTOF_Clone           = (TH1D*)inlistProj->FindObject("fHistTRDPIDKaonTPCTOF_Clone");
    TH1D *fHistTRDPIDProtonTPCTOF_Clone         = (TH1D*)inlistProj->FindObject("fHistTRDPIDProtonTPCTOF_Clone");
    TH1D *fHistTRDPIDDeuteronTPCTOF_Clone       = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPCTOF_Clone");
    TH1D *fHistTRDPIDDeuteronTPCTOF_clear_Clone = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPCTOF_clear_Clone");
    TH1D *fHistTRDPIDTritonTPCTOF_Clone         = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPCTOF_Clone");
    TH1D *fHistTRDPIDTritonTPCTOF_clear_Clone   = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPCTOF_clear_Clone");
    TH1D *fHistTRDPIDHelium3TPCTOF_Clone        = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3TPCTOF_Clone");
    TH1D *fHistTRDPIDAlphaTPCTOF_Clone          = (TH1D*)inlistProj->FindObject("fHistTRDPIDAlphaTPCTOF_Clone");
  }



  if(AllOrAnti == "Anti") {
  // PID histograms
    TH2D *fHistTRDpTvPID                        = (TH2D*)inlistPID->FindObject("fHistTRDpTvPID");
    // TPC preselection
    TH2D *fHistTRDpTvPIDvElectronTPC            = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvElectronTPC");
    TH2D *fHistTRDpTvPIDvPionTPC                = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvPionTPC");
    TH2D *fHistTRDpTvPIDvKaonTPC                = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvKaonTPC");
    TH2D *fHistTRDpTvPIDvProtonTPC              = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvProtonTPC");
    TH2D *fHistTRDpTvPIDvDeuteronTPC            = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvDeuteronTPC");
    TH2D *fHistTRDpTvPIDvDeuteronTPC_clear      = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiDeuteronTPC_clear");
    TH2D *fHistTRDpTvPIDvTritonTPC              = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvTritonTPC");
    TH2D *fHistTRDpTvPIDvTritonTPC_clear        = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiTritonTPC_clear");
    TH2D *fHistTRDpTvPIDvHelium3TPC             = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiHelium3TPC");
    TH2D *fHistTRDpTvPIDvAlphaTPC               = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiAlphaTPC");
    // TOF preselection
    TH2D *fHistTRDpTvPIDvElectronTOF            = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvElectronTOF");
    TH2D *fHistTRDpTvPIDvPionTOF                = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvPionTOF");
    TH2D *fHistTRDpTvPIDvKaonTOF                = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvKaonTOF");
    TH2D *fHistTRDpTvPIDvProtonTOF              = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvProtonTOF");
    TH2D *fHistTRDpTvPIDvDeuteronTOF            = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvDeuteronTOF");
    TH2D *fHistTRDpTvPIDvTritonTOF              = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvTritonTOF");
    TH2D *fHistTRDpTvPIDvHelium3TOF             = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvHelium3TOF");
    TH2D *fHistTRDpTvPIDvAlphaTOF               = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAlphaTOF");
    // TPC & TOF preselection
    TH2D *fHistTRDpTvPIDvElectronTPCTOF         = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvElectronTPCTOF");
    TH2D *fHistTRDpTvPIDvPionTPCTOF             = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvPionTPCTOF");
    TH2D *fHistTRDpTvPIDvKaonTPCTOF             = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvKaonTPCTOF");
    TH2D *fHistTRDpTvPIDvProtonTPCTOF           = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvProtonTPCTOF");
    TH2D *fHistTRDpTvPIDvDeuteronTPCTOF         = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvDeuteronTPCTOF");
    TH2D *fHistTRDpTvPIDvDeuteronTPCTOF_clear   = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiDeuteronTPCTOF_clear");
    TH2D *fHistTRDpTvPIDvTritonTPCTOF           = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvTritonTPCTOF");
    TH2D *fHistTRDpTvPIDvTritonTPCTOF_clear     = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiTritonTPCTOF_clear");
    TH2D *fHistTRDpTvPIDvHelium3TPCTOF          = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiHelium3TPCTOF");
    TH2D *fHistTRDpTvPIDvAlphaTPCTOF            = (TH2D*)inlistPID->FindObject("fHistTRDpTvPIDvAntiAlphaTPCTOF");

  // Projection histograms
    // TPC preselection
    TH1D *fHistTRDPIDallTPC                     = (TH1D*)inlistProj->FindObject("fHistTRDPIDallTPC");

    TH1D *fHistTRDPIDElectronTPC                = (TH1D*)inlistProj->FindObject("fHistTRDPIDElectronTPC");
    TH1D *fHistTRDPIDPionTPC                    = (TH1D*)inlistProj->FindObject("fHistTRDPIDPionTPC");
    TH1D *fHistTRDPIDKaonTPC                    = (TH1D*)inlistProj->FindObject("fHistTRDPIDKaonTPC");
    TH1D *fHistTRDPIDProtonTPC                  = (TH1D*)inlistProj->FindObject("fHistTRDPIDProtonTPC");
    TH1D *fHistTRDPIDDeuteronTPC                = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPC");
    TH1D *fHistTRDPIDDeuteronTPC_clear          = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuteronTPC_clear");
    TH1D *fHistTRDPIDTritonTPC                  = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPC");
    TH1D *fHistTRDPIDTritonTPC_clear            = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiTritonTPC_clear");
    TH1D *fHistTRDPIDHelium3TPC                 = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3TPC");
    TH1D *fHistTRDPIDAlphaTPC                   = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiAlphaTPC");

    TH1D *fHistTRDPIDElectronTPC_Clone          = (TH1D*)inlistProj->FindObject("fHistTRDPIDElectronTPC_Clone");
    TH1D *fHistTRDPIDPionTPC_Clone              = (TH1D*)inlistProj->FindObject("fHistTRDPIDPionTPC_Clone");
    TH1D *fHistTRDPIDKaonTPC_Clone              = (TH1D*)inlistProj->FindObject("fHistTRDPIDKaonTPC_Clone");
    TH1D *fHistTRDPIDProtonTPC_Clone            = (TH1D*)inlistProj->FindObject("fHistTRDPIDProtonTPC_Clone");
    TH1D *fHistTRDPIDDeuteronTPC_Clone          = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPC_Clone");
    TH1D *fHistTRDPIDDeuteronTPC_clear_Clone    = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuteronTPC_clear_Clone");
    TH1D *fHistTRDPIDTritonTPC_Clone            = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPC_Clone");
    TH1D *fHistTRDPIDTritonTPC_clear_Clone      = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiTritonTPC_clear_Clone");
    TH1D *fHistTRDPIDHelium3TPC_Clone           = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3TPC_Clone");
    TH1D *fHistTRDPIDAlphaTPC_Clone             = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiAlphaTPC_Clone");
    // TOF preselection
    TH1D *fHistTRDPIDallTOF                     = (TH1D*)inlistProj->FindObject("fHistTRDPIDallTOF");

    TH1D *fHistTRDPIDElectronTOF                = (TH1D*)inlistProj->FindObject("fHistTRDPIDElectronTOF");
    TH1D *fHistTRDPIDPionTOF                    = (TH1D*)inlistProj->FindObject("fHistTRDPIDPionTOF");
    TH1D *fHistTRDPIDKaonTOF                    = (TH1D*)inlistProj->FindObject("fHistTRDPIDKaonTOF");
    TH1D *fHistTRDPIDProtonTOF                  = (TH1D*)inlistProj->FindObject("fHistTRDPIDProtonTOF");
    TH1D *fHistTRDPIDDeuteronTOF                = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTOF");
    TH1D *fHistTRDPIDTritonTOF                  = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTOF");
    TH1D *fHistTRDPIDHelium3TOF                 = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3TOF");
    TH1D *fHistTRDPIDAlphaTOF                   = (TH1D*)inlistProj->FindObject("fHistTRDPIDAlphaTOF");

    TH1D *fHistTRDPIDElectronTOF_Clone          = (TH1D*)inlistProj->FindObject("fHistTRDPIDElectronTOF_Clone");
    TH1D *fHistTRDPIDPionTOF_Clone              = (TH1D*)inlistProj->FindObject("fHistTRDPIDPionTOF_Clone");
    TH1D *fHistTRDPIDKaonTOF_Clone              = (TH1D*)inlistProj->FindObject("fHistTRDPIDKaonTOF_Clone");
    TH1D *fHistTRDPIDProtonTOF_Clone            = (TH1D*)inlistProj->FindObject("fHistTRDPIDProtonTOF_Clone");
    TH1D *fHistTRDPIDDeuteronTOF_Clone          = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTOF_Clone");
    TH1D *fHistTRDPIDTritonTOF_Clone            = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTOF_Clone");
    TH1D *fHistTRDPIDHelium3TOF_Clone           = (TH1D*)inlistProj->FindObject("fHistTRDPIDHelium3TOF_Clone");
    TH1D *fHistTRDPIDAlphaTOF_Clone             = (TH1D*)inlistProj->FindObject("fHistTRDPIDAlphaTOF_Clone");
    // TPC & TOF preselection
    TH1D *fHistTRDPIDallTPCTOF                  = (TH1D*)inlistProj->FindObject("fHistTRDPIDallTPCTOF");

    TH1D *fHistTRDPIDElectronTPCTOF             = (TH1D*)inlistProj->FindObject("fHistTRDPIDElectronTPCTOF");
    TH1D *fHistTRDPIDPionTPCTOF                 = (TH1D*)inlistProj->FindObject("fHistTRDPIDPionTPCTOF");
    TH1D *fHistTRDPIDKaonTPCTOF                 = (TH1D*)inlistProj->FindObject("fHistTRDPIDKaonTPCTOF");
    TH1D *fHistTRDPIDProtonTPCTOF               = (TH1D*)inlistProj->FindObject("fHistTRDPIDProtonTPCTOF");
    TH1D *fHistTRDPIDDeuteronTPCTOF             = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPCTOF");
    TH1D *fHistTRDPIDDeuteronTPCTOF_clear       = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuteronTPCTOF_clear");
    TH1D *fHistTRDPIDTritonTPCTOF               = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPCTOF");
    TH1D *fHistTRDPIDTritonTPCTOF_clear         = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiTritonTPCTOF_clear");
    TH1D *fHistTRDPIDHelium3TPCTOF              = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3TPCTOF");
    TH1D *fHistTRDPIDAlphaTPCTOF                = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiAlphaTPCTOF");

    TH1D *fHistTRDPIDElectronTPCTOF_Clone       = (TH1D*)inlistProj->FindObject("fHistTRDPIDElectronTPCTOF_Clone");
    TH1D *fHistTRDPIDPionTPCTOF_Clone           = (TH1D*)inlistProj->FindObject("fHistTRDPIDPionTPCTOF_Clone");
    TH1D *fHistTRDPIDKaonTPCTOF_Clone           = (TH1D*)inlistProj->FindObject("fHistTRDPIDKaonTPCTOF_Clone");
    TH1D *fHistTRDPIDProtonTPCTOF_Clone         = (TH1D*)inlistProj->FindObject("fHistTRDPIDProtonTPCTOF_Clone");
    TH1D *fHistTRDPIDDeuteronTPCTOF_Clone       = (TH1D*)inlistProj->FindObject("fHistTRDPIDDeuteronTPCTOF_Clone");
    TH1D *fHistTRDPIDDeuteronTPCTOF_clear_Clone = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiDeuteronTPCTOF_clear_Clone");
    TH1D *fHistTRDPIDTritonTPCTOF_Clone         = (TH1D*)inlistProj->FindObject("fHistTRDPIDTritonTPCTOF_Clone");
    TH1D *fHistTRDPIDTritonTPCTOF_clear_Clone   = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiTritonTPCTOF_clear_Clone");
    TH1D *fHistTRDPIDHelium3TPCTOF_Clone        = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiHelium3TPCTOF_Clone");
    TH1D *fHistTRDPIDAlphaTPCTOF_Clone          = (TH1D*)inlistProj->FindObject("fHistTRDPIDAntiAlphaTPCTOF_Clone");
  }

  fHistTRDpTvPID->SetTitle("");

  fHistTRDPIDallTPC->SetStats(0);
  fHistTRDPIDallTPCTOF->SetStats(0);

  fHistTRDpTvPID->SetTitle("");
  fHistTRDpTvPID->GetXaxis()->SetTitle("#it{p}_{T}/#it{Z} (GeV/#it{c})");
  fHistTRDpTvPID->GetYaxis()->SetTitle("TRD  <#it{Q}_{s}> (arb. units)");
  fHistTRDpTvPID->SetStats(0);

  fHistTRDpTvPIDvHelium3TPC->SetTitle("");
  fHistTRDpTvPIDvHelium3TPC->GetXaxis()->SetTitle("#it{p}_{T}/#it{Z} (GeV/#it{c})");
  fHistTRDpTvPIDvHelium3TPC->GetYaxis()->SetTitle("TRD  <#it{Q}_{s}> (arb. units)");
  fHistTRDpTvPIDvHelium3TPC->SetStats(0);
  fHistTRDpTvPIDvDeuteronTPC_clear->SetTitle("");
  fHistTRDpTvPIDvDeuteronTPC_clear->GetXaxis()->SetTitle("#it{p}_{T}/#it{Z} (GeV/#it{c})");
  fHistTRDpTvPIDvDeuteronTPC_clear->GetYaxis()->SetTitle("TRD  <#it{Q}_{s}> (arb. units)");
  fHistTRDpTvPIDvDeuteronTPC_clear->SetStats(0);
  fHistTRDpTvPIDvTritonTPC_clear->SetTitle("");
  fHistTRDpTvPIDvTritonTPC_clear->GetXaxis()->SetTitle("#it{p}_{T}/#it{Z} (GeV/#it{c})");
  fHistTRDpTvPIDvTritonTPC_clear->GetYaxis()->SetTitle("TRD  <#it{Q}_{s}> (arb. units)");
  fHistTRDpTvPIDvTritonTPC_clear->SetStats(0);
  fHistTRDpTvPIDvAlphaTPC->SetTitle("");
  fHistTRDpTvPIDvAlphaTPC->GetXaxis()->SetTitle("#it{p}_{T}/#it{Z} (GeV/#it{c})");
  fHistTRDpTvPIDvAlphaTPC->GetYaxis()->SetTitle("TRD  <#it{Q}_{s}> (arb. units)");
  fHistTRDpTvPIDvAlphaTPC->SetStats(0);

  fHistTRDpTvPIDvHelium3TPCTOF->SetTitle("");
  fHistTRDpTvPIDvHelium3TPCTOF->GetXaxis()->SetTitle("#it{p}_{T}/#it{Z} (GeV/#it{c})");
  fHistTRDpTvPIDvHelium3TPCTOF->GetYaxis()->SetTitle("TRD  <#it{Q}_{s}> (arb. units)");
  fHistTRDpTvPIDvHelium3TPCTOF->SetStats(0);
  fHistTRDpTvPIDvHelium3TPCTOF->GetZaxis()->SetRangeUser(0.5,800);
  fHistTRDpTvPIDvDeuteronTPCTOF_clear->SetTitle("");
  fHistTRDpTvPIDvDeuteronTPCTOF_clear->GetXaxis()->SetTitle("#it{p}_{T}/#it{Z} (GeV/#it{c})");
  fHistTRDpTvPIDvDeuteronTPCTOF_clear->GetYaxis()->SetTitle("TRD  <#it{Q}_{s}> (arb. units)");
  fHistTRDpTvPIDvDeuteronTPCTOF_clear->SetStats(0);
  fHistTRDpTvPIDvDeuteronTPCTOF_clear->GetZaxis()->SetRangeUser(0.5,800);
  fHistTRDpTvPIDvTritonTPCTOF_clear->SetTitle("");
  fHistTRDpTvPIDvTritonTPCTOF_clear->GetXaxis()->SetTitle("#it{p}_{T}/#it{Z} (GeV/#it{c})");
  fHistTRDpTvPIDvTritonTPCTOF_clear->GetYaxis()->SetTitle("TRD  <#it{Q}_{s}> (arb. units)");
  fHistTRDpTvPIDvTritonTPCTOF_clear->SetStats(0);
  fHistTRDpTvPIDvTritonTPCTOF_clear->GetZaxis()->SetRangeUser(0.5,800);
  fHistTRDpTvPIDvAlphaTPCTOF->SetTitle("");
  fHistTRDpTvPIDvAlphaTPCTOF->GetXaxis()->SetTitle("#it{p}_{T}/#it{Z} (GeV/#it{c})");
  fHistTRDpTvPIDvAlphaTPCTOF->GetYaxis()->SetTitle("TRD  <#it{Q}_{s}> (arb. units)");
  fHistTRDpTvPIDvAlphaTPCTOF->SetStats(0);
  fHistTRDpTvPIDvAlphaTPCTOF->GetZaxis()->SetRangeUser(0.5,800);


  fHistTRDPIDallTPC->GetXaxis()->SetTitle("TRD  <#it{Q}_{s}> (arb. units)");
  fHistTRDPIDallTPCTOF->GetXaxis()->SetTitle("TRD  <#it{Q}_{s}> (arb. units)");

  fHistTRDPIDTritonTPC_clear->SetLineColor(6);
  fHistTRDPIDTritonTPCTOF_clear->SetLineColor(6);

  fHistTRDPIDKaonTPC          ->SetLineWidth(2);
  fHistTRDPIDAlphaTPC         ->SetLineWidth(2);
  fHistTRDPIDElectronTPC      ->SetLineWidth(2);
  fHistTRDPIDProtonTPC        ->SetLineWidth(2);
  fHistTRDPIDPionTPC          ->SetLineWidth(2);
  fHistTRDPIDHelium3TPC       ->SetLineWidth(2);
  fHistTRDPIDDeuteronTPC_clear->SetLineWidth(2);
  fHistTRDPIDDeuteronTPC_clear->SetLineColor(9);
  fHistTRDPIDTritonTPC_clear  ->SetLineWidth(2);
  fHistTRDPIDTritonTPC_clear  ->SetLineColor(8);

  fHistTRDPIDElectronTPCTOF      ->SetLineWidth(2);
  fHistTRDPIDProtonTPCTOF        ->SetLineWidth(2);
  fHistTRDPIDKaonTPCTOF          ->SetLineWidth(2);
  fHistTRDPIDPionTPCTOF          ->SetLineWidth(2);
  fHistTRDPIDDeuteronTPCTOF_clear ->SetStats(0);
  fHistTRDPIDDeuteronTPCTOF_clear->SetLineWidth(2);
  fHistTRDPIDDeuteronTPCTOF_clear->SetLineColor(4);
  fHistTRDPIDDeuteronTPCTOF_clear->SetTitle("");
  fHistTRDPIDDeuteronTPCTOF_clear->GetXaxis()->SetTitle("TRD  <#it{Q}_{s}> (arb. units)");
  fHistTRDPIDTritonTPCTOF_clear->SetStats(0);
  fHistTRDPIDTritonTPCTOF_clear->SetLineWidth(2);
  fHistTRDPIDTritonTPCTOF_clear->SetLineColor(4);
  fHistTRDPIDTritonTPCTOF_clear->SetTitle("");
  fHistTRDPIDTritonTPCTOF_clear->GetXaxis()->SetTitle("TRD  <#it{Q}_{s}> (arb. units)");
  fHistTRDPIDHelium3TPCTOF->SetStats(0);
  fHistTRDPIDHelium3TPCTOF->SetLineWidth(2);
  fHistTRDPIDHelium3TPCTOF->SetLineColor(4);
  fHistTRDPIDHelium3TPCTOF->SetTitle("");
  fHistTRDPIDHelium3TPCTOF->GetXaxis()->SetTitle("TRD  <#it{Q}_{s}> (arb. units)");
  fHistTRDPIDAlphaTPCTOF->SetStats(0);
  fHistTRDPIDAlphaTPCTOF->SetLineWidth(2);
  fHistTRDPIDAlphaTPCTOF->SetLineColor(4);
  fHistTRDPIDAlphaTPCTOF->SetTitle("");
  fHistTRDPIDAlphaTPCTOF->GetXaxis()->SetTitle("TRD  <#it{Q}_{s}> (arb. units)");



  fHistTRDPIDHelium3TPC->SetStats(0);
  fHistTRDPIDHelium3TPC->SetLineWidth(2);
  fHistTRDPIDHelium3TPC->SetLineColor(1);
  fHistTRDPIDHelium3TPC->SetTitle("");
  fHistTRDPIDHelium3TPC->GetXaxis()->SetTitle("TRD  <#it{Q}_{s}> (arb. units)");
  fHistTRDPIDAlphaTPC->SetStats(0);
  fHistTRDPIDAlphaTPC->SetLineWidth(2);
  fHistTRDPIDAlphaTPC->SetLineColor(16);
  fHistTRDPIDAlphaTPC->SetTitle("");
  fHistTRDPIDAlphaTPC->GetXaxis()->SetTitle("TRD  <#it{Q}_{s}> (arb. units)");

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Create legends
  ///_________________________________________________
  TLegend *legTracklets = new TLegend(.57,.52,.95,.6);
  legTracklets->SetFillStyle(0);
  legTracklets->SetBorderSize(0);
  if(Tracklets != "") legTracklets->AddEntry((TObject*)0,Form("%.f tracklet tracks",fnTrkl),"");
  TLegend *legTracklets2 = new TLegend(.44,.47,.82,.55);
  legTracklets2->SetFillStyle(0);
  legTracklets2->SetBorderSize(0);
  if(Tracklets != "") legTracklets2->AddEntry((TObject*)0,Form("%.f tracklet tracks",fnTrkl),"");
  TLegend *legTracklets3 = new TLegend(.44,.43,.82,.51);
  legTracklets3->SetFillStyle(0);
  legTracklets3->SetBorderSize(0);
  if(Tracklets != "") legTracklets3->AddEntry((TObject*)0,Form("%.f tracklet tracks",fnTrkl),"");
  TLegend *legTracklets4 = new TLegend(.12,.47,.5,.55);
  legTracklets4->SetFillStyle(0);
  legTracklets4->SetBorderSize(0);
  if(Tracklets != "") legTracklets4->AddEntry((TObject*)0,Form("%.f tracklet tracks",fnTrkl),"");

  TLegend *legpTCut = new TLegend(.57,.6,.95,.66);
  legpTCut->SetFillStyle(0);
  legpTCut->SetBorderSize(0);
  if(pTCut == 2.) legpTCut->AddEntry((TObject*)0,Form("#frac{#it{p}_{T}}{#it{Z}} < %.f GeV/#it{c}",pTCut),"");
  TLegend *legpTCut2 = new TLegend(.44,.55,.82,.61);
  legpTCut2->SetFillStyle(0);
  legpTCut2->SetBorderSize(0);
  if(pTCut == 2.) legpTCut2->AddEntry((TObject*)0,Form("#frac{#it{p}_{T}}{#it{Z}} < %.f GeV/#it{c}",pTCut),"");
  TLegend *legpTCut3 = new TLegend(.44,.51,.82,.57);
  legpTCut3->SetFillStyle(0);
  legpTCut3->SetBorderSize(0);
  if(pTCut == 2.) legpTCut3->AddEntry((TObject*)0,Form("#frac{#it{p}_{T}}{#it{Z}} < %.f GeV/#it{c}",pTCut),"");
  TLegend *legpTCut4 = new TLegend(.12,.55,.5,.61);
  legpTCut4->SetFillStyle(0);
  legpTCut4->SetBorderSize(0);
  if(pTCut == 2.) legpTCut4->AddEntry((TObject*)0,Form("#frac{#it{p}_{T}}{#it{Z}} < %.f GeV/#it{c}",pTCut),"");

  // Particles
  TLegend *legDeu = new TLegend(.05,.73,.4,.85);
  legDeu->SetFillStyle(0);
  legDeu->SetBorderSize(0);
  //legDeu->AddEntry((TObject*)0,"^{2}H","");
  legDeu->AddEntry((TObject*)0,"d","");
  TLegend *legTri = new TLegend(.05,.73,.4,.85);
  legTri->SetFillStyle(0);
  legTri->SetBorderSize(0);
  //legTri->AddEntry((TObject*)0,"^{3}H","");
  legTri->AddEntry((TObject*)0,"t","");
  TLegend *legHe3 = new TLegend(.05,.73,.4,.85);
  legHe3->SetFillStyle(0);
  legHe3->SetBorderSize(0);
  legHe3->AddEntry((TObject*)0,"^{3}He","");
  TLegend *legAlp = new TLegend(.05,.73,.4,.85);
  //TLegend *legAlp = new TLegend(.4,.73,.75,.85);
  legAlp->SetFillStyle(0);
  legAlp->SetBorderSize(0);
  legAlp->AddEntry((TObject*)0,"^{4}He","");


  // PID
  TLegend *legPID;
  if((pTCut == 2. && Tracklets == "") || (pTCut != 2. && Tracklets != "")) legPID = new TLegend(.05,.63,.42,.85);
  else if (pTCut == 2. && Tracklets != "") legPID = new TLegend(.05,.55,.42,.85);
  else                                     legPID = new TLegend(.05,.7,.42,.85);
  legPID->SetFillStyle(0);
  legPID->SetBorderSize(0);
  if(Collisions == "PbPb_11h" || Collisions == "PbPb_11hStd" || Collisions == "PbPb_11h_tracklets") {
    legPID->AddEntry((TObject*)0,"Pb-Pb, #sqrt{#it{s}_{NN}} = 2.76 TeV","");
    //legPID->AddEntry((TObject*)0,"run LHC11h pass2","");
    legPID->AddEntry((TObject*)0," ","");
  }
  else if(Collisions == "pPb_13b") {
    legPID->AddEntry((TObject*)0,"p-Pb, #sqrt{#it{s}_{NN}} = 5.023 TeV","");
    legPID->AddEntry((TObject*)0,"run LHC13b pass3","");
  }
  else if(Collisions == "pPb_13c") {
    legPID->AddEntry((TObject*)0,"p-Pb, #sqrt{#it{s}_{NN}} = 5.023 TeV","");
    legPID->AddEntry((TObject*)0,"run LHC13c pass2","");
  }
  else { //else if(Collisions == "pp_15f") {
    legPID->AddEntry((TObject*)0,"p-p, #sqrt{#it{s}} = 13 TeV","");
    legPID->AddEntry((TObject*)0,"run LHC15f pass2","");
  }
  if(pTCut == 2. && Tracklets == "")       legPID->AddEntry((TObject*)0,"#frac{#it{p}_{T}}{#it{Z}} < 2 GeV/#it{c}","");
  else if (pTCut != 2. && Tracklets != "") legPID->AddEntry((TObject*)0,Form("%.f tracklet tracks",fnTrkl),"");
  else if (pTCut == 2. && Tracklets != "") {
    legPID->AddEntry((TObject*)0,"#frac{#it{p}_{T}}{#it{Z}} < 2 GeV/#it{c}","");
    legPID->AddEntry((TObject*)0,Form("%.f tracklet tracks",fnTrkl),"");
  }

  TLegend *legPIDTPCnSigma = new TLegend(.05,.61,.4,.73);
  //TLegend *legPIDTPCnSigma = new TLegend(.4,.61,.75,.73);
  legPIDTPCnSigma->SetFillStyle(0);
  legPIDTPCnSigma->SetBorderSize(0);
  legPIDTPCnSigma->AddEntry((TObject*)0,"with preselection:","");
  legPIDTPCnSigma->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  TLegend *legPIDTPCnSigmaDeuteron = new TLegend(.05,.54,.4,.73);
  legPIDTPCnSigmaDeuteron->SetFillStyle(0);
  legPIDTPCnSigmaDeuteron->SetBorderSize(0);
  legPIDTPCnSigmaDeuteron->AddEntry((TObject*)0,"with preselection:","");
  legPIDTPCnSigmaDeuteron->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  legPIDTPCnSigmaDeuteron->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutDeuteronTPC),"");
  TLegend *legPIDTPCnSigmaTriton = new TLegend(.05,.54,.4,.73);
  legPIDTPCnSigmaTriton->SetFillStyle(0);
  legPIDTPCnSigmaTriton->SetBorderSize(0);
  legPIDTPCnSigmaTriton->AddEntry((TObject*)0,"with preselection:","");
  legPIDTPCnSigmaTriton->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  legPIDTPCnSigmaTriton->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutTritonTPC),"");
  TLegend *legPIDTOFnSigma = new TLegend(.05,.73,.4,.85);
  legPIDTOFnSigma->SetFillStyle(0);
  legPIDTOFnSigma->SetBorderSize(0);
  legPIDTOFnSigma->AddEntry((TObject*)0,"with preselection:","");
  legPIDTOFnSigma->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");
  TLegend *legPIDTPCTOFnSigma = new TLegend(.05,.55,.4,.73);
  legPIDTPCTOFnSigma->SetFillStyle(0);
  legPIDTPCTOFnSigma->SetBorderSize(0);
  legPIDTPCTOFnSigma->AddEntry((TObject*)0,"with preselection:","");
  legPIDTPCTOFnSigma->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  legPIDTPCTOFnSigma->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");
  TLegend *legPIDTPCTOFnSigmaDeuteron = new TLegend(.05,.48,.4,.73);
  legPIDTPCTOFnSigmaDeuteron->SetFillStyle(0);
  legPIDTPCTOFnSigmaDeuteron->SetBorderSize(0);
  legPIDTPCTOFnSigmaDeuteron->AddEntry((TObject*)0,"with preselection:","");
  legPIDTPCTOFnSigmaDeuteron->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  legPIDTPCTOFnSigmaDeuteron->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");
  legPIDTPCTOFnSigmaDeuteron->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutDeuteronTPCTOF),"");
  TLegend *legPIDTPCTOFnSigmaDeuteronright = new TLegend(.5,.6,.85,.85);
  legPIDTPCTOFnSigmaDeuteronright->SetFillStyle(0);
  legPIDTPCTOFnSigmaDeuteronright->SetBorderSize(0);
  legPIDTPCTOFnSigmaDeuteronright->AddEntry((TObject*)0,"with preselection:","");
  legPIDTPCTOFnSigmaDeuteronright->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  legPIDTPCTOFnSigmaDeuteronright->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");
  legPIDTPCTOFnSigmaDeuteronright->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutDeuteronTPCTOF),"");
  TLegend *legPIDTPCTOFnSigmaTriton = new TLegend(.05,.48,.4,.73);
  legPIDTPCTOFnSigmaTriton->SetFillStyle(0);
  legPIDTPCTOFnSigmaTriton->SetBorderSize(0);
  legPIDTPCTOFnSigmaTriton->AddEntry((TObject*)0,"with preselection:","");
  legPIDTPCTOFnSigmaTriton->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  legPIDTPCTOFnSigmaTriton->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");
  legPIDTPCTOFnSigmaTriton->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutTritonTPCTOF),"");
  TLegend *legPIDTPCTOFnSigmaTritonright = new TLegend(.5,.6,.85,.85);
  legPIDTPCTOFnSigmaTritonright->SetFillStyle(0);
  legPIDTPCTOFnSigmaTritonright->SetBorderSize(0);
  legPIDTPCTOFnSigmaTritonright->AddEntry((TObject*)0,"with preselection:","");
  legPIDTPCTOFnSigmaTritonright->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  legPIDTPCTOFnSigmaTritonright->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");
  legPIDTPCTOFnSigmaTritonright->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutTritonTPCTOF),"");

  // Proj legends
  TLegend *legProj;
  if((pTCut == 2. && Tracklets == "") || (pTCut != 2. && Tracklets != "")) legProj = new TLegend(.2,.63,.57,.85);
  else if (pTCut == 2. && Tracklets != "") legProj = new TLegend(.2,.55,.57,.85);
  //else                                     legProj = new TLegend(.2,.7,.57,.85);
  else                                     legProj = new TLegend(.4,.7,.77,.85);
  legProj->SetFillStyle(0);
  legProj->SetBorderSize(0);
  if(Collisions == "PbPb_11h" || Collisions == "PbPb_11hStd" || Collisions == "PbPb_11h_tracklets") {
    legProj->AddEntry((TObject*)0,"Pb-Pb, #sqrt{#it{s}_{NN}} = 2.76 TeV","");
    //legProj->AddEntry((TObject*)0,"run LHC11h pass2","");
    legProj->AddEntry((TObject*)0," ","");
  }
  else if(Collisions == "pPb_13b") {
    legProj->AddEntry((TObject*)0,"p-Pb, #sqrt{#it{s}_{NN}} = 5.023 TeV","");
    legProj->AddEntry((TObject*)0,"run LHC13b pass3","");
  }
  else if(Collisions == "pPb_13c") {
    legProj->AddEntry((TObject*)0,"p-Pb, #sqrt{#it{s}_{NN}} = 5.023 TeV","");
    legProj->AddEntry((TObject*)0,"run LHC13c pass2","");
  }
  else { //(Collisions == "pp_15f")
    legProj->AddEntry((TObject*)0,"p-p, #sqrt{#it{s}}} = 13 TeV","");
    legProj->AddEntry((TObject*)0,"run LHC15f pass2","");
  }
  if(pTCut == 2. && Tracklets == "")      legProj->AddEntry((TObject*)0,"#frac{#it{p}_{T}}{#it{Z}} < 2 GeV/#it{c}","");
  else if(pTCut != 2. && Tracklets != "") legProj->AddEntry((TObject*)0,Form("%.f tracklet tracks",fnTrkl),"");
  else if(pTCut == 2. && Tracklets != "") {
    legProj->AddEntry((TObject*)0,"#frac{#it{p}_{T}}{#it{Z}} < 2 GeV/#it{c}","");
    legProj->AddEntry((TObject*)0,Form("%.f tracklet tracks",fnTrkl),"");
  }

  // TPC legends
  TLegend *legTPC = new TLegend(.6,.45,.75,.85);
  legTPC->SetFillStyle(0);
  legTPC->SetBorderSize(0);
  legTPC->AddEntry(fHistTRDPIDElectronTPC,"Electron","l");
  legTPC->AddEntry(fHistTRDPIDPionTPC,"Pion","l");
  legTPC->AddEntry(fHistTRDPIDKaonTPC,"Kaon","l");
  legTPC->AddEntry(fHistTRDPIDProtonTPC,"Proton","l");
  legTPC->AddEntry(fHistTRDPIDHelium3TPC,"Helium3","l");
  legTPC->AddEntry(fHistTRDPIDDeuteronTPC_clear,"Deuteron","l");
  legTPC->AddEntry(fHistTRDPIDTritonTPC_clear,"Triton","l");
  legTPC->AddEntry(fHistTRDPIDAlphaTPC,"Alpha","l");
  TLegend *legTPC1 = new TLegend(.5,.51,1,.78);
  legTPC1->SetFillStyle(0);
  legTPC1->SetBorderSize(0);
  legTPC1->AddEntry(fHistTRDPIDDeuteronTPC_clear," d","l");
  legTPC1->AddEntry(fHistTRDPIDTritonTPC_clear," t","l");
  legTPC1->AddEntry(fHistTRDPIDallTPC,"all other tracks","l");
  TLegend *legTPC2 = new TLegend(.69,.6,1.09,.78);
  legTPC2->SetFillStyle(0);
  legTPC2->SetBorderSize(0);
  legTPC2->AddEntry(fHistTRDPIDHelium3TPC,"^{3}He","l");
  legTPC2->AddEntry(fHistTRDPIDAlphaTPC,"^{4}He","l");
  TLegend *legElectronTPC = new TLegend(.45,.67,.8,.85);
  legElectronTPC->SetFillStyle(0);
  legElectronTPC->SetBorderSize(0);
  legElectronTPC->AddEntry(fHistTRDPIDElectronTPC_Clone,"Electron","l");
  legElectronTPC->AddEntry((TObject*)0,"with preselection:","");
  legElectronTPC->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  TLegend *legPionTPC = new TLegend(.45,.67,.8,.85);
  legPionTPC->SetFillStyle(0);
  legPionTPC->SetBorderSize(0);
  legPionTPC->AddEntry(fHistTRDPIDPionTPC_Clone,"Pion","l");
  legPionTPC->AddEntry((TObject*)0,"with preselection:","");
  legPionTPC->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  TLegend *legKaonTPC = new TLegend(.45,.67,.8,.85);
  legKaonTPC->SetFillStyle(0);
  legKaonTPC->SetBorderSize(0);
  legKaonTPC->AddEntry(fHistTRDPIDKaonTPC_Clone,"Kaon","l");
  legKaonTPC->AddEntry((TObject*)0,"with preselection:","");
  legKaonTPC->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  TLegend *legProtonTPC = new TLegend(.45,.67,.8,.85);
  legProtonTPC->SetFillStyle(0);
  legProtonTPC->SetBorderSize(0);
  legProtonTPC->AddEntry(fHistTRDPIDProtonTPC_Clone,"Proton","l");
  legProtonTPC->AddEntry((TObject*)0,"with preselection:","");
  legProtonTPC->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  TLegend *legDeuteronTPC = new TLegend(.45,.63,.8,.85);
  legDeuteronTPC->SetFillStyle(0);
  legDeuteronTPC->SetBorderSize(0);
  legDeuteronTPC->AddEntry(fHistTRDPIDDeuteronTPC_clear_Clone,"Deuteron","l");
  legDeuteronTPC->AddEntry((TObject*)0,"with preselection:","");
  legDeuteronTPC->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  legDeuteronTPC->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutDeuteronTPC),"");
  TLegend *legTritonTPC = new TLegend(.45,.63,.8,.85);
  legTritonTPC->SetFillStyle(0);
  legTritonTPC->SetBorderSize(0);
  legTritonTPC->AddEntry(fHistTRDPIDTritonTPC_clear_Clone,"Triton","l");
  legTritonTPC->AddEntry((TObject*)0,"with preselection:","");
  legTritonTPC->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  legTritonTPC->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutTritonTPC),"");
  TLegend *legHelium3TPC = new TLegend(.12,.67,.47,.85);
  legHelium3TPC->SetFillStyle(0);
  legHelium3TPC->SetBorderSize(0);
  legHelium3TPC->AddEntry(fHistTRDPIDHelium3TPC_Clone,"Helium3","l");
  legHelium3TPC->AddEntry((TObject*)0,"with preselection:","");
  legHelium3TPC->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  TLegend *legAlphaTPC = new TLegend(.45,.67,.8,.85);
  legAlphaTPC->SetFillStyle(0);
  legAlphaTPC->SetBorderSize(0);
  legAlphaTPC->AddEntry(fHistTRDPIDAlphaTPC_Clone,"Alpha","l");
  legAlphaTPC->AddEntry((TObject*)0,"with preselection:","");
  legAlphaTPC->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");

  // TOF legends
  TLegend *legTOF = new TLegend(.6,.45,.8,.85);
  legTOF->SetFillStyle(0);
  legTOF->SetBorderSize(0);
  legTOF->AddEntry(fHistTRDPIDElectronTOF,"Electron","l");
  legTOF->AddEntry(fHistTRDPIDPionTOF,"Pion","l");
  legTOF->AddEntry(fHistTRDPIDKaonTOF,"Kaon","l");
  legTOF->AddEntry(fHistTRDPIDProtonTOF,"Proton","l");
  legTOF->AddEntry(fHistTRDPIDHelium3TOF,"Helium3","l");
  legTOF->AddEntry(fHistTRDPIDDeuteronTOF,"Deuteron","l");
  legTOF->AddEntry(fHistTRDPIDTritonTOF,"Triton","l");
  legTOF->AddEntry(fHistTRDPIDAlphaTOF,"Alpha","l");
  TLegend *legElectronTOF = new TLegend(.45,.67,.8,.85);
  legElectronTOF->SetFillStyle(0);
  legElectronTOF->SetBorderSize(0);
  legElectronTOF->AddEntry(fHistTRDPIDElectronTOF_Clone,"Electron","l");
  legElectronTOF->AddEntry((TObject*)0,"with preselection:","");
  legElectronTOF->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");
  TLegend *legPionTOF = new TLegend(.45,.67,.8,.85);
  legPionTOF->SetFillStyle(0);
  legPionTOF->SetBorderSize(0);
  legPionTOF->AddEntry(fHistTRDPIDPionTOF_Clone,"Pion","l");
  legPionTOF->AddEntry((TObject*)0,"with preselection:","");
  legPionTOF->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");
  TLegend *legKaonTOF = new TLegend(.45,.67,.8,.85);
  legKaonTOF->SetFillStyle(0);
  legKaonTOF->SetBorderSize(0);
  legKaonTOF->AddEntry(fHistTRDPIDKaonTOF_Clone,"Kaon","l");
  legKaonTOF->AddEntry((TObject*)0,"with preselection:","");
  legKaonTOF->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");
  TLegend *legProtonTOF = new TLegend(.45,.67,.8,.85);
  legProtonTOF->SetFillStyle(0);
  legProtonTOF->SetBorderSize(0);
  legProtonTOF->AddEntry(fHistTRDPIDProtonTOF_Clone,"Proton","l");
  legProtonTOF->AddEntry((TObject*)0,"with preselection:","");
  legProtonTOF->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");
  TLegend *legDeuteronTOF = new TLegend(.45,.67,.8,.85);
  legDeuteronTOF->SetFillStyle(0);
  legDeuteronTOF->SetBorderSize(0);
  legDeuteronTOF->AddEntry(fHistTRDPIDDeuteronTOF_Clone,"Deuteron","l");
  legDeuteronTOF->AddEntry((TObject*)0,"with preselection:","");
  legDeuteronTOF->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");
  TLegend *legTritonTOF = new TLegend(.45,.67,.8,.85);
  legTritonTOF->SetFillStyle(0);
  legTritonTOF->SetBorderSize(0);
  legTritonTOF->AddEntry(fHistTRDPIDTritonTOF_Clone,"Triton","l");
  legTritonTOF->AddEntry((TObject*)0,"with preselection:","");
  legTritonTOF->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");
  TLegend *legHelium3TOF = new TLegend(.45,.67,.8,.85);
  legHelium3TOF->SetFillStyle(0);
  legHelium3TOF->SetBorderSize(0);
  legHelium3TOF->AddEntry(fHistTRDPIDHelium3TOF_Clone,"Helium3","l");
  legHelium3TOF->AddEntry((TObject*)0,"with preselection:","");
  legHelium3TOF->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");
  TLegend *legAlphaTOF = new TLegend(.45,.67,.8,.85);
  legAlphaTOF->SetFillStyle(0);
  legAlphaTOF->SetBorderSize(0);
  legAlphaTOF->AddEntry(fHistTRDPIDAlphaTOF_Clone,"Alpha","l");
  legAlphaTOF->AddEntry((TObject*)0,"with preselection:","");
  legAlphaTOF->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");

  // TPC & TOF legends
  TLegend *legTPCTOF = new TLegend(.6,.45,.8,.85);
  legTPCTOF->SetFillStyle(0);
  legTPCTOF->SetBorderSize(0);
  legTPCTOF->AddEntry(fHistTRDPIDElectronTPCTOF,"Electron","l");
  legTPCTOF->AddEntry(fHistTRDPIDPionTPCTOF,"Pion","l");
  legTPCTOF->AddEntry(fHistTRDPIDKaonTPCTOF,"Kaon","l");
  legTPCTOF->AddEntry(fHistTRDPIDProtonTPCTOF,"Proton","l");
  legTPCTOF->AddEntry(fHistTRDPIDHelium3TPCTOF,"Helium3","l");
  legTPCTOF->AddEntry(fHistTRDPIDDeuteronTPCTOF_clear,"Deuteron","l");
  legTPCTOF->AddEntry(fHistTRDPIDTritonTPCTOF_clear,"Triton","l");
  legTPCTOF->AddEntry(fHistTRDPIDAlphaTPCTOF,"Alpha","l");
  TLegend *legElectronTPCTOF = new TLegend(.45,.61,.8,.85);
  legElectronTPCTOF->SetFillStyle(0);
  legElectronTPCTOF->SetBorderSize(0);
  legElectronTPCTOF->AddEntry(fHistTRDPIDElectronTPCTOF_Clone,"Electron","l");
  legElectronTPCTOF->AddEntry((TObject*)0,"with preselection:","");
  legElectronTPCTOF->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  legElectronTPCTOF->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");
  TLegend *legPionTPCTOF = new TLegend(.45,.61,.8,.85);
  legPionTPCTOF->SetFillStyle(0);
  legPionTPCTOF->SetBorderSize(0);
  legPionTPCTOF->AddEntry(fHistTRDPIDPionTPCTOF_Clone,"Pion","l");
  legPionTPCTOF->AddEntry((TObject*)0,"with preselection:","");
  legPionTPCTOF->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  legPionTPCTOF->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");
  TLegend *legKaonTPCTOF = new TLegend(.45,.61,.8,.85);
  legKaonTPCTOF->SetFillStyle(0);
  legKaonTPCTOF->SetBorderSize(0);
  legKaonTPCTOF->AddEntry(fHistTRDPIDKaonTPCTOF_Clone,"Kaon","l");
  legKaonTPCTOF->AddEntry((TObject*)0,"with preselection:","");
  legKaonTPCTOF->AddEntry((TObject*)0,Form("|TPCnSigma| = %i#sigma",Sigma),"");
  legKaonTPCTOF->AddEntry((TObject*)0,Form("|TOFnSigma| = %i#sigma",Sigma),"");
  TLegend *legProtonTPCTOF = new TLegend(.45,.61,.8,.85);
  legProtonTPCTOF->SetFillStyle(0);
  legProtonTPCTOF->SetBorderSize(0);
  legProtonTPCTOF->AddEntry(fHistTRDPIDProtonTPCTOF_Clone,"Proton","l");
  legProtonTPCTOF->AddEntry((TObject*)0,"with preselection:","");
  legProtonTPCTOF->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  legProtonTPCTOF->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");
  TLegend *legDeuteronTPCTOF = new TLegend(.45,.58,.8,.85);
  legDeuteronTPCTOF->SetFillStyle(0);
  legDeuteronTPCTOF->SetBorderSize(0);
  legDeuteronTPCTOF->AddEntry(fHistTRDPIDDeuteronTPCTOF_clear_Clone,"Deuteron","l");
  legDeuteronTPCTOF->AddEntry((TObject*)0,"with preselection:","");
  legDeuteronTPCTOF->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  legDeuteronTPCTOF->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");
  legDeuteronTPCTOF->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutDeuteronTPCTOF),"");
  TLegend *legTritonTPCTOF = new TLegend(.45,.58,.8,.85);
  legTritonTPCTOF->SetFillStyle(0);
  legTritonTPCTOF->SetBorderSize(0);
  legTritonTPCTOF->AddEntry(fHistTRDPIDTritonTPCTOF_clear_Clone,"Triton","l");
  legTritonTPCTOF->AddEntry((TObject*)0,"with preselection:","");
  legTritonTPCTOF->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  legTritonTPCTOF->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");
  legTritonTPCTOF->AddEntry((TObject*)0,Form("|#frac{#it{p}}{#it{z}}| < %.1f GeV/#it{c}",ClearCutTritonTPCTOF),"");
  TLegend *legHelium3TPCTOF = new TLegend(.12,.61,.47,.85);
  legHelium3TPCTOF->SetFillStyle(0);
  legHelium3TPCTOF->SetBorderSize(0);
  legHelium3TPCTOF->AddEntry(fHistTRDPIDHelium3TPCTOF_Clone,"Helium3","l");
  legHelium3TPCTOF->AddEntry((TObject*)0,"with preselection:","");
  legHelium3TPCTOF->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  legHelium3TPCTOF->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");
  TLegend *legAlphaTPCTOF = new TLegend(.12,.61,.47,.85);
  legAlphaTPCTOF->SetFillStyle(0);
  legAlphaTPCTOF->SetBorderSize(0);
  legAlphaTPCTOF->AddEntry(fHistTRDPIDAlphaTPCTOF_Clone,"Alpha","l");
  legAlphaTPCTOF->AddEntry((TObject*)0,"with preselection:","");
  legAlphaTPCTOF->AddEntry((TObject*)0,Form("|TPCnSigma| < %i#sigma",Sigma),"");
  legAlphaTPCTOF->AddEntry((TObject*)0,Form("|TOFnSigma| < %i#sigma",Sigma),"");


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Create and fill canvases
  ///__________
  // PID vs. pT
    TCanvas *c1 = new TCanvas("c1",Form("PIDpT w TPC preselection (Sigma=%i) for separated particles",Sigma));
    c1->Divide(3,3);
    c1->cd(1)->SetLogz();
      fHistTRDpTvPID->Draw("colz");
      legPID->Draw("same");
    c1->cd(2)->SetLogz();
      fHistTRDpTvPIDvElectronTPC->Draw("colz");
      legPIDTPCnSigma->Draw("same");
    c1->cd(3)->SetLogz();
      fHistTRDpTvPIDvPionTPC->Draw("colz");
      legPIDTPCnSigma->Draw("same");
    c1->cd(4)->SetLogz();
      fHistTRDpTvPIDvKaonTPC->Draw("colz");
      legPIDTPCnSigma->Draw("same");
    c1->cd(5)->SetLogz();
      fHistTRDpTvPIDvProtonTPC->Draw("colz");
      legPIDTPCnSigma->Draw("same");
    c1->cd(6)->SetLogz();
    //fHistTRDpTvPIDvDeuteronTPC->Draw("colz");
      fHistTRDpTvPIDvDeuteronTPC_clear->Draw("colz");
      legPIDTPCnSigmaDeuteron->Draw("same");
    c1->cd(7)->SetLogz();
    //fHistTRDpTvPIDvTritonTPC->Draw("colz");
      fHistTRDpTvPIDvTritonTPC_clear->Draw("colz");
      legPIDTPCnSigmaTriton->Draw("same");
    c1->cd(8)->SetLogz();
      fHistTRDpTvPIDvHelium3TPC->Draw("colz");
      legPIDTPCnSigma->Draw("same");
    c1->cd(9)->SetLogz();
      fHistTRDpTvPIDvAlphaTPC->Draw("colz");
      legPIDTPCnSigma->Draw("same");
    for(Int_t i = 2; i <= 9; i++) {
      c1->cd(i);
      legTracklets->Draw("same");
      legpTCut->Draw("same");
    }

    TCanvas *c2 = new TCanvas("c2",Form("PIDpT w TOF preselection (Sigma=%i) for separated particles",Sigma));
    c2->Divide(3,3);
    c2->cd(1)->SetLogz();
      fHistTRDpTvPID->Draw("colz");
      legPID->Draw("same");
    c2->cd(2)->SetLogz();
      fHistTRDpTvPIDvElectronTOF->Draw("colz");
      legPIDTOFnSigma->Draw("same");
    c2->cd(3)->SetLogz();
      fHistTRDpTvPIDvPionTOF->Draw("colz");
      legPIDTOFnSigma->Draw("same");
    c2->cd(4)->SetLogz();
      fHistTRDpTvPIDvKaonTOF->Draw("colz");
      legPIDTOFnSigma->Draw("same");
    c2->cd(5)->SetLogz();
      fHistTRDpTvPIDvProtonTOF->Draw("colz");
      legPIDTOFnSigma->Draw("same");
    c2->cd(6)->SetLogz();
      fHistTRDpTvPIDvDeuteronTOF->Draw("colz");
      legPIDTOFnSigma->Draw("same");
    c2->cd(7)->SetLogz();
      fHistTRDpTvPIDvTritonTOF->Draw("colz");
      legPIDTOFnSigma->Draw("same");
    c2->cd(8)->SetLogz();
      fHistTRDpTvPIDvHelium3TOF->Draw("colz");
      legPIDTOFnSigma->Draw("same");
    c2->cd(9)->SetLogz();
      fHistTRDpTvPIDvAlphaTOF->Draw("colz");
      legPIDTOFnSigma->Draw("same");
    for(Int_t i = 2; i <= 9; i++) {
      c2->cd(i);
      legTracklets->Draw("same");
      legpTCut->Draw("same");
    }

    TCanvas *c3 = new TCanvas("c3",Form("PIDpT w TPC & TOF preselection (Sigma=%i) for separated particles",Sigma));
    c3->Divide(3,3);
    c3->cd(1)->SetLogz();
      fHistTRDpTvPID->Draw("colz");
      legPID->Draw("same");
    c3->cd(2)->SetLogz();
      fHistTRDpTvPIDvElectronTPCTOF->Draw("colz");
      legPIDTPCTOFnSigma->Draw("same");
    c3->cd(3)->SetLogz();
      fHistTRDpTvPIDvPionTPCTOF->Draw("colz");
      legPIDTPCTOFnSigma->Draw("same");
    c3->cd(4)->SetLogz();
      fHistTRDpTvPIDvKaonTPCTOF->Draw("colz");
      legPIDTPCTOFnSigma->Draw("same");
    c3->cd(5)->SetLogz();
      fHistTRDpTvPIDvProtonTPCTOF->Draw("colz");
      legPIDTPCTOFnSigma->Draw("same");
    c3->cd(6)->SetLogz();
    //fHistTRDpTvPIDvDeuteronTPCTOF->Draw("colz");
      fHistTRDpTvPIDvDeuteronTPCTOF_clear->Draw("colz");
      legPIDTPCTOFnSigmaDeuteron->Draw("same");
    c3->cd(7)->SetLogz();
    //fHistTRDpTvPIDvTritonTPCTOF->Draw("colz");
      fHistTRDpTvPIDvTritonTPCTOF_clear->Draw("colz");
      legPIDTPCTOFnSigmaTriton->Draw("same");
    c3->cd(8)->SetLogz();
      fHistTRDpTvPIDvHelium3TPCTOF->Draw("colz");
      legPIDTPCTOFnSigma->Draw("same");
    c3->cd(9)->SetLogz();
      fHistTRDpTvPIDvAlphaTPCTOF->Draw("colz");
      legPIDTPCTOFnSigma->Draw("same");
   for(Int_t i = 2; i <= 9; i++) {
      c3->cd(i);
      legTracklets->Draw("same");
      legpTCut->Draw("same");
    }

  // PID vs. pT projection histograms
    TCanvas *c4 = new TCanvas("c4",Form("Proj PIDpT w TPC preselection (Sigma=%i) for separated particles",Sigma));
    c4->Divide(3,3);
    c4->cd(1);
      fHistTRDPIDallTPC->Draw(); // clone of fHistTRDPIDPionTPC because of most entries
      fHistTRDPIDAlphaTPC->Draw("same");
      fHistTRDPIDElectronTPC->Draw("same");
      //fHistTRDPIDDeuteronTPC->Draw("same");
      fHistTRDPIDDeuteronTPC_clear->Draw("same");
      //fHistTRDPIDTritonTPC->Draw("same");
      fHistTRDPIDTritonTPC_clear->Draw("same");
      fHistTRDPIDProtonTPC->Draw("same");
      fHistTRDPIDKaonTPC->Draw("same");
      fHistTRDPIDPionTPC->Draw("same");
      fHistTRDPIDHelium3TPC->Draw("same");
      legTPC->Draw("same");
      legProj->Draw("same");
    c4->cd(2);
      fHistTRDPIDElectronTPC_Clone->Draw();
      legElectronTPC->Draw("same");
    c4->cd(3);
      fHistTRDPIDPionTPC_Clone->Draw();
      legPionTPC->Draw("same");
    c4->cd(4);
      fHistTRDPIDKaonTPC_Clone->Draw();
      legKaonTPC->Draw("same");
    c4->cd(5);
      fHistTRDPIDProtonTPC_Clone->Draw();
      legProtonTPC->Draw("same");
    c4->cd(6);
    //fHistTRDPIDDeuteronTPC_Clone->Draw();
      fHistTRDPIDDeuteronTPC_clear_Clone->Draw();
      legDeuteronTPC->Draw("same");
    c4->cd(7);
    //fHistTRDPIDTritonTPC_Clone->Draw();
      fHistTRDPIDTritonTPC_clear_Clone->Draw();
      legTritonTPC->Draw("same");
    c4->cd(8);
      fHistTRDPIDHelium3TPC_Clone->Draw();
      legHelium3TPC->Draw("same");
    c4->cd(9);
      fHistTRDPIDAlphaTPC_Clone->Draw();
      legAlphaTPC->Draw("same");
    for(Int_t i = 2; i <= 9; i++) {
      c4->cd(i);
      if(i==6 || i==7) {
        legTracklets3->Draw("same");
        legpTCut3->Draw("same");
	continue;
      }
      if(i==8) {
        legTracklets4->Draw("same");
        legpTCut4->Draw("same");
	continue;
      }
      legTracklets2->Draw("same");
      legpTCut2->Draw("same");
    }

    TCanvas *c5 = new TCanvas("c5",Form("Proj PIDpT w TOF preselection (Sigma=%i) for separated particles",Sigma));
    c5->Divide(3,3);
    c5->cd(1);
      fHistTRDPIDallTOF->Draw(); // clone of fHistTRDPIDPionTOF because of most entries
      fHistTRDPIDAlphaTOF->Draw("same");
      fHistTRDPIDElectronTOF->Draw("same");
      fHistTRDPIDDeuteronTOF->Draw("same");
      fHistTRDPIDTritonTOF->Draw("same");
      fHistTRDPIDProtonTOF->Draw("same");
      fHistTRDPIDKaonTOF->Draw("same");
      fHistTRDPIDHelium3TOF->Draw("same");
      legTOF->Draw("same");
      legProj->Draw("same");
    c5->cd(2);
      fHistTRDPIDElectronTOF_Clone->Draw();
      legElectronTOF->Draw("same");
    c5->cd(3);
      fHistTRDPIDPionTOF_Clone->Draw();
      legPionTOF->Draw("same");
    c5->cd(4);
      fHistTRDPIDKaonTOF_Clone->Draw();
      legKaonTOF->Draw("same");
    c5->cd(5);
      fHistTRDPIDProtonTOF_Clone->Draw();
      legProtonTOF->Draw("same");
    c5->cd(6);
      fHistTRDPIDDeuteronTOF_Clone->Draw();
      legDeuteronTOF->Draw("same");
    c5->cd(7);
      fHistTRDPIDTritonTOF_Clone->Draw();
      legTritonTOF->Draw("same");
    c5->cd(8);
      fHistTRDPIDHelium3TOF_Clone->Draw();
      legHelium3TOF->Draw("same");
    c5->cd(9);
      fHistTRDPIDAlphaTOF_Clone->Draw();
      legAlphaTOF->Draw("same");
    for(Int_t i = 2; i <= 9; i++) {
      c5->cd(i);
      legTracklets2->Draw("same");
      legpTCut2->Draw("same");
    }

    TCanvas *c6 = new TCanvas("c6",Form("Proj PIDpT w TPC & TOF preselection (Sigma=%i) for separated particles",Sigma));
    c6->Divide(3,3);
    c6->cd(1);
      fHistTRDPIDallTPCTOF->Draw(); // clone of fHistTRDPIDPionTPCTOF because of most entries
      fHistTRDPIDAlphaTPCTOF->Draw("same");
      fHistTRDPIDElectronTPCTOF->Draw("same");
      //fHistTRDPIDDeuteronTPCTOF->Draw("same");
      fHistTRDPIDDeuteronTPCTOF_clear->Draw("same");
      //fHistTRDPIDTritonTPCTOF->Draw("same");
      fHistTRDPIDTritonTPCTOF_clear->Draw("same");
      fHistTRDPIDProtonTPCTOF->Draw("same");
      fHistTRDPIDKaonTPCTOF->Draw("same");
      fHistTRDPIDPionTPCTOF->Draw("same");
      fHistTRDPIDHelium3TPCTOF->Draw("same");
      legTPCTOF->Draw("same");
      legProj->Draw("same");
    c6->cd(2);
      fHistTRDPIDElectronTPCTOF_Clone->Draw();
      legElectronTPCTOF->Draw("same");
    c6->cd(3);
      fHistTRDPIDPionTPCTOF_Clone->Draw();
      legPionTPCTOF->Draw("same");
    c6->cd(4);
      fHistTRDPIDKaonTPCTOF_Clone->Draw();
      legKaonTPCTOF->Draw("same");
    c6->cd(5);
      fHistTRDPIDProtonTPCTOF_Clone->Draw();
      legProtonTPCTOF->Draw("same");
    c6->cd(6);
    //fHistTRDPIDDeuteronTPCTOF_Clone->Draw();
      fHistTRDPIDDeuteronTPCTOF_clear_Clone->Draw();
      legDeuteronTPCTOF->Draw("same");
    c6->cd(7);
    //fHistTRDPIDTritonTPCTOF_Clone->Draw();
      fHistTRDPIDTritonTPCTOF_clear_Clone->Draw();
      legTritonTPCTOF->Draw("same");
    c6->cd(8);
      fHistTRDPIDHelium3TPCTOF_Clone->Draw();
      legHelium3TPCTOF->Draw("same");
    c6->cd(9);
      fHistTRDPIDAlphaTPCTOF_Clone->Draw();
      legAlphaTPCTOF->Draw("same");
    for(Int_t i = 2; i <= 9; i++) {
      c6->cd(i);
      if(i==6 || i==7) {
        legTracklets3->Draw("same");
        legpTCut3->Draw("same");
	continue;
      }
      if(i==8 || i==9) {
        legTracklets4->Draw("same");
        legpTCut4->Draw("same");
	continue;
      }
      legTracklets2->Draw("same");
      legpTCut2->Draw("same");
    }



    fHistTRDPIDallTPC->SetTitle("");
    fHistTRDPIDallTPC->SetLineWidth(2);
    TCanvas *c7 = new TCanvas("c7",Form("wo_dthe4"));
    c7->Divide(1,1);
    c7->cd(1);
      fHistTRDPIDallTPC->Draw(); // clone of fHistTRDPIDPionTPC because of most entries
      //fHistTRDPIDAlphaTPC->Draw("same");
      fHistTRDPIDElectronTPC->Draw("same");
      //fHistTRDPIDDeuteronTPC->Draw("same");
      //fHistTRDPIDDeuteronTPC_clear->Draw("same");
      //fHistTRDPIDTritonTPC->Draw("same");
      //fHistTRDPIDTritonTPC_clear->Draw("same");
      fHistTRDPIDProtonTPC->Draw("same");
      fHistTRDPIDKaonTPC->Draw("same");
      fHistTRDPIDPionTPC->Draw("same");
      fHistTRDPIDHelium3TPC->Draw("same");
      //legTPC->Draw("same");
      //legProj->Draw("same");

    TCanvas *c8 = new TCanvas("c8",Form("w_dthe4"));
    c8->Divide(1,1);
    c8->cd(1);
      fHistTRDPIDallTPC->Draw(); // clone of fHistTRDPIDPionTPC because of most entries
      fHistTRDPIDAlphaTPC->Draw("same");
      //fHistTRDPIDElectronTPC->Draw("same");
      //fHistTRDPIDDeuteronTPC->Draw("same");
      fHistTRDPIDDeuteronTPC_clear->Draw("same");
      //fHistTRDPIDTritonTPC->Draw("same");
      fHistTRDPIDTritonTPC_clear->Draw("same");
      //fHistTRDPIDProtonTPC->Draw("same");
      //fHistTRDPIDKaonTPC->Draw("same");
      //fHistTRDPIDPionTPC->Draw("same");
      fHistTRDPIDHelium3TPC->Draw("same");
      legTPC1->Draw("same");
      legTPC2->Draw("same");
      //legProj->Draw("same");
      c8->SaveAs(Form("~/TRDProjTPCdt.pdf"));


    TCanvas *c9 = new TCanvas("c9",Form("w_dthe4TPCTOF"));
    c9->Divide(1,1);
    c9->cd(1);
      fHistTRDPIDallTPCTOF->Draw(); // clone of fHistTRDPIDPionTPC because of most entries
      fHistTRDPIDAlphaTPCTOF->Draw("same");
      fHistTRDPIDElectronTPCTOF->Draw("same");
      //fHistTRDPIDDeuteronTPC->Draw("same");
      fHistTRDPIDDeuteronTPCTOF_clear->Draw("same");
      //fHistTRDPIDTritonTPC->Draw("same");
      fHistTRDPIDTritonTPCTOF_clear->Draw("same");
      fHistTRDPIDProtonTPCTOF->Draw("same");
      fHistTRDPIDKaonTPCTOF->Draw("same");
      fHistTRDPIDPionTPCTOF->Draw("same");
      fHistTRDPIDHelium3TPCTOF->Draw("same");
      legTPC->Draw("same");
      legProj->Draw("same");


    TCanvas *TRDChDep = new TCanvas("TRDChDep",Form(""));
    TRDChDep->SetLogz();
      fHistTRDpTvPID->Draw("colz");
      legPID->Draw("same");
    TRDChDep->SaveAs(Form("~/TRDChDep.pdf"));

    TCanvas *TRDChDepDeu = new TCanvas("TRDChDepDeu",Form(""));
    TRDChDepDeu->SetLogz();
      fHistTRDpTvPIDvDeuteronTPC_clear->Draw("colz");
      legPIDTPCnSigmaDeuteron->Draw("same");
      legDeu->Draw("same");
    TRDChDepDeu->SaveAs(Form("~/TRDChDepDeud.pdf"));
    TCanvas *TRDChDepDeuTOF = new TCanvas("TRDChDepDeuTOF",Form(""));
    TRDChDepDeuTOF->SetLogz();
      fHistTRDpTvPIDvDeuteronTPCTOF_clear->Draw("colz");
      legPIDTPCTOFnSigmaDeuteron->Draw("same");
      legDeu->Draw("same");
    TRDChDepDeuTOF->SaveAs(Form("~/TRDChDepDeuTOFd.pdf"));

    TCanvas *TRDChDepTri = new TCanvas("TRDChDepTri",Form(""));
    TRDChDepTri->SetLogz();
      fHistTRDpTvPIDvTritonTPC_clear->Draw("colz");
      legPIDTPCnSigmaTriton->Draw("same");
      legTri->Draw("same");
    TRDChDepTri->SaveAs(Form("~/TRDChDepTrit.pdf"));
    TCanvas *TRDChDepTriTOF = new TCanvas("TRDChDepTriTOF",Form(""));
    TRDChDepTriTOF->SetLogz();
      fHistTRDpTvPIDvTritonTPCTOF_clear->Draw("colz");
      legPIDTPCTOFnSigmaTriton->Draw("same");
      legTri->Draw("same");
    TRDChDepTriTOF->SaveAs(Form("~/TRDChDepTriTOFt.pdf"));

    TCanvas *TRDChDepHe3 = new TCanvas("TRDChDepHe3",Form(""));
    TRDChDepHe3->SetLogz();
      fHistTRDpTvPIDvHelium3TPC->Draw("colz");
      legPIDTPCnSigma->Draw("same");
      legHe3->Draw("same");
    TRDChDepHe3->SaveAs(Form("~/TRDChDepHe3.pdf"));
    TCanvas *TRDChDepHe3TOF = new TCanvas("TRDChDepHe3TOF",Form(""));
    TRDChDepHe3TOF->SetLogz();
      fHistTRDpTvPIDvHelium3TPCTOF->Draw("colz");
      legPIDTPCTOFnSigma->Draw("same");
      legHe3->Draw("same");
    TRDChDepHe3TOF->SaveAs(Form("~/TRDChDepHe3TOF.pdf"));

    TCanvas *TRDChDepAlp = new TCanvas("TRDChDepAlp",Form(""));
    TRDChDepAlp->SetLogz();
      fHistTRDpTvPIDvAlphaTPC->Draw("colz");
      legPIDTPCnSigma->Draw("same");
      legAlp->Draw("same");
    TRDChDepAlp->SaveAs(Form("~/TRDChDepAlp.pdf"));
    TCanvas *TRDChDepAlpTOF = new TCanvas("TRDChDepAlpTOF",Form(""));
    TRDChDepAlpTOF->SetLogz();
      fHistTRDpTvPIDvAlphaTPCTOF->Draw("colz");
      legPIDTPCTOFnSigma->Draw("same");
      legAlp->Draw("same");
    TRDChDepAlpTOF->SaveAs(Form("~/TRDChDepAlpTOF.pdf"));


    TCanvas *TRDProjDeuTOF = new TCanvas("TRDProjDeuTOF",Form(""));
    TRDProjDeuTOF->SetLogz();
      fHistTRDPIDDeuteronTPCTOF_clear->Draw();
      legPIDTPCTOFnSigmaDeuteronright->Draw("same");
      legDeu->Draw("same");
    TRDProjDeuTOF->SaveAs(Form("~/TRDProjDeuTOFd.pdf"));

    TCanvas *TRDProjTriTOF = new TCanvas("TRDProjTriTOF",Form(""));
    TRDProjTriTOF->SetLogz();
      fHistTRDPIDTritonTPCTOF_clear->Draw();
      legPIDTPCTOFnSigmaTritonright->Draw("same");
      legTri->Draw("same");
    TRDProjTriTOF->SaveAs(Form("~/TRDProjTriTOFt.pdf"));

    TCanvas *TRDProjHe3TOF = new TCanvas("TRDProjHe3TOF",Form(""));
    TRDProjHe3TOF->SetLogz();
      fHistTRDPIDHelium3TPCTOF->Draw();
      legPIDTPCTOFnSigma->Draw("same");
      legHe3->Draw("same");
    TRDProjHe3TOF->SaveAs(Form("~/TRDProjHe3TOF.pdf"));

    TCanvas *TRDProjAlpTOF = new TCanvas("TRDProjAlpTOF",Form(""));
    TRDProjAlpTOF->SetLogz();
      fHistTRDPIDAlphaTPCTOF->Draw();
      legPIDTPCTOFnSigma->Draw("same");
      legAlp->Draw("same");
    TRDProjAlpTOF->SaveAs(Form("~/TRDProjAlpTOF.pdf"));


     TCanvas *TRDProjHe3TPC = new TCanvas("TRDProjHe3TPC",Form(""));
    TRDProjHe3TPC->SetLogz();
      fHistTRDPIDHelium3TPC->Draw();
      legPIDTPCnSigma->Draw("same");
      legHe3->Draw("same");
    TRDProjHe3TPC->SaveAs(Form("~/TRDProjHe3TPC.pdf"));

    TCanvas *TRDProjAlpTPC = new TCanvas("TRDProjAlpTPC",Form(""));
    TRDProjAlpTPC->SetLogz();
      fHistTRDPIDAlphaTPC->Draw();
      legPIDTPCnSigma->Draw("same");
      legAlp->Draw("same");
    TRDProjAlpTPC->SaveAs(Form("~/TRDProjAlpTPC.pdf"));

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Save canvases
  ///________________________________
  for(Int_t f = 0; f < nFiles; f++) {
    if(AllOrAnti == "All") {
      c1 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma%s/bbrudnyj_PIDTPC%s%iSigma%s%s.%s",Collisions,Tracklets,Sigma,pT2,Tracklets,Sigma,Collisions,pT2,File[f]));
      c2 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma%s/bbrudnyj_PIDTOF%s%iSigma%s%s.%s",Collisions,Tracklets,Sigma,pT2,Tracklets,Sigma,Collisions,pT2,File[f]));
      c3 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma%s/bbrudnyj_PIDTPCTOF%s%iSigma%s%s.%s",Collisions,Tracklets,Sigma,pT2,Tracklets,Sigma,Collisions,pT2,File[f]));
      c4 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma%s/bbrudnyj_ProjPIDTPC%s%iSigma%s%s.%s",Collisions,Tracklets,Sigma,pT2,Tracklets,Sigma,Collisions,pT2,File[f]));
      c5 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma%s/bbrudnyj_ProjPIDTOF%s%iSigma%s%s.%s",Collisions,Tracklets,Sigma,pT2,Tracklets,Sigma,Collisions,pT2,File[f]));
      c6 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma%s/bbrudnyj_ProjPIDTPCTOF%s%iSigma%s%s.%s",Collisions,Tracklets,Sigma,pT2,Tracklets,Sigma,Collisions,pT2,File[f]));
    }
    else if(AllOrAnti == "Anti") {
      c1 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma%s/bbrudnyj_PIDTPC%s%iSigma%s%s%s.%s",Collisions,Tracklets,Sigma,pT2,Tracklets,Sigma,Collisions,AllOrAnti,pT2,File[f]));
      c2 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma%s/bbrudnyj_PIDTOF%s%iSigma%s%s%s.%s",Collisions,Tracklets,Sigma,pT2,Tracklets,Sigma,Collisions,AllOrAnti,pT2,File[f]));
      c3 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma%s/bbrudnyj_PIDTPCTOF%s%iSigma%s%s%s.%s",Collisions,Tracklets,Sigma,pT2,Tracklets,Sigma,Collisions,AllOrAnti,pT2,File[f]));
      c4 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma%s/bbrudnyj_ProjPIDTPC%s%iSigma%s%s%s.%s",Collisions,Tracklets,Sigma,pT2,Tracklets,Sigma,Collisions,AllOrAnti,pT2,File[f]));
      c5 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma%s/bbrudnyj_ProjPIDTOF%s%iSigma%s%s%s.%s",Collisions,Tracklets,Sigma,pT2,Tracklets,Sigma,Collisions,AllOrAnti,pT2,File[f]));
      c6 ->SaveAs(Form("/lustre/nyx/alice/users/bbrudnyj/codes/%s/%s/bbrudnyj_ReadTreeTestTRDTrigger2/plots%iSigma%s/bbrudnyj_ProjPIDTPCTOF%s%iSigma%s%s%s.%s",Collisions,Tracklets,Sigma,pT2,Tracklets,Sigma,Collisions,AllOrAnti,pT2,File[f]));
    }
  }

  printf("\ndone!\n\n");
}
