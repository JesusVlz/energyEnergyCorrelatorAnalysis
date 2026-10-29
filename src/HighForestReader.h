// Reader for jet trees from CMS data
//
//===========================================================
// HighForestReader.h
//
// Author: Jussi Viinikainen
//===========================================================

#ifndef HIGHFORESTREADER_H
#define HIGHFORESTREADER_H

// Own includes
#include "ForestReader.h"

using namespace std;

class HighForestReader : public ForestReader{
  
private:
  static const Int_t fnMaxJet = 250;        // Maximum number of jets in an event
  static const Int_t fnMaxTrack = 60000;    // Maximum number of tracks in an event
  
public:
  
  // Constructors and destructors
  HighForestReader();                                              // Default constructor
  HighForestReader(Int_t dataType, Int_t useJetTrigger, Int_t jetType, Int_t jetAxis, Int_t matchJets, Bool_t readTrackTree = true, Bool_t mixingMode = false, Bool_t megaSkimMode = false); // Custom constructor
  HighForestReader(const HighForestReader& in);                    // Copy constructor
  virtual ~HighForestReader();                                     // Destructor
  HighForestReader& operator=(const HighForestReader& obj);        // Equal sign operator
  
  // Methods
  void ReadForestFromFile(TFile* inputFile);  // Read the forest from a file
  void ReadForestFromFileList(std::vector<TString> fileList);  // Read the forest from a file list
  void BurnForest();                           // Burn the forest  
  void GetEvent(Int_t nEvent);                 // Get the nEventh event from the file
  bool CheckFileProblems();                    // Check if there are problems in the file list
  
  // Getters for leaves in jet tree
  Float_t GetJetPt(Int_t iJet) const;         // Getter for jet pT
  Float_t GetJetPhi(Int_t iJet) const;        // Getter for jet phi
  Float_t GetJetEta(Int_t iJet) const;        // Getter for jet eta
  Float_t GetJetRawPt(Int_t iJet) const;      // Getter for jet raw pT
  Float_t GetJetMaxTrackPt(Int_t iJet) const; // Getter for maximum track pT inside a jet
  
  void EnableRun3JetSelection(Bool_t enabled = true) override;
  EECJetSelection::JetInputs GetJetSelectionInputs(Int_t iJet) const override;

  // Getters for leaves in the track tree
  Float_t GetTrackPt(Int_t iTrack) const;                    // Getter for track pT
  Float_t GetTrackPtError(Int_t iTrack) const;               // Getter for track pT error
  Float_t GetTrackPhi(Int_t iTrack) const;                   // Getter for track phi
  Float_t GetTrackEta(Int_t iTrack) const;                   // Getter for track eta
  Bool_t GetTrackHighPurity(Int_t iTrack) const;             // Getter for the high purity of the track
  Float_t GetTrackVertexDistanceZ(Int_t iTrack) const;       // Getter for track distance from primary vertex in z-direction
  Float_t GetTrackVertexDistanceZError(Int_t iTrack) const;  // Getter for error of track distance from primary vertex in z-direction
  Float_t GetTrackVertexDistanceXY(Int_t iTrack) const;      // Getter for track distance from primary vertex in xy-direction
  Float_t GetTrackVertexDistanceXYError(Int_t iTrack) const; // Getter for error of track distance from primary vertex in xy-direction
  Float_t GetTrackNormalizedChi2(Int_t iTrack) const;        // Getter for normalized track chi2 value from reconstruction fit
  Float_t GetTrackChi2(Int_t iTrack) const;                  // Getter for track chi2 value from reconstruction fit
  Int_t GetNTrackDegreesOfFreedom(Int_t iTrack) const;       // Getter for number of degrees of freedom in reconstruction fit
  Int_t GetNHitsTrackerLayer(Int_t iTrack) const;            // Getter for number of hits in tracker layers
  Int_t GetNHitsTrack(Int_t iTrack) const;                   // Getter for number of hits for the track
  Float_t GetTrackEnergyEcal(Int_t iTrack) const;            // Getter for track energy in ECal
  Float_t GetTrackEnergyHcal(Int_t iTrack) const;            // Getter for track energy in HCal
  Int_t GetTrackCharge(Int_t iTrack) const;                  // Getter for track charge (relevant only for generator level tracks)
  Int_t GetTrackSubevent(Int_t iTrack) const;                // Getter for track subevent index (relevant only for generator level tracks)
  Int_t GetTrackMCStatus(Int_t iTrack) const;                // Getter for track MC status (only for generator level tracks)
  
  Int_t GetTrackAlgorithm(Int_t iTrack) const;               // Getter for track algorithm
  Int_t GetTrackOriginalAlgorithm(Int_t iTrack) const;       // Getter for track original algorithm
  Float_t GetTrackMVA(Int_t iTrack) const;                   // Getter for track MVA

  // Getters for universally generator level quantities
  Int_t GetNParticles() const;                               // Getter for number of generator level particles
  Float_t GetParticlePt(Int_t iParticle) const;              // Getter for generator level particle pT
  Float_t GetParticlePhi(Int_t iParticle) const;             // Getter for generator level particle phi
  Float_t GetParticleEta(Int_t iParticle) const;             // Getter for generator level particle eta
  Int_t GetParticleCharge(Int_t iParticle) const;            // Getter for generator level particle charge
  Int_t GetParticleSubevent(Int_t iParticle) const;          // Getter for generator level particle subevent index
  
  // Check if generator level jet has a matching reconstructed jet
  Bool_t HasMatchingJet(Int_t iJet) const;      // Check if generator level jet has a matching reconstructed jet
  Float_t GetMatchedPt(Int_t iJet) const;       // Getter for matched generator level jet pT
  Float_t GetMatchedEta(Int_t iJet) const;      // Getter for matched generator level jet eta
  Float_t GetMatchedPhi(Int_t iJet) const;      // Getter for matched generator level jet phi
  Int_t GetPartonFlavor(Int_t iJet) const;      // Parton flavor for the parton initiating the jet
  
private:
  
  // Methods
  void Initialize();       // Connect the branches to the tree
  Int_t GetMatchingIndex(Int_t iJet) const; // Get the matching generator level jet index for the given reconstructed jet

  // Trees in the forest
  TChain* fHeavyIonTree;    // Tree for heavy ion event information
  TChain* fJetTree;         // Tree for jet information
  TChain* fHltTree;         // Tree for high level trigger information
  TChain* fSkimTree;        // Tree for event selection information
  TChain* fTrackTree;       // Tree for tracks.  miniAOD: ppTrack/trackTree. AOD: PbPbTracks/trackTree
  TChain* fParticleTree;    // Tree for generator level particles. Needed for event mixing in MC
  
  // Non-common branches for all types of trees
  TBranch* fnJetsBranch;         // Branch for number of jets in an event
  TBranch* fnMatchedJetsBranch;  // Branch for the number of generator level jets in an event
  TBranch* fnTracksBranch;       // Branch for number of tracks
  TBranch* fMatchedJetWTAEtaBranch; // Branch for matched jet eta with respect to WTA axis
  TBranch* fMatchedJetWTAPhiBranch; // Branch for matched jet phi with respect to WTA axis
  
  TBranch* fTrackAlgorithmBranch;         // Branch for track algorithm
  TBranch* fTrackOriginalAlgorithmBranch; // Branch for track original algorithm
  TBranch* fTrackMVABranch;               // Branch for track MVA

  TBranch* fParticlePtBranch;             // Branch for particle pT:s
  TBranch* fParticleEtaBranch;            // Branch for particle phis
  TBranch* fParticlePhiBranch;            // Branch for particle etas
  TBranch* fParticleChargeBranch;         // Branch for particle charges
  TBranch* fParticleSubeventBranch;       // Branch for particle subevent indices (0 = PYTHIA, (>0) = HYDJET)
  
  EECJetSelection::Detail::NativePFJetInputs<fnMaxJet> fRun3PFJetInputs;
  Bool_t fRun3JetInputsBound = false;
  Int_t fRun3ValidatedJetTreeNumber = -1;
  Float_t fDetectorJetEtaArray[fnMaxJet] = {0};
  Float_t fDetectorJetPhiArray[fnMaxJet] = {0};

  // Leaves for jet tree
  Float_t fJetPtArray[fnMaxJet] = {0};            // pT:s of all the jets in an event
  Float_t fJetPhiArray[fnMaxJet] = {0};           // phis of all the jets in an event
  Float_t fJetEtaArray[fnMaxJet] = {0};           // etas of all the jets in an event
  Float_t fJetRawPtArray[fnMaxJet] = {0};         // raw jet pT for all the jets in an event
  Float_t fJetMaxTrackPtArray[fnMaxJet] = {0};    // maximum track pT inside a jet for all the jets in an event
  Float_t fJetRefPtArray[fnMaxJet] = {0};         // reference generator level pT for a reconstructed jet
  Float_t fJetRefEtaArray[fnMaxJet] = {0};        // reference generator level eta for a reconstructed jet
  Float_t fJetRefPhiArray[fnMaxJet] = {0};        // reference generator level phi for a reconstructed jet
  Int_t fJetRefFlavorArray[fnMaxJet] = {0};       // flavor for initiating parton for the reference gen jet
  Float_t fMatchedJetPtArray[fnMaxJet] = {0};     // pT:s of all the generator level jets in an event
  Float_t fMatchedJetPhiArray[fnMaxJet] = {0};    // phis of all the generator level jets in an event
  Float_t fMatchedJetWTAPhiArray[fnMaxJet] = {0}; // phis of all the generator level jets in an event with respect to WTA axis
  Float_t fMatchedJetEtaArray[fnMaxJet] = {0};    // etas of all the generator level jets in an event
  Float_t fMatchedJetWTAEtaArray[fnMaxJet] = {0}; // etas of all the generator level jets in an event with respect to WTA axis
  
  // Leaves for the track tree in AOD forests
  Float_t fTrackPtArray[fnMaxTrack] = {0};                    // Array for track pT:s
  Float_t fTrackPtErrorArray[fnMaxTrack] = {0};               // Array for track pT errors
  Float_t fTrackPhiArray[fnMaxTrack] = {0};                   // Array for track phis
  Float_t fTrackEtaArray[fnMaxTrack] = {0};                   // Array for track etas
  Bool_t fHighPurityTrackArray[fnMaxTrack] = {0};             // Array for the high purity of tracks
  Float_t fTrackVertexDistanceZArray[fnMaxTrack] = {0};       // Array for track distance from primary vertex in z-direction
  Float_t fTrackVertexDistanceZErrorArray[fnMaxTrack] = {0};  // Array for error for track distance from primary vertex in z-direction
  Float_t fTrackVertexDistanceXYArray[fnMaxTrack] = {0};      // Array for track distance from primary vertex in xy-direction
  Float_t fTrackVertexDistanceXYErrorArray[fnMaxTrack] = {0}; // Array for error for track distance from primary vertex in xy-direction
  Float_t fTrackChi2Array[fnMaxTrack] = {0};                  // Array for track chi2 value from reconstruction fit
  UChar_t fnTrackDegreesOfFreedomArray[fnMaxTrack] = {0};     // Array for number of degrees of freedom in reconstruction fit
  UChar_t fnHitsTrackerLayerArray[fnMaxTrack] = {0};          // Array for number of hits in tracker layers
  UChar_t fnHitsTrackArray[fnMaxTrack] = {0};                 // Array for number of hits for the track
  Float_t fTrackEnergyEcalArray[fnMaxTrack] = {0};            // Array for track energy in ECal
  Float_t fTrackEnergyHcalArray[fnMaxTrack] = {0};            // Array for track energy in HCal
  
  UChar_t fTrackAlgorithmArray[fnMaxTrack] = {0};             // Array for track algorithm
  UChar_t fTrackOriginalAlgorithmArray[fnMaxTrack] = {0};     // Array for track original algorithm
  Float_t fTrackMVAArray[fnMaxTrack] = {0};                   // Array for track MVA
  
  // Leaves for the track tree in MiniAOD forests
  vector<float>* fTrackPtVector;                    // Vector for track pT:s
  vector<float>* fTrackPtErrorVector;               // Vector for track pT errors
  vector<float>* fTrackPhiVector;                   // Vector for track phis
  vector<float>* fTrackEtaVector;                   // Vector for track etas
  vector<bool>* fHighPurityTrackVector;             // Vector for the high purity of tracks
  vector<float>* fTrackVertexDistanceZVector;       // Vector for track distance from primary vertex in z-direction
  vector<float>* fTrackVertexDistanceZErrorVector;  // Vector for error for track distance from primary vertex in z-direction
  vector<float>* fTrackVertexDistanceXYVector;      // Vector for track distance from primary vertex in xy-direction
  vector<float>* fTrackVertexDistanceXYErrorVector; // Vector for error for track distance from primary vertex in xy-direction
  vector<float>* fTrackNormalizedChi2Vector;        // Vector for normalized track chi2 value from reconstruction fit
  vector<char>* fnHitsTrackerLayerVector;           // Vector for number of hits in tracker layers
  vector<char>* fnHitsTrackVector;                  // Vector for number of hits for the track
  vector<float>* fTrackEnergyEcalVector;            // Vector for track energy in ECal
  vector<float>* fTrackEnergyHcalVector;            // Vector for track energy in HCal

  // Leaves for the generator level particle tree
  vector<float>* fParticlePtVector;       // Vector for particle pT:s
  vector<float>* fParticlePhiVector;      // Vector for particle phis
  vector<float>* fParticleEtaVector;      // Vector for particle etas
  vector<int>* fParticleChargeVector;     // Vector for particle charges
  vector<int>* fParticleSubeventVector;   // Vector for particle subevent indices (0 = PYTHIA, (>0) = HYDJET)

};

#endif
