// Reader for jet and track trees from CMS data
//
//===========================================================
// UnfoldingForestReader.h
//
// Author: Jussi Viinikainen
//===========================================================

#ifndef UNFOLDINGFORESTREADER_H
#define UNFOLDINGFORESTREADER_H

#include "Run3JetSelectionInputs.h"

// C++ includes
#include <iostream>
#include <assert.h>
#include <vector>

// Root includes
#include <TString.h>
#include <TTree.h>
#include <TChain.h>
#include <TBranch.h>
#include <TFile.h>

using namespace std;

class UnfoldingForestReader{
  
private:
  static const Int_t fnMaxJet = 250;        // Maximum number of jets in an event
  static const Int_t fnMaxTrack = 60000;    // Maximum number of tracks in an event
  
public:
  
  // Possible data types to be read with the reader class
  // Keep the values synchronized with ForestReader. pPb is not supported here.
  enum enumDataTypes{kPp = 0, kPbPb = 1, kPpMC = 2, kPbPbMC = 3,
                     kOO = 9, kOOMC = 10, kPpRef5p36TeV = 11,
                     kPpRef5p36TeVMC = 12, knDataTypes = 13};
  
  // Constructors and destructors
  UnfoldingForestReader();                                          // Default constructor
  UnfoldingForestReader(Int_t dataType, Int_t jetType, Int_t jetAxis); // Custom constructor
  UnfoldingForestReader(const UnfoldingForestReader& in);                    // Copy constructor
  ~UnfoldingForestReader();                                 // Destructor
  UnfoldingForestReader& operator=(const UnfoldingForestReader& obj);        // Equal sign operator
  
  // Methods
  void GetEvent(Int_t nEvent);                 // Get the nth event in tree
  Int_t GetNEvents() const;                        // Get the number of events
  void ReadForestFromFile(TFile *inputFile);   // Read the forest from a file
  void ReadForestFromFileList(std::vector<TString> fileList);   // Read the forest from a file list
  void BurnForest();                           // Burn the forest
  
  // Getters for leaves in heavy ion tree
  Float_t GetVz() const;              // Getter for vertex z position
  Float_t GetCentrality() const;      // Getter for centrality
  Int_t GetHiBin() const;             // Getter for CMS hiBin
  Float_t GetPtHat() const;           // Getter for pT hat
  Float_t GetEventWeight() const;     // Getter for event weight in MC
  
  // Getters for leaves in reconstructed jet tree
  Int_t GetNJets() const;                     // Getter for number of jets
  Float_t GetJetPt(Int_t iJet) const;         // Getter for jet pT
  Float_t GetJetPhi(Int_t iJet) const;        // Getter for jet phi
  Float_t GetJetEta(Int_t iJet) const;        // Getter for jet eta
  Float_t GetJetRawPt(Int_t iJet) const;      // Getter for jet raw pT
  Float_t GetJetMaxTrackPt(Int_t iJet) const; // Getter for maximum track pT inside a jet
  
  // Reco to gen matching
  Bool_t HasMatchingGenJet(Int_t iJet) const; // Check if reconstructed jet has a matching generator level jet
  Int_t GetMatchingGenIndex(Int_t iJet) const; // Get the matching generator level jet index for the given reconstructed jet
  Float_t GetMatchedGenPt(Int_t iJet) const;  // Getter for matched generator level jet pT
  Float_t GetMatchedGenEta(Int_t iJet) const; // Getter for matched generator level jet eta
  Float_t GetMatchedGenPhi(Int_t iJet) const; // Getter for matched generator level jet phi

  void EnableRun3JetSelection(Bool_t enabled = true);
  EECJetSelection::JetInputs GetJetSelectionInputs(Int_t iJet) const;

  // Getters for leaves in generator level jet tree
  Int_t GetNGeneratorJets() const;                   // Getter for number of generator level jets
  Int_t GetMatchingRecoIndex(Int_t iJet) const; // Get the matching reconstructed jet index for the given generator level jet
  Float_t GetGeneratorJetPt(Int_t iJet) const;       // Getter for generator level jet pT
  Float_t GetGeneratorJetPhi(Int_t iJet) const;      // Getter for generator level jet phi
  Float_t GetGeneratorJetEta(Int_t iJet) const;      // Getter for generator level jet eta

  // Gen to reco matching
  Bool_t HasMatchingRecoJet(Int_t iJet) const; // Check if generator level has a matching reconstructed jet
  Float_t GetMatchedRecoPt(Int_t iJet) const;  // Getter for matched reconstructed jet pT
  Float_t GetMatchedRecoEta(Int_t iJet) const; // Getter for matched reconstructed jet eta
  Float_t GetMatchedRecoPhi(Int_t iJet) const; // Getter for matched reconstructed jet phi
  
  // Getters for leaves in track tree
  Int_t GetNTracks() const;                                  // Getter for number of tracks
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
  Int_t GetTrackCharge(Int_t iTrack) const;                  // Getter for track charge
  
  // Getters for leaves in generator level particle tree
  Int_t GetNGenParticles() const;                            // Getter for number of generator level particles
  Float_t GetGenParticlePt(Int_t iTrack) const;              // Getter for generator level particle pT
  Float_t GetGenParticlePhi(Int_t iTrack) const;             // Getter for generator level particle phi
  Float_t GetGenParticleEta(Int_t iTrack) const;             // Getter for generator level particle eta
  Int_t GetGenParticleCharge(Int_t iTrack) const;            // Getter for generator level particle charge
  Int_t GetGenParticleSubevent(Int_t iTrack) const;          // Getter for generator level particle subevent index
  
  // Setter for data type
  void SetDataType(Int_t dataType); // Setter for data type
  
private:
  
  // Methods
  void Initialize();      // Connect the branches to the tree
  Bool_t IsMC() const;
  Bool_t IsRun3() const;
    
  Int_t fDataType;        // Same type identifiers as ForestReader (0-3 and 9-12)
  Int_t fJetType;         // Choose the type of jets usedfor analysis. 0 = Calo jets, 1 = PF jets
  Int_t fJetAxis;         // Jet axis used for the jets. 0 = Anti-kT, 1 = WTA
  Bool_t fIsMiniAOD;      // Flag for type of the forest True = MiniAOD forest, False = AOD forest
  
  // Trees in the forest
  TTree* fHeavyIonTree;    // Tree for heavy ion event information
  TTree* fJetTree;         // Tree for jet information
  TTree* fTrackTree;       // Tree for reconstructed tracks
  TTree* fGenParticleTree;  // Tree for generator level particles

  // Branches for heavy ion tree
  TBranch* fHiVzBranch;                   // Branch for vertex z-position
  TBranch* fHiBinBranch;                  // Branch for centrality
  TBranch* fPtHatBranch;                  // Branch for pT hat
  TBranch* fEventWeightBranch;            // Branch for event weight
  
  // Branches for jet tree
  TBranch* fnJetsBranch;          // Branch for number of jets
  TBranch* fJetPtBranch;          // Branch for jet pT
  TBranch* fJetPhiBranch;         // Branch for jet phi
  TBranch* fJetWTAPhiBranch;      // Branch for jet phi
  TBranch* fJetEtaBranch;         // Branch for jet eta
  TBranch* fJetWTAEtaBranch;      // Branch for jet eta
  TBranch* fJetRawPtBranch;       // Branch for raw jet pT
  TBranch* fJetMaxTrackPtBranch;  // Maximum pT for a track inside a jet

  TBranch* fJetRefPtBranch;      // Branch for reference generator level pT for a reconstructed jet
  TBranch* fJetRefEtaBranch;     // Branch for reference generator level eta for a reconstructed jet
  TBranch* fJetRefPhiBranch;     // Branch for reference generator level phi for a reconstructed jet

  TBranch* fnGenJetsBranch;      // Branch for number of generator level jets
  TBranch* fGenJetPtBranch;      // Branch for generator level jet pT
  TBranch* fGenJetPhiBranch;     // Branch for generator level jet phi
  TBranch* fGenJetWTAPhiBranch;  // Branch for generator level jet WTA phi
  TBranch* fGenJetEtaBranch;     // Branch for generator level jet eta
  TBranch* fGenJetWTAEtaBranch;  // Branch for generator level jet eta

  // Branches for track tree
  TBranch* fnTracksBranch;                     // Branch for number of tracks
  TBranch* fTrackPtBranch;                     // Branch for track pT
  TBranch* fTrackPtErrorBranch;                // Branch for track pT error
  TBranch* fTrackPhiBranch;                    // Branch for track phi
  TBranch* fTrackEtaBranch;                    // Branch for track eta
  TBranch* fHighPurityTrackBranch;             // Branch for high purity of the track
  TBranch* fTrackVertexDistanceZBranch;        // Branch for track distance from primary vertex in z-direction
  TBranch* fTrackVertexDistanceZErrorBranch;   // Branch for error for track distance from primary vertex in z-direction
  TBranch* fTrackVertexDistanceXYBranch;       // Branch for track distance from primary vertex in xy-direction
  TBranch* fTrackVertexDistanceXYErrorBranch;  // Branch for error for track distance from primary vertex in xy-direction
  TBranch* fTrackChi2Branch;                   // Branch for track chi2 value from reconstruction fit
  TBranch* fnTrackDegreesOfFreedomBranch;      // Branch for number of degrees of freedom in reconstruction fit
  TBranch* fnHitsTrackerLayerBranch;           // Branch for number of hits in tracker layers
  TBranch* fnHitsTrackBranch;                  // Branch for number of hits for the track
  TBranch* fTrackEnergyEcalBranch;             // Branch for track energy in ECal
  TBranch* fTrackEnergyHcalBranch;             // Branch for track energy in HCal
  TBranch* fTrackChargeBranch;                 // Branch for track charge

  // Branches for genenerator level particle tree
  TBranch* fGenParticlePtBranch;        // Branch for generator level particle pT:s
  TBranch* fGenParticlePhiBranch;       // Branch for generator level particle phis
  TBranch* fGenParticleEtaBranch;       // Branch for generator level particle etas
  TBranch* fGenParticleChargeBranch;    // Branch for generator level particle charges
  TBranch* fGenParticleSubeventBranch;  // Branch for generator level particle subevent indices (0 = PYTHIA, (>0) = HYDJET)

  // Leaves for heavy ion tree
  Float_t fVertexZ;    // Vertex z-position
  Int_t fHiBin;        // HiBin = Centrality percentile * 2
  Float_t fPtHat;      // pT hat
  
  Bool_t fRun3JetSelectionEnabled = false;
  Bool_t fRun3JetInputsBound = false;
  EECJetSelection::Detail::NativePFJetInputs<fnMaxJet> fRun3PFJetInputs;

  // Leaves for jet tree
  Int_t fnJets;          // number of jets in an event
  Int_t fnGenJets;       // Number of generator level jets in an event
  Float_t fEventWeight;  // jet weight in the MC tree
  
  Float_t fJetPtArray[fnMaxJet] = {0};         // pT:s of all the jets in an event
  Float_t fJetPhiArray[fnMaxJet] = {0};        // phis of all the jets in an event
  Float_t fJetWTAPhiArray[fnMaxJet] = {0};     // WTA phis of all the jets in an event
  Float_t fJetEtaArray[fnMaxJet] = {0};        // etas of all the jets in an event
  Float_t fJetWTAEtaArray[fnMaxJet] = {0};     // WTA etas of all the jets in an event
  Float_t fJetRawPtArray[fnMaxJet] = {0};      // raw jet pT for all the jets in an event
  Float_t fJetMaxTrackPtArray[fnMaxJet] = {0}; // maximum track pT inside a jet for all the jets in an event
   
  Float_t fJetRefPtArray[fnMaxJet] = {0};      // reference generator level pT for a reconstructed jet
  Float_t fJetRefEtaArray[fnMaxJet] = {0};     // reference generator level eta for a reconstructed jet
  Float_t fJetRefPhiArray[fnMaxJet] = {0};     // reference generator level phi for a reconstructed jet

  Float_t fGenJetPtArray[fnMaxJet] = {0};       // pT:s of the generator level jets in an event
  Float_t fGenJetPhiArray[fnMaxJet] = {0};      // phis of the generator level jets in an event
  Float_t fGenJetWTAPhiArray[fnMaxJet] = {0};   // WTA phis of the generator level jets in an event
  Float_t fGenJetEtaArray[fnMaxJet] = {0};      // etas of the generator level jets in an event
  Float_t fGenJetWTAEtaArray[fnMaxTrack] = {0}; // WTA etas of the generator level jets in an event
  
  // Leaves for the track tree regardless of forest type
  Int_t fnTracks;  // Number of tracks
  
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
  Int_t fTrackChargeArray[fnMaxTrack] = {0};                  // Array for track charge

  // Leaves for the track tree in MiniAOD forests
  vector<float>* fTrackPtVector;                     // Vector for track pT:s
  vector<float>* fTrackPtErrorVector;                // Vector for track pT errors
  vector<float>* fTrackPhiVector;                    // Vector for track phis
  vector<float>* fTrackEtaVector;                    // Vector for track etas
  vector<bool>* fHighPurityTrackVector;              // Vector for the high purity of tracks
  vector<float>* fTrackVertexDistanceZVector;        // Vector for track distance from primary vertex in z-direction
  vector<float>* fTrackVertexDistanceZErrorVector;   // Vector for error for track distance from primary vertex in z-direction
  vector<float>* fTrackVertexDistanceXYVector;       // Vector for track distance from primary vertex in xy-direction
  vector<float>* fTrackVertexDistanceXYErrorVector;  // Vector for error for track distance from primary vertex in xy-direction
  vector<float>* fTrackNormalizedChi2Vector;         // Vector for normalized track chi2 value from reconstruction fit
  vector<char>* fnHitsTrackerLayerVector;            // Vector for number of hits in tracker layers
  vector<char>* fnHitsTrackVector;                   // Vector for number of hits for the track
  vector<float>* fTrackEnergyEcalVector;             // Vector for track energy in ECal
  vector<float>* fTrackEnergyHcalVector;             // Vector for track energy in HCal
  vector<char>* fTrackChargeVector;                  // Vector for track charges

  // Leaves for the generator level particle tree
  Int_t fnGenParticles;                    // Number of generator level particles
  vector<float>* fGenParticlePtArray;      // Array for generator level particle pT:s
  vector<float>* fGenParticlePhiArray;     // Array for generator level particle phis
  vector<float>* fGenParticleEtaArray;     // Array for generator level particle etas
  vector<int>* fGenParticleChargeArray;    // Array for generator level particle charges
  vector<int>* fGenParticleSubeventArray;  // Array for generator level particle subevent indices (0 = PYTHIA, (>0) = HYDJET)
};

#endif
