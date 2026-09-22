#include "EECHistogramManager.h" R__LOAD_LIBRARY(plotting/DrawingClasses.so)
#include "EECCard.h"
#include "JDrawer.h"

#include <TLine.h>

/*
 * Compare raw EEC in ppRef and OO.
 *
 * File 0 = ppRef reference
 * File 1 = OO
 *
 * Lower panel:
 *   ppRef / ppRef = 1
 *   OO / ppRef
 */
void compareEECinCustomBins(){

  // ====================================================
  //                    Input files
  // ====================================================

  const int nComparisonFiles = 2;

  TString fileName[nComparisonFiles];

  fileName[0] =
    "/afs/cern.ch/user/j/jvelazqu/EECAnalysis/"
    "CMSSW_15_0_9_patch4/src/"
    "energyEnergyCorrelatorAnalysis/"
    "merged_eecppRef_2024-09-10-processed.root";

  fileName[1] =
    "/afs/cern.ch/user/j/jvelazqu/EECAnalysis/"
    "CMSSW_15_0_9_patch4/src/"
    "energyEnergyCorrelatorAnalysis/"
    "merged_eecOO_2025-09-10-processed.root";


  // ====================================================
  //                    Open files
  // ====================================================

  TFile* inputFile[nComparisonFiles];
  EECCard* card[nComparisonFiles];

  for(int iFile = 0;
      iFile < nComparisonFiles;
      iFile++){

    inputFile[iFile] =
      TFile::Open(fileName[iFile]);

    if(inputFile[iFile] == NULL ||
       inputFile[iFile]->IsZombie()){

      cout
        << "ERROR: Could not open "
        << fileName[iFile]
        << endl;

      return;
    }

    card[iFile] =
      new EECCard(inputFile[iFile]);

    cout
      << "File " << iFile
      << ": "
      << card[iFile]->GetDataType()
      << endl;
  }


  // ====================================================
  //                 Binning configuration
  // ====================================================

  // Centrality:
  //
  // ppRef has no centrality selection -> C0.
  // This value is therefore just a placeholder for file 0.
  //
  // OO uses the requested centrality bin.
  std::pair<double,double>
    comparedCentralityBin[nComparisonFiles];

  comparedCentralityBin[0] =
    std::make_pair(-1,100);  // ppRef inclusive

  comparedCentralityBin[1] =
    std::make_pair(0,30);    // OO 0-30%


  // Jet pT
  std::pair<double,double>
    comparedJetPtBin[nComparisonFiles];

  comparedJetPtBin[0] =
    std::make_pair(60,80);

  comparedJetPtBin[1] =
    std::make_pair(60,80);


  // Track pT threshold
  double comparedTrackPtBin[nComparisonFiles];

  comparedTrackPtBin[0] = 2.0;
  comparedTrackPtBin[1] = 2.0;


  // ====================================================
  //                 Drawing configuration
  // ====================================================

  const bool saveFigures = true;

  const char* figureFormat = "pdf";

  std::pair<double,double>
    drawingRange =
      std::make_pair(0.008,0.39);

  std::pair<double,double>
    ratioZoom =
      std::make_pair(0.8,1.2);

  const double epsilon = 0.0001;


  // ====================================================
  //                Histogram managers
  // ====================================================

  EECHistogramManager*
    histograms[nComparisonFiles];

  for(int iFile = 0;
      iFile < nComparisonFiles;
      iFile++){

    histograms[iFile] =
      new EECHistogramManager(
        inputFile[iFile],
        card[iFile]
      );
  }


  // ====================================================
  //                  Read raw EEC
  // ====================================================

  TH1D*
    hEnergyEnergyCorrelator[nComparisonFiles];

  TH1D*
    hEnergyEnergyCorrelatorRatio[nComparisonFiles];

  for(int iFile = 0;
      iFile < nComparisonFiles;
      iFile++){

    hEnergyEnergyCorrelator[iFile] = NULL;
    hEnergyEnergyCorrelatorRatio[iFile] = NULL;
  }


  for(int iFile = 0;
      iFile < nComparisonFiles;
      iFile++){

    // ------------------------------------------
    // Centrality treatment
    // ------------------------------------------

    const bool isPpLike =
      card[iFile]
        ->GetDataType()
        .Contains("pp") ||
      card[iFile]
        ->GetDataType()
        .Contains("pPb");

    int iCentrality = 0;

    if(!isPpLike){

      iCentrality =
        card[iFile]
          ->FindBinIndexCentrality(
            comparedCentralityBin[iFile]
          );
    }


    // ------------------------------------------
    // Jet and track pT bins
    // ------------------------------------------

    const int iJetPt =
      card[iFile]
        ->FindBinIndexJetPtEEC(
          comparedJetPtBin[iFile]
        );

    const int iTrackPt =
      card[iFile]
        ->GetBinIndexTrackPtEEC(
          comparedTrackPtBin[iFile]
        );


    cout << endl;
    cout
      << "Loading "
      << card[iFile]->GetDataType()
      << endl;

    cout
      << "  C = "
      << iCentrality
      << endl;

    cout
      << "  J = "
      << iJetPt
      << "  ("
      << comparedJetPtBin[iFile].first
      << "-"
      << comparedJetPtBin[iFile].second
      << " GeV)"
      << endl;

    cout
      << "  T = "
      << iTrackPt
      << "  (> "
      << comparedTrackPtBin[iFile]
      << " GeV)"
      << endl;


    if(iCentrality < 0 ||
       iJetPt < 0 ||
       iTrackPt < 0){

      cout
        << "ERROR: Requested bin does not exist!"
        << endl;

      return;
    }


    // ------------------------------------------
    // RAW same-jet EEC
    // ------------------------------------------

    TH1D* sourceHistogram =
      histograms[iFile]
        ->GetHistogramEnergyEnergyCorrelator(
          EECHistogramManager::
            kEnergyEnergyCorrelator,
          iCentrality,
          iJetPt,
          iTrackPt,
          EECHistograms::kSameJetPair
        );


    if(sourceHistogram == NULL){

      cout
        << "ERROR: Raw EEC histogram "
        << "was not found for file "
        << iFile
        << endl;

      return;
    }


    // Clone before scaling.
    // Do not modify the histogram owned by the input file.
    hEnergyEnergyCorrelator[iFile] =
      (TH1D*) sourceHistogram->Clone(
        Form(
          "rawEEC_file%d",
          iFile
        )
      );

    hEnergyEnergyCorrelator[iFile]
      ->SetDirectory(NULL);


    // ------------------------------------------
    // Normalize distribution
    // ------------------------------------------

    const int lowNormalizationBin =
      hEnergyEnergyCorrelator[iFile]
        ->GetXaxis()
        ->FindBin(
          drawingRange.first +
          epsilon
        );

    const int highNormalizationBin =
      hEnergyEnergyCorrelator[iFile]
        ->GetXaxis()
        ->FindBin(
          drawingRange.second -
          epsilon
        );


    const double normalization =
      hEnergyEnergyCorrelator[iFile]
        ->Integral(
          lowNormalizationBin,
          highNormalizationBin,
          "width"
        );


    if(normalization <= 0){

      cout
        << "ERROR: EEC normalization "
        << "is zero for file "
        << iFile
        << endl;

      return;
    }


    hEnergyEnergyCorrelator[iFile]
      ->Scale(
        1.0 / normalization
      );

  } // File loop


  // ====================================================
  //                    Calculate ratios
  // ====================================================

  for(int iFile = 0;
      iFile < nComparisonFiles;
      iFile++){

    hEnergyEnergyCorrelatorRatio[iFile] =
      (TH1D*)
      hEnergyEnergyCorrelator[iFile]
        ->Clone(
          Form(
            "eecRatio%d",
            iFile
          )
        );

    hEnergyEnergyCorrelatorRatio[iFile]
      ->SetDirectory(NULL);

    hEnergyEnergyCorrelatorRatio[iFile]
      ->Divide(
        hEnergyEnergyCorrelator[0]
      );
  }


  // ====================================================
  //                       Drawing
  // ====================================================

  JDrawer* drawer =
    new JDrawer();

  drawer->SetDefaultAppearanceSplitCanvas();
  drawer->SetRelativeCanvasSize(1.1,1.1);
  drawer->SetLeftMargin(0.14);
  drawer->SetTopMargin(0.07);
  drawer->SetTitleOffsetY(1.7);
  drawer->SetTitleOffsetX(1.0);

  drawer->SetLogX(true);

  drawer->CreateSplitCanvas();


  // ====================================================
  //                     Marker styles
  // ====================================================

  // Same convention as reference figure:
  //
  // ppRef = black full circle
  // OO    = red open square
  const int markerStyle[nComparisonFiles] = {
    kFullCircle,
    kOpenSquare
  };

  const int color[nComparisonFiles] = {
    kBlack,
    kRed
  };


  for(int iFile = 0;
      iFile < nComparisonFiles;
      iFile++){

    hEnergyEnergyCorrelator[iFile]
      ->SetMarkerStyle(
        markerStyle[iFile]
      );

    hEnergyEnergyCorrelator[iFile]
      ->SetMarkerColor(
        color[iFile]
      );

    hEnergyEnergyCorrelator[iFile]
      ->SetLineColor(
        color[iFile]
      );


    hEnergyEnergyCorrelatorRatio[iFile]
      ->SetMarkerStyle(
        markerStyle[iFile]
      );

    hEnergyEnergyCorrelatorRatio[iFile]
      ->SetMarkerColor(
        color[iFile]
      );

    hEnergyEnergyCorrelatorRatio[iFile]
      ->SetLineColor(
        color[iFile]
      );
  }


  // ====================================================
  //                     Upper panel
  // ====================================================

  drawer->SetLogY(true);

  hEnergyEnergyCorrelator[0]
    ->GetXaxis()
    ->SetRangeUser(
      drawingRange.first,
      drawingRange.second
    );


  drawer->DrawHistogramToUpperPad(
    hEnergyEnergyCorrelator[0],
    "#Deltar",
    "EEC (raw)",
    " "
  );


  hEnergyEnergyCorrelator[1]
    ->Draw("same");


  // ====================================================
  //                         Legend
  // ====================================================

  TLegend* legend =
    new TLegend(
      0.18,
      0.04,
      0.55,
      0.55
    );

  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextSize(0.05);
  legend->SetTextFont(62);


  // Common jet pT selection
  legend->AddEntry(
    (TObject*)0,
    Form(
      "%.0f < jet p_{T} < %.0f GeV",
      comparedJetPtBin[0].first,
      comparedJetPtBin[0].second
    ),
    ""
  );


  // Common track pT selection
  legend->AddEntry(
    (TObject*)0,
    Form(
      "p_{T}^{ch} > %.1f GeV",
      comparedTrackPtBin[0]
    ),
    ""
  );


  TString energyWeightString =
    (
      card[0]->GetWeightExponent() == 1
    )
      ? "Nominal energy weight"
      : "Energy weight squared";


  legend->AddEntry(
    (TObject*)0,
    energyWeightString,
    ""
  );


  // ppRef
  legend->AddEntry(
    hEnergyEnergyCorrelator[0],
    "pp 5.36 TeV, raw",
    "p"
  );


  // OO
  legend->AddEntry(
    hEnergyEnergyCorrelator[1],
    Form(
      "OO 5.36 TeV, %.0f-%.0f%%, raw",
      comparedCentralityBin[1].first,
      comparedCentralityBin[1].second
    ),
    "p"
  );


  legend->Draw();


  // ====================================================
  //                     Lower panel
  // ====================================================

  drawer->SetLogY(false);


  hEnergyEnergyCorrelatorRatio[0]
    ->GetXaxis()
    ->SetRangeUser(
      drawingRange.first,
      drawingRange.second
    );

  hEnergyEnergyCorrelatorRatio[0]
    ->GetYaxis()
    ->SetRangeUser(
      ratioZoom.first,
      ratioZoom.second
    );


  drawer->SetGridY(false);


  // First draw ppRef / ppRef = 1
  drawer->DrawHistogramToLowerPad(
    hEnergyEnergyCorrelatorRatio[0],
    "#Deltar",
    "#frac{OO}{pp}",
    " "
  );


  // Reference horizontal lines
  const double referenceLines[3] = {
    0.9,
    1.0,
    1.1
  };


  for(int iLine = 0;
      iLine < 3;
      iLine++){

    TLine* line =
      new TLine(
        drawingRange.first,
        referenceLines[iLine],
        drawingRange.second,
        referenceLines[iLine]
      );

    line->SetLineStyle(3);
    line->Draw("same");
  }


  // Red OO / ppRef points
  hEnergyEnergyCorrelatorRatio[1]
    ->Draw("same");


  // Redraw black baseline points on top of lines
  hEnergyEnergyCorrelatorRatio[0]
    ->Draw("same");


  // ====================================================
  //                      Save figure
  // ====================================================

  if(saveFigures){

    gPad->GetCanvas()->SaveAs(
      Form(
        "figures/"
        "OO%.0fto%.0f_over_ppRef_"
        "jetPt%.0fto%.0f_"
        "trackPt%.1f_raw.%s",

        comparedCentralityBin[1].first,
        comparedCentralityBin[1].second,

        comparedJetPtBin[0].first,
        comparedJetPtBin[0].second,

        comparedTrackPtBin[0],

        figureFormat
      )
    );
  }
}