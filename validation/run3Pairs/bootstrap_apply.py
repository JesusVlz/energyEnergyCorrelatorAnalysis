#!/usr/bin/env python3
"""One-time, atomic preparation patch for the inspected 7225afff source revision.

Used by the repair-branch build. All expected source blobs and replacement
anchors are checked before any existing source file is written. After CI has
committed the repaired files this script is an idempotent no-op.
"""
import hashlib
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MARKER = '// RUN3_TRACK_PREPARATION_V1'
PINS = {
 'src/EECAnalyzer.cxx': '375471c5664b693b0e97c7eb5408997c2f0074d9',
 'src/HighForestReader.cxx': 'c1044c144270d977f52f3376765d9aeb45e45e4a',
 'src/HighForestReader.h': '614e217d737b224edce79aea44eb6c001c856ba8',
 'src/trackingEfficiency2025OO.cxx': 'b58a01bc71e2118354979381ed9dc890fa085bc5',
 'src/trackingEfficiency2024pp.cxx': '15e11ef22e3ac182abe776cf771e370c717d53ba',
 'validation/run3Pairs/export_forest.py': '57102069dfda35049ada436e0e5c142dd967bd11',
 'validation/run3Pairs/pair_counts.py': '0b457bfa14c15701a7df0e6e439aead0c2d029cc',
 'validation/run3Pairs/plot_closure.py': '90cf2c4bc2582c9aa15a0e4af6af3ea41b2e3973',
 'validation/run3Pairs/batch_pairs.py': 'b901d0cce729d19dc49dec74e7cd74dbbb5d5847',
 'Makefile': '50c3b5c836050bb38ead086e6853ec76c2df3d73',
}


def apply():
    if MARKER in (ROOT/'src/EECAnalyzer.cxx').read_text():
        print('Run-3 preparation already applied; existing files unchanged.')
        return
    texts = {}
    for name, expected in PINS.items():
        data = (ROOT/name).read_bytes()
        actual = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
        if actual != expected:
            raise RuntimeError('Source changed; do not apply blind patch: ' + name + ' ' + actual)
        texts[name] = data.decode()

    def replace(name, before, after, count=1):
        found = texts[name].count(before)
        if found != count:
            raise RuntimeError(f'{name}: expected {count} exact anchors, found {found}: {before[:100]}')
        texts[name] = texts[name].replace(before, after)

    a = 'src/EECAnalyzer.cxx'
    replace(a, '#include "Run3JetGeometry.h"', '#include "Run3JetGeometry.h"\n#include "Run3TrackSelection.h"\n' + MARKER)
    start = texts[a].index('    // 1. Apply the nominal 2025 OO tracking correction.')
    end = texts[a].index('    // Derived from eecAnalysis_OOMC_eventWeights', start)
    texts[a] = texts[a][:start] + '''    // Match the Run-3 selection to its OO table (never a Run-2 payload).
    fTrackEfficiencyCorrector = new TrkEff2025OO(false,
      EECRun3Tracking::Payload(true, fTrackSelectionVariation));

''' + texts[a][end:]
    replace(a, 'new TrkEff2024ppRef(false, "trackCorrectionTables/pp2024/Eff_ppref_2024_Pythia_QCD_pThat15_NopU_2D_Nominal_Official.root")',
            'new TrkEff2024ppRef(false, EECRun3Tracking::Payload(false, fTrackSelectionVariation))')
    begin = texts[a].index('  // Determine the correct efficiency table to use based on track selection cuts')
    end = texts[a].index('  // Select the jet pT weighting configuration', begin)
    old = texts[a][begin:end]
    texts[a] = texts[a][:begin] + '''  if(fIsOOData || fIsPp2024Data){
    const std::string requested = fCard->HasKey("Run3TrackSelection")
      ? std::string(fCard->GetStr("Run3TrackSelection").Data()) : std::string();
    const std::string selection = EECRun3Tracking::ResolveSelection(
      requested, fMaxTrackPtRelativeError, fMaxTrackDistanceToVertex);
    fTrackSelectionVariation = selection == "Nominal" ? "" : selection;
    cout << "Run-3 track selection: " << selection << ", " << EECRun3Tracking::Version()
         << "; pT>10: sigmaPt/pT <= 0.1; vertex significance <= " << fMaxTrackDistanceToVertex << endl;
  } else {
''' + old + '  }\n\n' + texts[a][end:]
    line = '  if(trackReader->GetTrackPtError(iTrack)/trackPt >= fMaxTrackPtRelativeError) return false; // Cut for track pT relative error'
    replace(a, line, '''  if(fIsOOData || fIsPp2024Data){
    if(!std::isfinite(trackPt) || !std::isfinite(trackEta) ||
       !EECRun3Tracking::PassMomentum(trackPt, trackReader->GetTrackPtError(iTrack))) return false;
  } else {
''' + line + '\n  }', 2)
    for suffix, label in [('Z', 'z-direction'), ('XY', 'xy-direction')]:
        # Replace the entire exact comparison, leaving the historical comments.
        expression = f'if(TMath::Abs(trackReader->GetTrackVertexDistance{suffix}(iTrack)/trackReader->GetTrackVertexDistance{suffix}Error(iTrack)) >= fMaxTrackDistanceToVertex) return false;'
        new = f'''if(fIsOOData || fIsPp2024Data){{
    if(!EECRun3Tracking::PassVertex(trackReader->GetTrackVertexDistance{suffix}(iTrack),
       trackReader->GetTrackVertexDistance{suffix}Error(iTrack), fMaxTrackDistanceToVertex)) return false;
  }} else {{
    {expression}
  }}'''
        replace(a, expression, new, 2)
    anchor = 'if(fCard->HasKey("MinPairTrackPt") && lowerTrackPt <= fCard->Get("MinPairTrackPt")) continue;'
    replace(a, anchor, anchor + '''
        if((fIsOOData || fIsPp2024Data) &&
           (!std::isfinite(trackEfficiencyCorrection1) || trackEfficiencyCorrection1<=0 ||
            !std::isfinite(trackEfficiencyCorrection2) || trackEfficiencyCorrection2<=0))
          throw std::runtime_error("Unsupported single-track correction in a selected Run-3 EEC pair; no silent zero/unity replacement");''')
    anchor = '    if(!fTrackPairEfficiencyCorrector->IsRun3Table())\n      throw std::invalid_argument("Run3TrackPairCorrectionFile is missing Run3PairCorrectionVersion");'
    replace(a, anchor, anchor + '''
    // Old pre-repair maps cannot be applied to the corrected Run-3 selections.
    TFile metadata(fCard->GetStr("Run3TrackPairCorrectionFile"), "READ");
    const auto requireTag = [&metadata](const char* key, const std::string& expected){
      TObject* object = metadata.Get(key);
      if(!object || expected != object->GetTitle())
        throw std::invalid_argument(std::string("Run-3 pair payload mismatch/missing tag: ") + key);
    };
    requireTag("Run3TrackSelectionVersion", EECRun3Tracking::Version());
    requireTag("Run3TrackSelection", EECRun3Tracking::CanonicalSelection(fTrackSelectionVariation));
    requireTag("Run3PairCollisionSystem", fIsOOData ? "OO" : "pp");''')

    h = 'src/HighForestReader.h'
    replace(h, '#include "ForestReader.h"', '#include "ForestReader.h"\n#include "Run3TrackCharge.h"')
    replace(h, '  // Leaves for the track tree in AOD forests', '  EECRun3Tracking::TrackCharge fRun3TrackCharge;\n\n  // Leaves for the track tree in AOD forests')
    c = 'src/HighForestReader.cxx'
    run3 = '(fDataType == kOO || fDataType == kOOMC || fDataType == kPpRef5p36TeV || fDataType == kPpRef5p36TeVMC)'
    replace(c, '    fTrackTree->SetBranchStatus("*", 0);',
            '    fTrackTree->SetBranchStatus("*", 0);\n    if(' + run3 + ' && !fMegaSkimMode) fRun3TrackCharge.Bind(fTrackTree);')
    replace(c, 'Int_t HighForestReader::GetTrackCharge(Int_t iTrack) const{\n  return 1;\n}',
            'Int_t HighForestReader::GetTrackCharge(Int_t iTrack) const{\n  if(' + run3 + ' && !fMegaSkimMode) return fRun3TrackCharge.Get(iTrack);\n  return 1; // Unchanged legacy/mega-skim path, not used for Run-3 matching.\n}')
    replace(c, '  fRun3PFJetInputs = in.fRun3PFJetInputs;', '  fRun3PFJetInputs = in.fRun3PFJetInputs;\n  fRun3TrackCharge = in.fRun3TrackCharge;', 2)

    oo = 'src/trackingEfficiency2025OO.cxx'
    pp = 'src/trackingEfficiency2024pp.cxx'
    for name in (oo, pp):
        texts[name] = '#include "TrackingTableValidation.h"\n' + texts[name]
    replace(oo, '      TMath::Abs(eta) <= 2.4 &&', '      std::isfinite(pt) && std::isfinite(eta) &&\n      TMath::Abs(eta) <= 2.4 &&')
    anchor = '    (pt > kMaximumLookupPt) ? kMaximumLookupPt : pt;'
    replace(oo, anchor, anchor + '''
  if(!EECRun3Tracking::InsideAxis(efficiency->GetXaxis(), eta) ||
     !EECRun3Tracking::InsideAxis(efficiency->GetYaxis(), lookupPt) ||
     !EECRun3Tracking::InsideAxis(efficiency->GetZaxis(), hiBin)) return 0;''')
    replace(oo, '  if (efficiencyValue <= 0.001) {', '''  if (!EECRun3Tracking::ValidContents(efficiencyValue, fakeRateValue, secondaryRateValue, multipleRateValue) || efficiencyValue <= 0.001) {''')
    anchor = '  if (!isQuiet) {\n\n    std::cout\n      << "TrkEff2025OO: using nominal OO tracking correction table: "'
    replace(oo, anchor, '''  EECRun3Tracking::ValidateTables({efficiency, fakeRate, secondaryRate, multipleRate}, 3, fileName);

  if (!isQuiet) {

    std::cout
      << "TrkEff2025OO: using selected OO tracking correction table: "''')
    replace(pp, '\n}\n\nTrkEff2024ppRef::~TrkEff2024ppRef()', '\n  EECRun3Tracking::ValidateTables({fEfficiency, fFakeRate, fSecondaryRate, fMultipleRecoRate}, 2, filePath);\n}\n\nTrkEff2024ppRef::~TrkEff2024ppRef()')
    replace(pp, '  if(fEfficiency == NULL){', '  if(fEfficiency == NULL || !std::isfinite(pt) || !std::isfinite(eta)){')
    replace(pp, '  const float denominator =', '''  if(!EECRun3Tracking::ValidContents(efficiency, fakeRate, secondaryRate, multipleRecoRate)) return 0;

  const float denominator =''')

    e = 'validation/run3Pairs/export_forest.py'
    replace(e, 'from pathlib import Path', 'from pathlib import Path\nfrom preparation_contract import selection_from_card, payload_path, export_contract')
    replace(e, "    payload = PAYLOADS[args.system]", '''    base = args.card or REPO/('cardEECOOMCInclusiveR04.input' if args.system == 'OO' else 'cardEECpp5v36TeVMCR04.input')
    selection = selection_from_card(base)
    payload = payload_path(args.system, selection)
    contract = export_contract(args.system, selection, payload)''')
    replace(e, "'MinJetPtCut': 40, 'MaxJetPtCut': 160, 'MinTrackPtCut': 2,", "'MinJetPtCut': 40, 'MaxJetPtCut': 160, 'MinTrackPtCut': 0.5,")
    replace(e, "'FilledHistograms': 275 if args.responses and population == 'reco' else 19,", "'FilledHistograms': 279 if args.responses and population == 'reco' else 23,")
    # Reusing unversioned native rows could relabel an old selection as new.
    replace(e, '    paths, commands = {}, []', '''    if args.reuse_native:
        raise ValueError('Native reuse requires a verified campaign receipt; regenerate this standalone export after the Run-3 cut repair')
    paths, commands = {}, []''')
    replace(e, "            result.write(json.dumps(reco, allow_nan=False)+'\\n')", "            reco['tracking_contract'] = contract\n            result.write(json.dumps(reco, allow_nan=False)+'\\n')")
    replace(e, "    metadata = {'status': 'NATIVE_SELECTION_EXPORT', 'system': args.system,", "    metadata = {'status': 'NATIVE_SELECTION_EXPORT', 'system': args.system, 'tracking_contract': contract,")
    replace(e, "'analyzer_sha256': digest(REPO/'src/EECAnalyzer.cxx'),", "'analyzer_sha256': digest(REPO/'src/EECAnalyzer.cxx'), 'native_binary_sha256': digest(REPO/'eecAnalysis'),")
    # Guard the auxiliary old-coordinate diagnostic without treating it as usable.
    replace(e, "                    relative = abs(track['correction']/track['legacy_relative_eta_correction']-1)", "                    old_factor = track['legacy_relative_eta_correction']\n                    relative = abs(track['correction']/old_factor-1) if math.isfinite(old_factor) and old_factor > 0 else 1.")

    p = 'validation/run3Pairs/pair_counts.py'
    replace(p, "VERSION = 'run3-an751-v1'", "VERSION = 'run3-an751-v2'")
    replace(p, '    events = JsonlEvents(paths)', '    from preparation_contract import check_export_metadata\n    tracking_contracts = check_export_metadata(paths)\n    events = JsonlEvents(paths)')
    replace(p, '    result = run(events, args.min_pairs, tuple(args.control))', '    result = run(events, args.min_pairs, tuple(args.control))\n    result[\'tracking_contracts\'] = tracking_contracts')
    # The prior generator function did not execute its identity check unless
    # consumed: validate duplicate track indices explicitly.
    replace(p, '    pairs(tracks)\n', "    indices = [track['index'] for track in tracks]\n    if len(indices) != len(set(indices)):\n        raise ValueError('Duplicate physical track index within one population')\n")

    plot = 'validation/run3Pairs/plot_closure.py'
    replace(plot, 'from pathlib import Path', 'from pathlib import Path\nfrom preparation_contract import write_contract')
    replace(plot, "ROOT.TNamed('Run3PairCorrectionVersion', 'run3-an751-v1').Write()", "ROOT.TNamed('Run3PairCorrectionVersion', result['correction_version']).Write()\n        write_contract(ROOT, result, rootfile, system)")
    # Corrected tables must never overwrite the historical pre-repair version.
    replace(plot, "'_run3-an751-v1.root'", "'_' + result['correction_version'] + '.root'")

    batch = 'validation/run3Pairs/batch_pairs.py'
    replace(batch, 'from pathlib import Path', 'from pathlib import Path\nfrom preparation_contract import validate_contract')
    replace(batch, "'payload_sha256','tool_sha256','event_weight_policy','jet_policy'", "'payload_sha256','tool_sha256','event_weight_policy','jet_policy','tracking_contract'")
    anchor = "    metadata=read_json(metadata_path) if metadata_path else {}"
    replace(batch, anchor, anchor + '''
    if metadata.get('tracking_contract') is not None:
        validate_contract(metadata['tracking_contract'])
    elif task and task.get('preparation_version'):
        raise ValueError('New campaign requires corrected, versioned native export')''')
    anchor = "                  batch_status='GLOBAL_COUNTS_NORMALIZED_HOLDOUT_NOT_YET_EVALUATED')"
    replace(batch, anchor, anchor + '''
    if native_signature and native_signature.get('tracking_contract'):
        contract = validate_contract(native_signature['tracking_contract'])
        result['tracking_contracts'] = {system: contract}''')

    texts['Makefile'] += '''
# Run-3 tracking changes must invalidate dependent objects in incremental builds.
src/EECAnalyzer.o: src/Run3TrackSelection.h
src/HighForestReader.o: src/Run3TrackCharge.h
src/trackingEfficiency2025OO.o src/trackingEfficiency2024pp.o: src/TrackingTableValidation.h
'''
    # Leave Run-2 and R8 cards alone. Preserve jet/track thresholds and response
    # guards; pair preparation has its own bounded 40-160 GeV cards.
    for path in ROOT.glob('card*.input'):
        text = path.read_text()
        if not re.search(r'^DataType\s+(9|10|11|12)\b',text,re.M) or not re.search(r'^JetRadius\s+0\.4\b',text,re.M):
            continue
        text = re.sub(r'^MaxTrackPtRelativeError.*$', 'MaxTrackPtRelativeError 0.1 # Run-3: <=0.1 only above 10 GeV, all selections', text, flags=re.M)
        if 'Run3TrackSelection ' not in text:
            text = text.replace('# Cuts for tracks', '# Cuts for tracks\nRun3TrackSelection Nominal # Change together with VertexMaxDistance and matching payload')
        text = re.sub(r'^HFEnergyShift.*$', 'HFEnergyShift 0 # OFF for Run-3; no inherited HF shift',text,flags=re.M)
        texts[str(path.relative_to(ROOT))] = text
    for name, text in texts.items():
        (ROOT/name).write_text(text)
    print('Applied checked changes to', len(texts), 'source/card files; no MC has run.')
    from prepare_inputs import write_cards
    write_cards(ROOT/'validation/run3Pairs/cards')


if __name__ == '__main__':
    apply()
