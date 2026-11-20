#!/usr/bin/env python3
"""Generate explicit preliminary pp R4 cards and relocate deferred R8 cards."""
import json
from build_r4_catalog import REPO, defer_r8, parse

def main():
    manifest=json.loads((REPO/'trackCorrectionTables/pp2024/pair/20261116/manifest.json').read_text())
    out=REPO/'cards/ppRef2024/R4';out.mkdir(parents=True,exist_ok=True)
    for sample,base in [('Data','cardEECpp5v36TeVR04.input'),('MC','cardEECpp5v36TeVMCR04.input')]:
        for wp,vertex in [('Nominal',3),('Loose',5),('Tight',2)]:
            c=parse((REPO/base).read_text())
            c.update(Run3TrackSelection=wp,VertexMaxDistance=str(vertex),MaxTrackPtRelativeError='0.1',
                Run3TrackPairCorrectionFile=manifest['tables'][wp]['path'],DisableTrackPairEfficiencyCorrection='0',
                Run3JetID='1',Run3UsePpL1FastJet='0',MinPairTrackPt='2',MinTrackPtCut='0.5',
                McCorrelationType='0',FilledHistograms='31',BackgroundMethods='0',
                TrackEfficiencyVariation='0',MaxEventsPerFile='1000',MinJetPtCut='15',MaxJetPtCut='200',
                JetPtBinEdgesEEC='15 20 30 40 80 120 160 200',FirstUnfoldedJetPtBin='0',LastUnfoldedJetPtBin='6',
                CentralityBinEdges='-1 100',DebugLevel='1')
            c.pop('MixingFileList',None)
            text=(f'# PRELIMINARY ppRef2024 R4 {sample} {wp}: supplied 20261116 residual pair table.\n'
                '# Requires checksum-checked ROOT installation; see README.md and manifest.json.\n'
                '# Jet ID ON by request; original table was derived with Jet ID OFF.\n'
                '# No new closure is claimed. Ordinary EEC only; no particle/unfolding responses.\n'
                '# 15-200 GeV application test: 15-20 uses the stored 10-20 correction bin.\n'
                '# Run from repository root. Trigger/sample applicability remains a user input.\n')
            text+='\n'.join(k+' '+v for k,v in c.items())+'\n'
            (out/f'ppRef_{sample}_{wp}_pairOn.input').write_text(text)
    moved=defer_r8()
    records=json.loads((REPO/'cards/R8/moved.json').read_text())
    p=REPO/'cards/README.md';text=p.read_text()
    for old,row in records.items():
        text=text.replace('](../'+old+')','](../'+row['path']+')')
    if '## ppRef R4 preliminary pair tables' not in text:
        text+=('\n## ppRef R4 preliminary pair tables\n\n'
            'See [ppRef2024/R4/README.md](ppRef2024/R4/README.md) for the checksum-checked\n'
            '16 November 2026 Nominal/Loose/Tight candidates and JetID-on application cards.\n'
            'The ROOT binaries require the separate verified installation; cards alone do not install them.\n'
            'These are preliminary tests, not qualified final-spectrum defaults.\n'
            'R8 cards have moved to [R8/pending/](R8/pending/); R4 remains the active scope.\n')
    p.write_text(text)
    p=REPO/'README.md';text=p.read_text().replace(
        'The repository root retains the 12 OO2025 and ppRef2024 5.36 TeV cards,\nincluding the existing separate R=0.8 workstream.',
        'The repository root retains the eight R4 OO2025 and ppRef2024 5.36 TeV cards.\nDeferred R8 cards are under [cards/R8/pending/](cards/R8/pending/).')
    p.write_text(text)
    print(json.dumps({'generated_pp_cards':6,'relocated_r8_cards':len(moved),
        'root_payloads':'SEPARATE_CHECKSUM_VERIFIED_INSTALLATION_REQUIRED',
        'production_submitted':False,'physics_closure':'PENDING_JETID_ON'}))

if __name__=='__main__':
    main()
