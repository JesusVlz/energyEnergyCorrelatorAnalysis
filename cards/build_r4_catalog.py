#!/usr/bin/env python3
"""Prepare configurations only. Never submit jobs or invent calibration inputs.

Run from the repository root. --write generates the catalog; --check verifies it.
R8 relocation is explicit and preserves original card bytes and provenance.
"""
import argparse
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
VERSION = 'run3-r4-cards-v2-slim-tracks'
UNUSED_TRACK_KEYS = {'CalorimeterSignalLimitPt', 'HighPtEtFraction', 'Chi2QualityCut', 'MinimumTrackHits'}
BASE_SHA = '6d4800b8c57747e1db761fc4fecf72a38a1a00c3'
EDGES = '5 10 20 40 80 120 160 200 300 500 5020'
VERTEX = {'Nominal': 3, 'Loose': 5, 'Tight': 2}
MODES = {'RecoReco': 0, 'RecoGen': 1, 'GenReco': 2, 'GenGen': 3}
CONES = {'none': 0, 'reflected': 1, 'perpendicular': 2, 'sameEventCones': 3,
         'mixed': 4, 'mixedPerpendicular': 6, 'allCones': 7}
TOKEN = re.compile(r'@([A-Z][A-Z0-9_]*)@')


def sha(data):
    return hashlib.sha256(data).hexdigest()


def parse(text):
    result = {}
    for line in text.splitlines():
        words = line.split('#', 1)[0].split()
        if not words:
            continue
        if len(words) < 2 or words[0] in result:
            raise ValueError('Invalid/duplicate configuration key: ' + words[0])
        result[words[0]] = ' '.join(words[1:])
    return result


def values(path):
    return parse((REPO / path).read_text())


def base(system, mc, activity):
    name = ('cardEECOOMCInclusiveR04.input' if mc else 'cardEECOOInclusiveR04.input') if system == 'OO' else (
        'cardEECpp5v36TeVMCR04.input' if mc else 'cardEECpp5v36TeVR04.input')
    result = values(name)
    for key in list(result):
        if key in UNUSED_TRACK_KEYS or key.startswith(('PairExport', 'OOClosure')) or key in ('MixingFileList', 'Run3TrackPairCorrectionFile'):
            del result[key]
    result.update({k: str(v) for k, v in {
        'DataType': (10 if mc else 9) if system == 'OO' else (12 if mc else 11),
        'McCorrelationType': 0, 'MatchJets': 0, 'JetType': 1, 'JetRadius': .4,
        'ChannelRadiusCode': 4, 'JetAxis': 1, 'JetEtaCut': 1.6,
        'MinJetPtCut': 40, 'MaxJetPtCut': 5020, 'JetPtBinEdgesEEC': EDGES,
        'FirstUnfoldedJetPtBin': 3, 'LastUnfoldedJetPtBin': 5,
        'TrackEtaCut': 2.4, 'MinTrackPtCut': .5, 'MinPairTrackPt': 2,
        'MaxTrackPtCut': 500, 'WeightExponent': '1 2', 'SubeventCut': 2,
        'Run3TrackSelection': 'Nominal', 'MaxTrackPtRelativeError': .1,
        'VertexMaxDistance': 3, 'DisableTrackPairEfficiencyCorrection': 1,
        'TrackEfficiencyVariation': 0, 'JetUncertainty': 0, 'JetPtWeight': 0,
        'SmearDeltaR': 0, 'SmearEnergyWeight': 0, 'BackgroundMethods': 0,
        'MegaSkimMixing': 0, 'BinnedMixing': 0, 'EventMatchForMixing': 0,
        'HFEnergyShift': 0, 'MinJetPtUnfoldingReco': 40, 'MinJetPtUnfoldingTruth': 5,
        'FilledHistograms': 31, 'SkipCovarianceMatrix': 0, 'MaxEventsPerFile': 0,
        'Run3JEC': 1, 'Run3JetID': 1, 'Run3JetIDTable': 'MiniPOG20261108',
        'Run3JetQA': 1, 'ApplyJetVeto': 1, 'DebugLevel': 0,
        'CentralityBinEdges': '-1 0 30 50 90 100' if system == 'OO' and activity == 'split' else '-1 100'
    }.items()})
    # Preserve the existing data-trigger choice, including its applicability caveat.
    if mc:
        result['TriggerSelection'] = '0'
    if system == 'pp':
        result['Run3UsePpL1FastJet'] = '0'
    return result


def validate(card):
    if UNUSED_TRACK_KEYS.intersection(card):
        raise ValueError('Run-3 cards must not request PbPb-only track cuts')
    def number(key):
        return float(card[key])
    if 'InclusiveMode' in card:
        raise ValueError('InclusiveMode is retired')
    for key, expected in {'JetType': 1, 'JetRadius': .4, 'ChannelRadiusCode': 4,
                          'JetAxis': 1, 'JetEtaCut': 1.6, 'TrackEtaCut': 2.4,
                          'MinPairTrackPt': 2, 'JetUncertainty': 0, 'JetPtWeight': 0,
                          'SmearDeltaR': 0, 'SmearEnergyWeight': 0, 'HFEnergyShift': 0}.items():
        if number(key) != expected:
            raise ValueError('Incompatible R4 catalog value: ' + key)
    if card['WeightExponent'].split() != ['1', '2']:
        raise ValueError('Both powers required')
    selection = card['Run3TrackSelection']
    if number('VertexMaxDistance') != VERTEX[selection] or number('MaxTrackPtRelativeError') != .1:
        raise ValueError('Track selection/table mismatch')
    jet_edges = list(map(float, card['JetPtBinEdgesEEC'].split()))
    first, last = int(number('FirstUnfoldedJetPtBin')), int(number('LastUnfoldedJetPtBin'))
    if jet_edges[first] != 40 or jet_edges[last + 1] != 160:
        raise ValueError('Wrong report-bin indices')
    if any(b <= a for a, b in zip(jet_edges, jet_edges[1:])):
        raise ValueError('Jet edges must increase')
    if len(jet_edges) > 11 or len(card['CentralityBinEdges'].split()) > 6:
        raise ValueError('Native covariance array capacity exceeded')
    if int(number('BackgroundMethods')) & 4 and 'MixingFileList' not in card:
        raise ValueError('Mixing requires an explicit pool')
    if number('DisableTrackPairEfficiencyCorrection') == 0 and 'Run3TrackPairCorrectionFile' not in card:
        raise ValueError('Pair-corrected template requires explicit payload')
    bits = int(number('FilledHistograms'))
    if int(number('DataType')) in (9, 11) and bits & (64 | 128 | 256):
        raise ValueError('Data cannot fill truth closure/response objects')
    if bits & 128:
        if number('McCorrelationType') != 1 or number('MinTrackPtCut') != 2 or number('SkipCovarianceMatrix') != 0:
            raise ValueError('Generator-particle response contract changed')
        if number('BackgroundMethods') != 0 or number('DisableTrackPairEfficiencyCorrection') != 1:
            raise ValueError('Responses must not inherit pair correction or mixing')
    if card.get('OOClosurePilot') == '1' and (number('BackgroundMethods') != 6 or number('McCorrelationType') != 3 or number('EventMatchForMixing') != 2):
        raise ValueError('Strict OO benchmark is mixed+perpendicular GenGen/NH only')
    return True


def catalog():
    files, records = {}, []
    def add(group, name, card, purpose, notes='', requirements=()):
        card = {k: str(v) for k, v in card.items()}
        validate(card)
        tokens = sorted(set(TOKEN.findall('\n'.join(card.values()))))
        suffix = '.input.in' if tokens else '.input'
        path = f'cards/R4/{group}/{name}{suffix}'
        text = (f'# {VERSION}: {purpose}\n'
                '# Preparation, not a production submission or a physics validation.\n'
                '# Run from repository root; all payload paths are repository-relative.\n'
                '# Report 40-80,80-120,120-160 GeV; both n=1,2; strict pair pT>2.\n'
                '# Stored activity is not an event cut or percentile certification.\n')
        if notes:
            text += '# ' + notes.replace('\n', '\n# ') + '\n'
        if tokens:
            text += '# TEMPLATE ONLY: resolve ' + ', '.join(tokens) + ' before use.\n'
        text += '\n'.join(k + ' ' + v for k, v in card.items()) + '\n'
        if path in files:
            raise ValueError('Duplicate catalog path')
        files[path] = text
        records.append({'path': path, 'purpose': purpose,
                        'status': 'TEMPLATE_REQUIRES_INPUTS' if tokens else 'CONFIGURED_NOT_PHYSICS_VALIDATED',
                        'required_tokens': tokens, 'requirements': list(requirements),
                        'sha256': sha(text.encode()), 'DataType': int(card['DataType']),
                        'correlation': int(card['McCorrelationType']),
                        'background_mask': int(card['BackgroundMethods']),
                        'pair_correction': card['DisableTrackPairEfficiencyCorrection'] == '0',
                        'selection': card['Run3TrackSelection']})
    for system in ('OO', 'pp'):
        scopes = ('integrated', 'split') if system == 'OO' else ('inclusive',)
        for activity in scopes:
            for mc in (False, True):
                label = ('OO' if system == 'OO' else 'ppRef') + ('_MC_' if mc else '_Data_') + activity
                c = base(system, mc, activity)
                for cone, mask in CONES.items():
                    cc = dict(c, BackgroundMethods=mask)
                    req = ['compatible native R4 forest; nominal jet/tracking/weight applicability']
                    if mask & 4:
                        cc['MixingFileList'] = '@MIXING_LIST@'
                        req += ['independent compatible donor pool and matching validation']
                    add('measurement', label + '_' + cone, cc,
                        'Single-track-corrected EEC and cone components; residual pair correction OFF.',
                        'No subtraction happens in the analyzer. pp nominal no-subtraction remains a downstream choice.\n'
                        'Ordinary matcher is the retained legacy path, not the strict OO NH pilot.\n'
                        'Data trigger is inherited from its current R4 base; check actual sample/trigger acceptance.', req)
                for selection in ('Loose', 'Tight'):
                    cc = dict(c, Run3TrackSelection=selection, VertexMaxDistance=VERTEX[selection], BackgroundMethods=3)
                    add('systematics/track-selection', label + '_' + selection, cc,
                        'Track-selection study with its corresponding existing single-track table.',
                        'Pair correction is OFF. Final variation also needs a compatible selection-specific pair map.')
                for selection in VERTEX:
                    for cone, mask in (('none', 0), ('allCones', 7)):
                        cc = dict(c, Run3TrackSelection=selection, VertexMaxDistance=VERTEX[selection],
                                  DisableTrackPairEfficiencyCorrection=0,
                                  Run3TrackPairCorrectionFile='@PAIR_PAYLOAD@', BackgroundMethods=mask)
                        if mask & 4:
                            cc['MixingFileList'] = '@MIXING_LIST@'
                        add('corrected-templates', label + '_' + selection + '_' + cone, cc,
                            'Candidate tracking+pair-corrected spectrum for downstream unfolding/subtraction.',
                            'Require system/selection-compatible independently qualified map over ALL filled jets and tracks.\n'
                            'A map limited to 40<=pTjet<160 cannot serve the wider migration-support spectrum here.\n'
                            'Do not narrow response guards merely to fit the map. Same-source pair topology is unchanged.',
                            ['qualified current selection-specific pair map including migration support',
                             'compatible data/MC response and downstream background adjustment'])
                cc = dict(c, FilledHistograms=63, BackgroundMethods=7,
                          MixingFileList='@MIXING_LIST@', DisableTrackPairEfficiencyCorrection=0,
                          Run3TrackPairCorrectionFile='@PAIR_PAYLOAD@',
                          TrackEfficiencyVariation='@TRACKING_RELATIVE_VARIATION@')
                add('systematics/templates', label + '_trackingFactors', cc,
                    'Fill existing nominal/no-efficiency/up/down EEC variation objects.',
                    'No Run-2 percentage is assumed. Table errors are NOT a complete pair systematic.\n'
                    'Use separately qualified varied pair payloads for pair-model and finite-map studies.',
                    ['Run-3 tracking variation prescription; covariance and data propagation'])
                if not mc:
                    continue
                for mode in ('RecoGen', 'GenReco', 'GenGen'):
                    add('mc-closure', label + '_' + mode,
                        dict(c, McCorrelationType=MODES[mode]),
                        mode + ' diagnostic on the standard R4 path; compare appropriate populations.',
                        'RecoGen holds reconstructed jets; GenGen changes jet population/axis. Not pure tracking closure.')
                for mode in ('RecoReco', 'RecoGen'):
                    add('mc-closure', label + '_' + mode + '_matchedJets',
                        dict(c, McCorrelationType=MODES[mode], MatchJets=1),
                        'Matched-jet-only diagnostic, not the nominal inclusive measurement.',
                        'Run both reconstructed-jet modes on identical events and preserve unmatched exposure.')
                add('jet-qa', label + '_jetResponse', dict(c, FilledHistograms=67),
                    'Jet QA and reconstructed/generated jet response; no EEC cone filling.')
                add('particle-response', label + '_trackPairResolution', dict(c, FilledHistograms=279),
                    'Matched-track angle and momentum-product diagnostics, separate from jet-pT unfolding.',
                    'No dedicated Run-3 smearing model is enabled by this response output.')
                for guard in (5, 10):
                    cc = dict(c, McCorrelationType=1, MinTrackPtCut=2, FilledHistograms=211,
                              MinJetPtUnfoldingTruth=guard, SkipCovarianceMatrix=0)
                    add('response', label + '_truthGuard' + str(guard), cc,
                        'Native generator-particle jet-pT response, misses/fakes and per-jet covariance.',
                        'Apply existing synchronized physical-event folds externally; filenames do not create independent halves.\n'
                        'truthGuard5 reproduces the current support; truthGuard10 is a sensitivity test, not an approved replacement.\n'
                        'Response production is not unfolding or a selected iteration count.')
    # Preserve the preparation cards' existing exact prescription, including their bounded cap.
    for path in sorted((REPO/'validation/run3Pairs/cards').glob('*.input')):
        cc = parse(path.read_text())
        cc['FirstUnfoldedJetPtBin'], cc['LastUnfoldedJetPtBin'] = '0', '2'
        cc['Run3JEC'], cc['Run3JetIDTable'] = '1', 'MiniPOG20261108'
        cc.pop('MixingFileList', None)
        add('pairs', path.stem, cc,
            'Bounded existing pair derivation card; use prepare_inputs.py plan/export/batch, not EEC division.',
            '10000 entries per file here is a pilot cap. Derive unweighted momentum pair counts;\n'
            'merge training counts before large-angle normalization and retain independent validation.\n'
            'This 40-160 map domain is NOT automatically sufficient for the full unfolding-support spectrum.')
    for activity in ('integrated', 'split'):
        cc = base('OO', True, activity)
        cc.update(McCorrelationType=3, BackgroundMethods=6, EventMatchForMixing=2,
                  OOClosurePilot=1, OOClosureMaxSignalEvents=10000, OOClosureMaxPoolEvents=100000,
                  OOClosureSeed=20261030, OOClosurePartitionPolicy=0,
                  OOClosureOutputPrefix='@LEDGER_PREFIX@', MixingFileList='@MIXING_LIST@')
        add('background-benchmark', 'OO_MC_' + activity + '_GenGen_NH', cc,
            'Strict two-donor OO GenGen NH/vz pilot, separate from ordinary all-cone configuration.',
            'Reflection is forbidden by the current strict-pilot implementation.\n'
            'Use the existing closure-pilot runner for overlap screening and full ledgers; do not relabel as full production.',
            ['pure-HIJING donor/source independence; unique writable ledger prefix'])
    manifests = {}
    for definition in ('response_marginal', 'ordinary_GenGen'):
        channels = []
        for system in ('OO', 'pp'):
            activities = [('integrated', [-1, 100]), ('C0_30', [0, 30]), ('C30_50', [30, 50]), ('C50_90', [50, 90])] if system == 'OO' else [('inclusive', [-1, 100])]
            for activity, bounds in activities:
                for power in (1, 2):
                    cid = f'{system}_{activity}_n{power}'
                    c = {'id': cid, 'system': system, 'activity': bounds, 'power': power,
                         'particle_pt': 2, 'truth_definition': definition,
                         'train': '@' + cid.upper() + '_TRAIN@',
                         'validation': '@' + cid.upper() + '_VALIDATION@',
                         'split_evidence': '@' + cid.upper() + '_SPLIT_EVIDENCE@'}
                    if definition == 'ordinary_GenGen':
                        for key in ('train_reco', 'train_gen', 'validation_reco', 'validation_gen'):
                            c[key] = '@' + cid.upper() + '_' + key.upper() + '@'
                    channels.append(c)
        path = 'cards/R4/unfolding/' + definition + '.json.in'
        manifests[path] = {'schema_version': 1, 'delta_r_region': [.008, .39], 'channels': channels}
    for path, data in manifests.items():
        files[path] = json.dumps(data, indent=2) + '\n'
    pending = {
        'schema_version': 1, 'status': 'NOT_ALL_SYSTEMATICS_ARE_EXECUTABLE',
        'items': [
            {'source': 'JES/JER', 'status': 'PENDING_INTERFACE_AND_INPUTS',
             'reason': 'Analyzer explicitly rejects nonzero Run-3 JetUncertainty; applicable inputs and implementation are required.'},
            {'source': 'jet-spectrum prior', 'status': 'PENDING_INTERFACE_AND_INPUTS',
             'reason': 'Run-3 JetPtWeight is guarded; use justified response/pseudo-data reweighting only after implementation/validation.'},
            {'source': 'DeltaR and momentum smearing', 'status': 'PENDING_RESPONSE_INPUTS',
             'reason': 'SmearDeltaR/SmearEnergyWeight cannot use the old Run-2 models.'},
            {'source': 'track efficiency', 'status': 'TEMPLATE_INPUT_REQUIRED',
             'reason': 'Use the applicable Run-3 fractional variation; 0.024 is not adopted by this catalog.'},
            {'source': 'track selection and residual pair maps', 'status': 'MC_DERIVATION_REQUIRED',
             'reason': 'Nominal/Loose/Tight single-track payloads exist; each corrected result needs compatible qualified pair maps.'},
            {'source': 'pair and response finite MC', 'status': 'DOWNSTREAM_RESAMPLING',
             'reason': 'Reuse independent counts/ledgers; repeated identical cards do not generate independent MC statistics.'},
            {'source': 'background method and level', 'status': 'COMPONENTS_THEN_DOWNSTREAM',
             'reason': '0/1/2/3/4/6/7 masks are available. Validate matching/donor reuse and post-unfolding adjustment. pp nominal remains unsubtracted.'},
            {'source': 'HF shift', 'status': 'OPTIONAL_NOT_NOMINAL',
             'reason': 'Keep zero with existing hiBin/vz ordinary matching or strict NH pilot. Zero alone does not certify an untransformed legacy HF-matching path.'},
            {'source': 'iteration/prior/model/refolding', 'status': 'UNFOLDING_STUDY',
             'reason': 'Use manifests and existing 1-20 two-direction scan; no iteration value selected here. Alternative detector MC must be supplied.'},
            {'source': 'trigger, calibration, weights, centrality', 'status': 'EXTERNAL_APPLICABILITY',
             'reason': 'Current base policy is retained. Validate source retention and data recommendations; no centrality recalibration or invented normalizations.'},
            {'source': 'full-chain residual and OO/pp covariance', 'status': 'AFTER_NOMINAL_QUALIFICATION',
             'reason': 'Apply variations to measured data too; component MC closure is not a completed measurement systematic.'}
        ]}
    files['cards/R4/systematics/pending.json'] = json.dumps(pending, indent=2) + '\n'
    inventory = {'version': VERSION, 'reviewed_base': BASE_SHA,
                 'scope': 'R4 only; R8 deferred and excluded from active catalog',
                 'source_cards': {p: sha((REPO/p).read_bytes()) for p in
                                  ('cardEECOOInclusiveR04.input', 'cardEECOOMCInclusiveR04.input',
                                   'cardEECpp5v36TeVR04.input', 'cardEECpp5v36TeVMCR04.input')},
                 'cards': records, 'unfolding_templates': list(manifests),
                 'physics_closure': 'NOT_RUN', 'production_submitted': False}
    files['cards/R4/catalog.json'] = json.dumps(inventory, indent=2) + '\n'
    files['cards/R4/README.md'] = documentation(records)
    return files


def documentation(records):
    ready = sum(not r['required_tokens'] for r in records)
    templated = len(records) - ready
    return f'''# R4 OO/ppRef configuration catalog

{ready} configured analyzer cards, {templated} input-dependent templates, and two unfolding-study manifest templates.
Generated by `python3 cards/build_r4_catalog.py --write`; verify with `--check`.
No data/MC jobs have been submitted. Configuration validity is not physics closure.

## Scope and minimal execution order

R4 WTA; jet |eta|<1.6; track |eta|<2.4; strict pair pT>2 GeV; n=1 and n=2.
Report 40-80, 80-120, 120-160 GeV. OO integrated selected MB and the three
0-30/30-50/50-90 bins; split storage also retains -1..0 and 90..100 bookkeeping.
There is no activity event cut or InclusiveMode. Validate any out-of-range flow.
The response keeps 5..5020 GeV support, reco >=40 GeV and report indices 3..5.
Do not confuse reporting bins with migration support.

| Folder | Use |
|---|---|
| measurement | Data/RecoReco MC, single-track correction ON, pair OFF, covariance filling requested |
| mc-closure | RecoGen/GenReco/GenGen comparisons and matched-only reconstructed-jet diagnostics |
| pairs | Copies of the existing 12 bounded Nominal/Loose/Tight RecoReco/RecoGen preparation cards |
| response | RecoGen generator-particle response, both powers; 5 GeV truth guard and named 10 GeV comparison |
| jet-qa | Jet closure/selection/calibration QA |
| particle-response | Track-pair angular/momentum response diagnostics |
| corrected-templates | Input-dependent Nominal/Loose/Tight pair-corrected data/MC, no-cone or all-cone |
| background-benchmark | Strict two-donor OO GenGen NH/vz closure, mixed+perpendicular only |
| systematics | Track-selection cards, factor-variation templates and explicit pending-source register |
| unfolding | Two-fold response-marginal and ordinary-GenGen truth manifest templates |

Start with native schema/jet QA and the nominal pair preparation. Derive maps
from unweighted momentum pair counts with the existing exporter/batch machinery,
not by dividing energy-weighted EEC histograms. Keep training and holdout disjoint.
Then produce compatible response and measured spectra, run the two-direction
1-20 iteration study, adjust/subtract background AFTER jet-pT unfolding,
normalize, and construct OO/pp. Unfolding is not performed by analyzer cards.

## Background masks: one all-cone run is usually enough

0 none; 1 reflected; 2 perpendicular; 3 reflected+perpendicular; 4 mixed;
6 mixed+perpendicular; 7 all. The single-mask cards are diagnostics, NOT an
instruction to process the full data seven times. All-cone components share
signal-event processing; check matched exposure before comparing their results.
The analyzer fills components; it does not perform final subtraction.
pp no-nominal-subtraction remains the reference convention.

Ordinary mixed cards retain EventMatchForMixing=0: stored hiBin/vz for OO,
vz-only for pp with the current constant-hiBin reader. Their legacy widening,
identity and counting limitations are NOT fixed by these cards. Supply a genuine
independent compatible MC/data donor pool, not a signal-QCD list renamed MB.
The strict OO NH pilot is separate, uses BackgroundMethods=6, fixed caps and
identity/background-overlap ledgers. Reflected cones are forbidden in that pilot.
HF shifting is optional and OFF in this suite. No guessed HF matching replacement.

## Inputs and unresolved corrections

Current R4 JEC/Jet ID/veto and event-weight policies are retained. pp L1FastJet
stays OFF; do not enable it on Vanderbilt forests missing rho/area inputs.
Data trigger is inherited from the dedicated base (currently selection 8,
OO MB / ppRef ZeroBias). Check trigger/sample acceptance; this is NOT a claim
that every HardProbes forest provides an unbiased usable reference.
All paths resolve from repository root. Existing root R4 and CRAB cards remain
compatibility entries with the same physics settings and obsolete track keys removed; select this catalog explicitly for new jobs.
Existing campaign tarballs were NOT rebuilt and their cards were NOT replaced.

An .input.in file is NOT a runnable card. Materialize using explicit values:

```sh
python3 cards/build_r4_catalog.py materialize \\
  cards/R4/measurement/OO_MC_split_allCones.input.in \\
  --values /path/to/real-bindings.json --output /new/path/OO_MC_split.input
```

Example bindings schema (replace values with real checked files):
`{{"MIXING_LIST":"/absolute/path/to/independent_pool.txt"}}`.
Corrected templates also need `PAIR_PAYLOAD`. The native consumer checks current
system/selection tags; independent closure and full jet/track domain coverage are
still required. Initial pair cards train only 40-160 GeV; those maps must not be
applied silently to the high-pT migration region. Extend/qualify the map before
materializing a final-spectrum campaign, without deleting response guards.

`TRACKING_RELATIVE_VARIATION` must come from an applicable Run-3 prescription.
The native systematic histograms alone do not supply a full residual-pair or
finite-training uncertainty. Separate varied pair payloads are needed when the
prescription requires them. See systematics/pending.json for JES/JER, priors,
smearing and other studies that cannot be enabled by card changes alone.

The two unfolding JSON templates use the existing run3_study.py schema and
require the actual projected train/validation products and original split evidence.
Integrated responses must be projected to the explicit -1..100 storage interval;
a split file is not silently substituted for an integrated projection. The
response worker owns synchronized event folds; duplicated card names do not
establish independence. Retain the existing per-jet covariance convention and
separately evaluate donor/training/shared-reference statistics.

R8 is deferred and excluded from this catalog. Its separate configurations
are not promoted to a calibrated production prescription. Do not run them or apply AK4 corrections to them.
'''


def defer_r8():
    tracked = subprocess.check_output(['git', 'ls-files', '*.input'], cwd=REPO, text=True).splitlines()
    moved = {}
    for old in tracked:
        if old.startswith(('cards/R8/', 'cards/archive/')) or not (REPO/old).exists():
            continue
        c = values(old)
        if float(c.get('JetRadius', '0')) != .8:
            continue
        new = 'cards/R8/pending/' + old
        dest = REPO/new
        if dest.exists():
            raise FileExistsError(dest)
        dest.parent.mkdir(parents=True, exist_ok=True)
        data = (REPO/old).read_bytes()
        shutil.move(str(REPO/old), str(dest))
        moved[old] = {'path': new, 'sha256': sha(data)}
    record = REPO/'cards/R8/moved.json'
    if moved:
        prior = json.loads(record.read_text()) if record.exists() else {}
        prior.update(moved)
        record.parent.mkdir(parents=True, exist_ok=True)
        record.write_text(json.dumps(prior, indent=2) + '\n')
        # Static references in active source/documentation follow their moved cards.
        paths = subprocess.check_output(['git', 'ls-files'], cwd=REPO, text=True).splitlines()
        for name in paths:
            p = REPO/name
            if not p.is_file() or name.startswith(('cards/', '.git/')) or p.suffix not in ('.md', '.py', '.sh', '.cxx', '.C', '.h'):
                continue
            text = p.read_text()
            original = text
            for old, info in sorted(moved.items(), key=lambda kv: -len(kv[0])):
                text = re.sub(r'(?<![A-Za-z0-9_/])' + re.escape(old) + r'(?![A-Za-z0-9_.])', info['path'], text)
            # Dynamic card names in retained tests still resolve archived R8 fixtures.
            if name == 'validation/run3Workflow/check_cards.py':
                text = text.replace('path = root / (prefix + code + ".input")',
                                    'path = root / (prefix + code + ".input")\n            if code == "R08":\n                path = root / "cards/R8/pending" / (prefix + code + ".input")')
                text = text.replace('"JetVetoPtConvention"}', '"JetVetoPtConvention", "Run3JetIDTable", "Run3TrackSelection", "Run3TrackPairCorrectionFile"}')
            if text != original:
                p.write_text(text)
    (REPO/'cards/R8').mkdir(parents=True, exist_ok=True)
    (REPO/'cards/R8/README.md').write_text(
        '# R8: deferred, not active production\n\n'
        'All relocated R8 cards retain their original bytes; see moved.json.\n'
        'The pending/ hierarchy mirrors original locations, including synthetic workflow fixtures.\n'
        'They are historical configurations, NOT validated all-raw cards. Native AK8 input,\n'
        'an explicit correction-free execution path and later calibration remain separate work.\n'
        'Do not run these in the R4 campaign or substitute AK4 calibration.\n')
    return moved


def materialize(args):
    source = Path(args.template).resolve(strict=True)
    text = source.read_text()
    bindings = json.loads(Path(args.values).read_text())
    needed = set(TOKEN.findall(text))
    if not needed or not needed <= set(bindings):
        raise ValueError('Missing bindings or input is not a template: ' + str(sorted(needed-set(bindings))))
    for key in needed:
        value = str(bindings[key])
        if not value or any(ch.isspace() for ch in value) or '#' in value or '@' in value:
            raise ValueError('Unsafe/empty card substitution for ' + key)
        if key in ('MIXING_LIST', 'PAIR_PAYLOAD'):
            path = Path(value).expanduser()
            path = path if path.is_absolute() else REPO/path
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(path)
            if key == 'MIXING_LIST' and not any(x.split('#',1)[0].strip() for x in path.read_text().splitlines()):
                raise ValueError('Mixing list contains no inputs')
            value = str(path.resolve())
        if key == 'TRACKING_RELATIVE_VARIATION' and not 0 < float(value) < 1:
            raise ValueError('Tracking variation must be a positive fractional value below one')
        text = text.replace('@' + key + '@', value)
    if TOKEN.search(text):
        raise ValueError('Unresolved placeholder')
    if source.name.endswith('.json.in'):
        json.loads(text)
    else:
        validate(parse(text))
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('x') as stream:
        stream.write(text)
    print(json.dumps({'output': str(out), 'sha256': sha(text.encode()),
                      'status': 'MATERIALIZED_NOT_PHYSICS_QUALIFIED'}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--defer-r8', action='store_true')
    sub = parser.add_subparsers(dest='command')
    m = sub.add_parser('materialize')
    m.add_argument('template'); m.add_argument('--values', required=True); m.add_argument('--output', required=True)
    args = parser.parse_args()
    if args.command == 'materialize':
        materialize(args); return
    if args.write == args.check:
        parser.error('Choose exactly one of --write or --check')
    if args.defer_r8:
        if not args.write:
            parser.error('--defer-r8 requires --write')
        print('R8 relocated:', json.dumps(defer_r8(), sort_keys=True))
    expected = catalog()
    for name, text in expected.items():
        p = REPO/name
        if args.write:
            p.parent.mkdir(parents=True, exist_ok=True); p.write_text(text)
        elif not p.is_file() or p.read_text() != text:
            raise ValueError('Missing/stale catalog file: ' + name)
    extra = set(str(p.relative_to(REPO)) for p in (REPO/'cards/R4').rglob('*') if p.is_file()) - set(expected)
    if extra:
        raise ValueError('Unexpected files in generated catalog: ' + str(sorted(extra)))
    inventory = json.loads(expected['cards/R4/catalog.json'])
    print(json.dumps({'status': 'CONFIGURATION_CHECK_ONLY', 'cards': len(inventory['cards']),
                      'configured': sum(not c['required_tokens'] for c in inventory['cards']),
                      'templates': sum(bool(c['required_tokens']) for c in inventory['cards']),
                      'production_submitted': False, 'physics_closure': 'NOT_RUN'}))


if __name__ == '__main__':
    main()
