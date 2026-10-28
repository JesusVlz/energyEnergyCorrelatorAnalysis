#!/usr/bin/env python3
"""Check dedicated inclusive card coverage; this does not approve physics inputs."""
from pathlib import Path


def parse(path):
    result = {}
    for line in path.read_text().splitlines():
        fields = line.split("#", 1)[0].split()
        if not fields:
            continue
        if fields[0] in result:
            raise ValueError(f"Duplicate key in {path}: {fields[0]}")
        result[fields[0]] = [float(x) for x in fields[1:]]
    return result


def main():
    root = Path(__file__).resolve().parents[2]
    samples = [("cardEECOOInclusive", 9), ("cardEECOOMCInclusive", 10),
               ("cardEECpp5v36TeV", 11), ("cardEECpp5v36TeVMC", 12)]
    for prefix, data_type in samples:
        for code, radius, jet_type in (("R04", .4, 1), ("R08", .8, 3)):
            path = root / (prefix + code + ".input")
            card = parse(path)
            expected = dict(DataType=[data_type], JetType=[jet_type], JetRadius=[radius],
                            JetAxis=[1], JetEtaCut=[1.6], TrackEtaCut=[2.4],
                            MinTrackPtCut=[.5], WeightExponent=[1, 2],
                            InclusiveMode=[1], ChannelRadiusCode=[round(radius*10)],
                            TriggerSelection=[8 if data_type in (9, 11) else 0])
            if data_type in (9, 10):
                expected.update(CentralityBinEdges=[0, 100], InclusiveCentralityTarget=[0, 100])
            for key, value in expected.items():
                if card.get(key) != value:
                    raise ValueError(f"{path}: {key} is {card.get(key)}, expected {value}")
            if not all(x in card["JetPtBinEdgesEEC"] for x in (40, 80, 120, 160)):
                raise ValueError(f"Missing primary jet boundaries: {path}")
            if card["MinJetPtUnfoldingReco"][0] >= 40 or card["MinJetPtUnfoldingTruth"][0] >= 40:
                raise ValueError(f"Missing lower migration guard: {path}")
            print("PASS card metadata:", path.name)
    print("R08 analyzer execution remains BLOCKED on applicable AK8 calibration; pp AK8 input is also required.")


if __name__ == "__main__":
    main()
