"""Supplied event-veto predicate; activation requires approved inputs.

No default map, Jet ID, pT convention or calibration coordinate is inferred.
Evaluate *all* eligible reconstructed jets, including those outside analysis pT
bins, and apply the same approved selection to signal and mixing events.
"""
import math


def reject_event(jets, map_value, passes_id, *, pt_key, eta_key, phi_key):
    for jet in jets:
        pt, eta, phi = (float(jet[key]) for key in (pt_key, eta_key, phi_key))
        cef, nef = float(jet['charged_em']), float(jet['neutral_em'])
        if not all(math.isfinite(x) for x in (pt, eta, phi, cef, nef)):
            raise ValueError('Invalid veto coordinate/fraction; no zero-content fallback')
        if pt <= 15 or not passes_id(jet) or cef+nef >= .9:
            continue
        value = float(map_value(eta, phi))
        if not math.isfinite(value):
            raise ValueError('Invalid map value; no zero-content fallback')
        if value != 0:
            return True
    return False
