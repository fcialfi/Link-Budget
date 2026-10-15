import json
import os

import pytest

import mission_schema as ms

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_rows_cover_the_eps_sterna_file():
    with open(os.path.join(HERE, "missions", "EPS-Sterna", "parameters.json"), encoding="utf-8") as f:
        payload = json.load(f)
    rows = ms.rows(payload)
    by_path = {r["path"]: r for r in rows}
    assert by_path[("orbit", "phasing_deg")]["value"] == 60.0
    assert by_path[("downlink_data_kbps", "instruments")]["unit"] == "kbps"
    assert by_path[("fixed_elevation_downlink", "required_ebno_db")]["group"] == "Static LB · Downlink"
    assert all(r["group"] != "Other" or r["path"][0] not in ("orbit", "fixed_elevation_downlink") for r in rows)
    # Every scalar top-level key of the file is shown somewhere.
    shown = {r["path"][0] for r in rows}
    assert all(k in shown for k, v in payload.items() if not isinstance(v, dict))


def test_parse_and_set_values():
    payload = {"orbit": {"ltan": ["15:30"]}, "memory_gb": 4.0}
    ms.set_path(payload, ("orbit", "ltan"), ms.parse_value("15:30, 19:30", "list"))
    ms.set_path(payload, ("downlink_data_kbps", "hktm"), ms.parse_value("3,105", "float"))
    ms.set_path(payload, ("memory_gb",), ms.parse_value("", "float"))
    ms.set_path(payload, ("latency_at_pass_end",), ms.parse_value("yes", "bool"))
    assert payload == {"orbit": {"ltan": ["15:30", "19:30"]}, "downlink_data_kbps": {"hktm": 3.105},
                       "latency_at_pass_end": True}
    with pytest.raises(ValueError):
        ms.parse_value("abc", "float")
    assert ms.format_value(True) == "yes" and ms.format_value(["a", "b"]) == "a, b"


def test_standardize_adds_every_key_as_null():
    payload = {"project_name": "AWS", "frequency_ghz": 1.707, "fixed_elevation_downlink": {"elevation_deg": 3.0},
               "custom_key": 1}
    std = ms.standardize(payload, custom_orbit=False)
    assert std["project_name"] == "AWS" and std["frequency_ghz"] == 1.707
    assert std["downlink_data_kbps"] == {"instruments": None, "hktm": None, "other": None, "realtime": None}
    assert std["memory_gb"] is None and std["latency_requirement"] == {"minutes": None, "percent": None}
    assert "orbit" not in std  # not a custom-orbit mission
    assert std["fixed_elevation_downlink"] == {"elevation_deg": 3.0} and std["custom_key"] == 1
    assert list(std)[0] == "project_name"
    assert ms.standardize(payload, custom_orbit=True)["orbit"]["altitude_km"] is None
    # Idempotent, and the values of an existing section are kept.
    assert ms.standardize(std, custom_orbit=False) == std
