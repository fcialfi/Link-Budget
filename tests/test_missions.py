import json
import os

import pytest

import missions

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_repository_catalogue():
    cat = missions.load_catalogue(os.path.join(REPO, "satellites.json"))
    assert cat.default == "AWS"
    aws = cat.find("AWS")
    assert aws.norad_id == 60543 and not aws.custom_orbit
    assert cat.parameters_path(aws).endswith(os.path.join("missions", "AWS", "parameters.json"))
    sterna = cat.find("EPS-Sterna")
    assert sterna.custom_orbit and sterna.norad_id is None
    assert cat.parameters_path(cat.find("Metop-B")) is None
    assert "custom orbit" in sterna.label and "NORAD 60543" in aws.label


def test_mission_parameters():
    with open(os.path.join(REPO, "missions", "AWS", "parameters.json"), encoding="utf-8") as f:
        aws = json.load(f)
    with open(os.path.join(REPO, "missions", "EPS-Sterna", "parameters.json"), encoding="utf-8") as f:
        sterna = json.load(f)
    assert aws["norad_id"] == 60543
    assert aws["fixed_elevation_downlink"]["sat_altitude_km"] == 600.0
    assert sterna["orbit"]["altitude_km"] == 595.6
    assert sterna["orbit"]["ltan"] == ["15:30", "19:30", "23:30"]
    assert sterna["fixed_elevation_downlink"]["sat_altitude_km"] == 628.0


def test_catalogue_errors(tmp_path):
    path = tmp_path / "satellites.json"
    path.write_text(json.dumps({"satellites": [{"name": "X"}]}))
    with pytest.raises(ValueError, match="norad_id"):
        missions.load_catalogue(str(path))
    path.write_text(json.dumps({"satellites": [{"name": "X", "norad_id": 1}, {"name": "X", "norad_id": 2}]}))
    with pytest.raises(ValueError, match="twice"):
        missions.load_catalogue(str(path))
    path.write_text(json.dumps({"default": "Y", "satellites": [{"name": "A/B", "norad_id": 1}]}))
    cat = missions.load_catalogue(str(path))
    assert cat.default is None
    assert cat.mission_dir(cat.find("A/B")).endswith(os.path.join("missions", "A_B"))
