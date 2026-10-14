"""Satellite catalogue and mission folders.

``satellites.json`` lists the satellites that can be selected in the
Orbit & Passes tab::

    {"default": "AWS",
     "satellites": [
        {"name": "AWS", "norad_id": 60543, "description": "Arctic Weather Satellite"},
        {"name": "EPS-Sterna", "description": "...", "custom_orbit": true}
     ]}

A satellite with ``norad_id`` uses the TLE downloaded from the internet; one
with ``custom_orbit`` (not launched yet) uses a synthetic TLE built from the
``orbit`` section of its parameters (altitude, LTAN, optional inclination).

The link budget parameters of a satellite are in
``missions/<name>/parameters.json`` (same format as the old top-level
parameters.json), with an optional ``antenna_pattern.txt``. The TLE last
downloaded is kept in ``missions/<name>/tle.txt``. A satellite without
parameters gets the contact analysis only.
"""
from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass

CATALOGUE_FILE = "satellites.json"
MISSIONS_DIR = "missions"
PARAMETERS_FILE = "parameters.json"
TLE_FILE = "tle.txt"
PATTERN_FILE = "antenna_pattern.txt"


@dataclass
class Satellite:
    name: str
    norad_id: int | None = None
    description: str = ""
    custom_orbit: bool = False

    @property
    def label(self) -> str:
        """Text of the satellite menu."""
        text = self.name
        if self.description:
            text += f" – {self.description}"
        if self.custom_orbit:
            text += " (custom orbit)"
        elif self.norad_id is not None:
            text += f" (NORAD {self.norad_id})"
        return text


@dataclass
class Catalogue:
    satellites: list[Satellite]
    default: str | None
    base_dir: str

    def find(self, name: str) -> Satellite | None:
        for sat in self.satellites:
            if sat.name == name:
                return sat
        return None

    def mission_dir(self, sat: Satellite) -> str:
        folder = re.sub(r'[<>:"/\\|?*]', "_", sat.name).strip() or "unnamed"
        return os.path.join(self.base_dir, MISSIONS_DIR, folder)

    def mission_file(self, sat: Satellite, filename: str) -> str:
        return os.path.join(self.mission_dir(sat), filename)

    def parameters_path(self, sat: Satellite) -> str | None:
        """The mission's parameters file, or None when it has none."""
        path = self.mission_file(sat, PARAMETERS_FILE)
        return path if os.path.isfile(path) else None


def load_catalogue(path: str) -> Catalogue:
    """Read ``satellites.json``; raise ValueError when it is not valid."""
    with open(path, "r", encoding="utf-8") as f:
        payload = json.load(f)
    items = payload.get("satellites") if isinstance(payload, dict) else None
    if not isinstance(items, list) or not items:
        raise ValueError(f"{path}: no 'satellites' list")
    satellites = []
    names = set()
    for item in items:
        name = str(item.get("name", "")).strip()
        if not name:
            raise ValueError(f"{path}: a satellite without name")
        if name in names:
            raise ValueError(f"{path}: satellite {name} listed twice")
        names.add(name)
        norad = item.get("norad_id")
        custom = bool(item.get("custom_orbit", False))
        if norad is None and not custom:
            raise ValueError(f"{path}: {name} needs a norad_id or \"custom_orbit\": true")
        satellites.append(
            Satellite(
                name=name,
                norad_id=int(norad) if norad is not None else None,
                description=str(item.get("description", "")),
                custom_orbit=custom,
            )
        )
    default = payload.get("default")
    return Catalogue(satellites, default if default in names else None, os.path.dirname(os.path.abspath(path)))
