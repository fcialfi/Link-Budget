"""Cesium/CesiumPy export helpers for the link budget tool."""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import json
import os
import tempfile
from typing import Iterable, Sequence

from skyfield.api import EarthSatellite, load, wgs84

try:  # Optional dependency
    import cesiumpy  # type: ignore
except Exception:  # pragma: no cover - optional dependency
    cesiumpy = None


DEFAULT_CESIUM_JS_URL = (
    "https://cesium.com/downloads/cesiumjs/releases/1.114/Build/Cesium/Cesium.js"
)
DEFAULT_CESIUM_CSS_URL = (
    "https://cesium.com/downloads/cesiumjs/releases/1.114/Build/Cesium/Widgets/widgets.css"
)


@dataclass(frozen=True)
class CesiumExportPaths:
    """Paths to the generated Cesium artifacts."""

    czml_path: str
    html_path: str


def _isoformat_z(dt: datetime) -> str:
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return dt.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")


def _build_czml(
    tle1: str,
    tle2: str,
    gs_name: str,
    gs_lat: float,
    gs_lon: float,
    gs_alt_m: float,
    times: Sequence[datetime],
) -> list[dict[str, object]]:
    if not times:
        raise ValueError("At least one timestamp is required to build CZML.")

    ts = load.timescale()
    sat = EarthSatellite(tle1, tle2, "SAT", ts)

    sky_times = ts.utc(
        [t.year for t in times],
        [t.month for t in times],
        [t.day for t in times],
        [t.hour for t in times],
        [t.minute for t in times],
        [t.second for t in times],
    )

    geocentric = sat.at(sky_times)
    subpoints = wgs84.subpoint(geocentric)
    lats = subpoints.latitude.degrees
    lons = subpoints.longitude.degrees
    alts = subpoints.elevation.m

    epoch = times[0]
    epoch_iso = _isoformat_z(epoch)
    position_samples: list[float] = []
    for idx, t in enumerate(times):
        offset = (t - epoch).total_seconds()
        position_samples.extend(
            [offset, float(lons[idx]), float(lats[idx]), float(alts[idx])]
        )

    start = _isoformat_z(times[0])
    end = _isoformat_z(times[-1])

    return [
        {
            "id": "document",
            "name": "Link Budget Cesium Export",
            "version": "1.0",
            "clock": {
                "interval": f"{start}/{end}",
                "currentTime": start,
                "multiplier": 20,
                "range": "LOOP_STOP",
                "step": "SYSTEM_CLOCK_MULTIPLIER",
            },
        },
        {
            "id": "satellite",
            "name": "Satellite",
            "availability": f"{start}/{end}",
            "position": {
                "epoch": epoch_iso,
                "cartographicDegrees": position_samples,
            },
            "path": {
                "material": {
                    "solidColor": {"color": {"rgba": [0, 255, 255, 200]}}
                },
                "width": 2,
                "leadTime": 0,
                "trailTime": 3600,
            },
            "point": {
                "color": {"rgba": [255, 255, 0, 220]},
                "pixelSize": 8,
            },
        },
        {
            "id": "ground-station",
            "name": gs_name,
            "position": {
                "cartographicDegrees": [float(gs_lon), float(gs_lat), float(gs_alt_m)],
            },
            "label": {
                "text": gs_name,
                "font": "14px sans-serif",
                "fillColor": {"rgba": [255, 255, 255, 255]},
                "style": "FILL_AND_OUTLINE",
                "outlineWidth": 2,
                "verticalOrigin": "BOTTOM",
            },
            "point": {
                "color": {"rgba": [255, 0, 0, 220]},
                "pixelSize": 10,
            },
        },
        {
            "id": "link",
            "name": "Link Line",
            "polyline": {
                "positions": {
                    "references": [
                        "satellite#position",
                        "ground-station#position",
                    ]
                },
                "material": {
                    "solidColor": {"color": {"rgba": [0, 255, 0, 150]}}
                },
                "width": 1.5,
                "followSurface": False,
            },
        },
    ]


def _write_html(
    html_path: str,
    czml_filename: str,
    title: str,
    cesium_js_url: str,
    cesium_css_url: str,
) -> None:
    html = f"""<!DOCTYPE html>
<html lang=\"en\">
  <head>
    <meta charset=\"utf-8\" />
    <title>{title}</title>
    <script src=\"{cesium_js_url}\"></script>
    <link rel=\"stylesheet\" href=\"{cesium_css_url}\" />
    <style>
      html, body, #cesiumContainer {{
        width: 100%;
        height: 100%;
        margin: 0;
        padding: 0;
        overflow: hidden;
      }}
    </style>
  </head>
  <body>
    <div id=\"cesiumContainer\"></div>
    <script>
      const viewer = new Cesium.Viewer('cesiumContainer', {{
        timeline: true,
        animation: true,
        shouldAnimate: true,
      }});
      const dataSource = new Cesium.CzmlDataSource();
      viewer.dataSources.add(dataSource);
      dataSource.load('{czml_filename}').then(function() {{
        viewer.zoomTo(dataSource);
      }});
    </script>
  </body>
</html>
"""
    with open(html_path, "w", encoding="utf-8") as handle:
        handle.write(html)


def export_cesium_bundle(
    tle1: str,
    tle2: str,
    gs_name: str,
    gs_lat: float,
    gs_lon: float,
    gs_alt_m: float,
    times: Sequence[datetime],
    output_base: str,
    title: str = "Link Budget Cesium View",
    cesium_js_url: str = DEFAULT_CESIUM_JS_URL,
    cesium_css_url: str = DEFAULT_CESIUM_CSS_URL,
) -> CesiumExportPaths:
    """Write CZML + HTML bundle for the selected time window."""

    if not output_base:
        raise ValueError("Output base path cannot be empty.")

    base = os.path.splitext(output_base)[0]
    czml_path = f"{base}.czml"
    html_path = f"{base}.html"

    czml = _build_czml(tle1, tle2, gs_name, gs_lat, gs_lon, gs_alt_m, times)
    with open(czml_path, "w", encoding="utf-8") as handle:
        json.dump(czml, handle, indent=2)

    _write_html(
        html_path,
        os.path.basename(czml_path),
        title,
        cesium_js_url,
        cesium_css_url,
    )

    if cesiumpy is not None and hasattr(cesiumpy, "Cesium"):
        try:  # Best-effort CesiumPy integration
            viewer = cesiumpy.Cesium(czml_path)
            viewer.save(html_path)
        except Exception:
            pass

    return CesiumExportPaths(czml_path=czml_path, html_path=html_path)


def preview_cesium_view(
    tle1: str,
    tle2: str,
    gs_name: str,
    gs_lat: float,
    gs_lon: float,
    gs_alt_m: float,
    times: Sequence[datetime],
    title: str = "Link Budget Cesium Preview",
) -> CesiumExportPaths:
    """Create a temporary Cesium bundle and return the paths."""

    if cesiumpy is None or not hasattr(cesiumpy, "Cesium"):
        raise RuntimeError("CesiumPy is not available.")

    tmp_dir = tempfile.mkdtemp(prefix="link_budget_cesium_")
    output_base = os.path.join(tmp_dir, "cesium_preview")
    paths = export_cesium_bundle(
        tle1,
        tle2,
        gs_name,
        gs_lat,
        gs_lon,
        gs_alt_m,
        times,
        output_base,
        title=title,
    )
    return paths


def slice_times(
    times: Iterable[datetime],
    start: datetime,
    end: datetime,
) -> list[datetime]:
    """Return a sorted list of timestamps between ``start`` and ``end``."""

    return sorted([t for t in times if start <= t <= end])
