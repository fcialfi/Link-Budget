"""Result rows and PDF report for the fixed-elevation (static) link budget.

The row builders turn the dict returned by
:func:`calculations.calculate_fixed_elevation_link_budget` into labelled
rows, and are shared by the GUI results tables and the PDF report so both
always show the same values. The PDF itself is built with ReportLab
(imported lazily, so the rest of the application works without it).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Sequence, Tuple


@dataclass
class ResultRow:
    """One row of the link budget table (the same value in both columns
    unless the clear-sky and rain-faded conditions differ)."""

    group: str
    label: str
    clear: float
    rain: float
    # Pass/fail threshold (value >= threshold passes); None when unchecked.
    threshold: Optional[float] = None
    emphasis: bool = False

    def status(self, value: float) -> Optional[str]:
        if self.threshold is None or value != value:  # NaN is never checked
            return None
        return "pass" if value >= self.threshold else "fail"


@dataclass
class PfdRow:
    label: str
    value: float
    threshold: Optional[float] = None
    digits: int = 2

    def status(self) -> Optional[str]:
        if self.threshold is None or self.value != self.value:
            return None
        return "pass" if self.value >= self.threshold else "fail"


def link_budget_rows(budget: Dict[str, Any], required_margin_db: Optional[float] = None) -> List[ResultRow]:
    """Rows of the two-condition (clear sky / rain faded) link budget table."""

    clear = budget["clear"]
    rain = budget["rain_faded"]
    rows: List[ResultRow] = []

    def both(group: str, label: str, value: float, **kwargs) -> None:
        rows.append(ResultRow(group, label, value, value, **kwargs))

    def per_condition(group: str, label: str, key: str, **kwargs) -> None:
        rows.append(ResultRow(group, label, clear[key], rain[key], **kwargs))

    tx, path, rx, demod, margins = "Transmitter", "Propagation", "Receiver", "Demodulation", "Margins"
    if "eirp_dbw" in budget:
        both(tx, "EIRP [dBW]", budget["eirp_dbw"], emphasis=True)
    both(path, "Free Space Loss [dB]", budget["path_loss_db"])
    per_condition(path, "Gas Attenuation [dB]", "gas_attenuation_db")
    per_condition(path, "Cloud Attenuation [dB]", "cloud_attenuation_db")
    per_condition(path, "Rain Attenuation [dB]", "rain_attenuation_db")
    per_condition(path, "Scintillation [dB]", "scintillation_db")
    per_condition(path, "Total Atmospheric Attenuation [dB]", "atmospheric_attenuation_db")
    if budget.get("ionospheric_loss_db"):
        both(path, "Ionospheric Loss [dB]", budget["ionospheric_loss_db"])
    if budget.get("polarisation_loss_db"):
        both(path, "Polarisation Mismatch Loss [dB]", budget["polarisation_loss_db"])
    if budget.get("multipath_loss_db"):
        both(path, "Multipath Losses [dB]", budget["multipath_loss_db"])
    if budget.get("other_attenuation_db"):
        both(path, "Other Attenuations [dB]", budget["other_attenuation_db"])
    per_condition(path, "Total Propagation Loss [dB]", "total_propagation_loss_db", emphasis=True)
    model = budget.get("gt_model")
    if budget.get("gt_diameter_scaling_db"):
        both(rx, "G/T of the Reference Antenna [dB/K]", budget["gt_reference_dbk"])
        both(rx, "Antenna Diameter Scaling [dB]", budget["gt_diameter_scaling_db"])
    if model:
        both(rx, f"G/T at {model['ref_elevation_deg']:g}° El, clear sky [dB/K]", budget["gt_dbk"])
        if model.get("radome_loss_db"):
            both(rx, "Radome Loss [dB]", model["radome_loss_db"])
        per_condition(rx, "Sky Noise Temperature [K]", "sky_noise_k")
        per_condition(rx, "System Noise Temperature [K]", "system_noise_k")
        per_condition(rx, "Rx G/T [dB/K]", "gt_dbk", emphasis=True)
    elif "gt_dbk" in clear:
        per_condition(rx, "Rx G/T [dB/K]", "gt_dbk")
    elif "gt_dbk" in budget:
        both(rx, "Rx G/T [dB/K]", budget["gt_dbk"])
    if "rx_depointing_deg" in budget:
        both(rx, "Rx Antenna Pointing Loss [dB]", budget["rx_pointing_loss_db"])
    per_condition(rx, "Rx Power [dBW]", "rx_power_dbw")
    if "received_cno_dbhz" in clear:
        per_condition(rx, "Received S/No [dBHz]", "received_cno_dbhz", emphasis=True)
    if budget.get("modulation_degradation_db"):
        both(demod, "Modulation Degradation [dB]", budget["modulation_degradation_db"])
    if "receiver_degradation_db" in budget:
        both(demod, "Receiver Degradation [dB]", budget["receiver_degradation_db"])
    per_condition(demod, "C/No after Degradations [dBHz]", "cno_dbhz")
    if "bitrate_dbhz" in budget:
        both(demod, "Transfer Frame Bit Rate [dBHz]", budget["bitrate_dbhz"])
    per_condition(demod, "Received Eb/No [dB]", "ebno_db", emphasis=True)
    if "required_ebno_db" in budget:
        both(demod, "Required Eb/No [dB]", budget["required_ebno_db"])
    if "margin_db" in clear:
        per_condition(margins, "Data Recovery Margin [dB]", "margin_db", threshold=required_margin_db, emphasis=True)
    clear_stats = clear.get("margin_statistics")
    rain_stats = rain.get("margin_statistics")
    if clear_stats and rain_stats:
        for key, label, threshold, emphasis in (
            ("favourable_db", "Margin – All Favourable [dB]", None, False),
            ("adverse_db", "Margin – All Adverse [dB]", None, False),
            ("mean_db", "Margin Mean [dB]", None, False),
            ("sigma_db", "Margin σ [dB]", None, False),
            ("mean_minus_3sigma_db", "MEAN – 3σ [dB]", 0.0, True),
            ("worst_case_rss_db", "MARGIN – Worst Case RSS [dB]", 0.0, True),
        ):
            rows.append(ResultRow(margins, label, clear_stats[key], rain_stats[key], threshold, emphasis))
    return rows


@dataclass
class ToleranceRow:
    """One row of the link budget with tolerances (ECSS sheet layout)."""

    label: str
    nominal: float
    favourable: float
    adverse: float
    mean: float
    variance: float
    distribution: str
    emphasis: bool = False


# Which tolerance contributors (GUI keys) affect which link budget quantity, and
# the sign of a margin improvement on that quantity (+1: the quantity grows
# when the margin improves, e.g. EIRP; -1: it shrinks, e.g. a loss).
_TOLERANCE_PATH = (
    # (label, nominal key or callable, contributors, sign, emphasis)
    ("Transmitter Power [dBW]", "tx_power_dbw", ("tx_power",), +1, False),
    ("Tx Antenna Effective Gain [dBi]", "effective_gain_dbi", ("tx_gain",), +1, False),
    ("EIRP [dBW]", "eirp_dbw", ("tx_power", "tx_gain"), +1, True),
    ("Total Propagation Loss [dB]", "total_propagation_loss_db", ("propagation",), -1, False),
    ("Rx G/T [dB/K]", "gt_dbk", ("rx_gt",), +1, False),
    ("Received S/No [dBHz]", "received_cno_dbhz", ("tx_power", "tx_gain", "propagation", "rx_gt"), +1, True),
    ("Receiver Degradation [dB]", "receiver_degradation_db", ("rx_degradation",), -1, False),
    (
        "Received Eb/No [dB]", "ebno_db",
        ("tx_power", "tx_gain", "propagation", "rx_gt", "rx_degradation"), +1, True,
    ),
    (
        "Data Recovery Margin [dB]", "margin_db",
        ("tx_power", "tx_gain", "propagation", "rx_gt", "rx_degradation"), +1, True,
    ),
)


def tolerance_rows(budget: Dict[str, Any], condition: str = "rain_faded") -> List[ToleranceRow]:
    """Link budget rows with favourable / adverse / mean / variance columns.

    Each toleranced parameter and every quantity derived from it, down to
    the Data Recovery Margin, of ``condition`` ("clear" or "rain_faded").
    Tolerances are effects on the margin, so they add (with the sign of the
    quantity) along the chain; means and variances of independent
    contributors add. Rows of contributors without a tolerance, and rows
    whose nominal value is not available, are left out.
    """

    cond = budget[condition]
    stats = cond.get("margin_statistics")
    if not stats or not stats["contributions"]:
        return []
    by_key = {c.get("key") or c["name"]: c for c in stats["contributions"]}
    rows: List[ToleranceRow] = []
    for label, nominal_key, contributors, sign, emphasis in _TOLERANCE_PATH:
        involved = [by_key[k] for k in contributors if k in by_key]
        if not involved:
            continue
        nominal = cond.get(nominal_key, budget.get(nominal_key))
        if nominal is None:
            continue
        fav = sum(c["favourable_db"] for c in involved)
        adv = sum(c["adverse_db"] for c in involved)
        mean = sum(c["mean_db"] for c in involved)
        variance = sum(c["variance_db2"] for c in involved)
        distribution = involved[0]["distribution"] if len(involved) == 1 else "combined"
        rows.append(
            ToleranceRow(
                label,
                float(nominal),
                float(nominal + sign * fav),
                float(nominal - sign * adv),
                float(nominal + sign * mean),
                float(variance),
                distribution,
                emphasis,
            )
        )
    return rows


TOLERANCE_TABLE_EXPLANATION = (
    "Each toleranced parameter and the quantities derived from it, down to the margin: nominal "
    "value, value with every tolerance at its favourable / adverse limit, mean and variance. "
    "PDF: TRI triangular peaked at nominal, UNI uniform, GAU Gaussian with the limits at ±3σ; "
    "'combined' where several contributors add (means and variances add). The statistical margins "
    "(MEAN − 3σ, Worst Case RSS) are in the Link Budget table."
)


def pfd_rows(budget: Dict[str, Any], extra_info: Dict[str, Any]) -> List[PfdRow]:
    """Rows of the power flux density check."""

    rows: List[PfdRow] = []
    if "pfd_dbm_m2" in budget:
        rows.append(PfdRow("Power Flux Density [dBm/m²]", budget["pfd_dbm_m2"]))
    rows.append(PfdRow("Power Flux Density [dBW/m²]", budget["pfd_dbw_m2"]))
    if extra_info.get("occupied_bandwidth_khz") is not None:
        rows.append(PfdRow("Occupied Bandwidth [kHz]", extra_info["occupied_bandwidth_khz"], digits=1))
    if "pfd_dbw_m2_per_4khz" in budget:
        rows.append(PfdRow("Power Flux Density in 4 kHz [dBW/m²/4kHz]", budget["pfd_dbw_m2_per_4khz"]))
    if "pfd_limit_high_dbw_m2_4khz" in budget:
        rows.append(PfdRow("PFD Limit 0–5° El [dBW/m²/4kHz]", budget["pfd_limit_low_dbw_m2_4khz"]))
        rows.append(PfdRow("PFD Limit 25–90° El [dBW/m²/4kHz]", budget["pfd_limit_high_dbw_m2_4khz"]))
        rows.append(PfdRow("PFD Limit at this Elevation [dBW/m²/4kHz]", budget["pfd_limit_dbw_m2_4khz"]))
    elif "pfd_limit_dbw_m2_4khz" in budget:
        rows.append(PfdRow("PFD Limit [dBW/m²/4kHz]", budget["pfd_limit_dbw_m2_4khz"]))
    if "pfd_margin_db" in budget:
        rows.append(PfdRow("PFD Margin [dB]", budget["pfd_margin_db"], threshold=0.0))
    return rows


def pfd_notes(budget: Dict[str, Any], extra_info: Dict[str, Any]) -> List[str]:
    """Why part of the PFD check is missing, if it is."""

    if extra_info.get("bandwidth_missing"):
        return [
            "PFD in 4 kHz not computed: set Roll-off and Spectral Efficiency in the "
            "parameters panel to get the occupied bandwidth."
        ]
    if "pfd_margin_db" not in budget:
        return ["Set a PFD Limit (Data Rate & PFD Limit section) to get the PFD margin."]
    return []


# Other Attenuations at or above this value, on top of the computed
# atmospheric losses, look like a total atmospheric allowance counted twice.
OTHER_ATTENUATION_WARNING_DB = 1.0


def analysis_notes(extra_info: Dict[str, Any], budget: Optional[Dict[str, Any]] = None) -> List[str]:
    """Warnings about how far the results can be trusted."""

    notes: List[str] = []
    elevation = extra_info.get("elevation_deg")
    if elevation is not None and elevation < 5.0:
        notes.append(
            f"Elevation {elevation:g}° is below 5°: the ITU-R P.676 gaseous attenuation model used "
            "here is only recommended for elevations of 5–90°, so the atmospheric attenuation is "
            "approximate."
        )
    if extra_info.get("tolerances_ignored"):
        notes.append(
            "Tolerances were entered but not used: set a Required Eb/No to get the statistical "
            "margins (MEAN – 3σ, Worst Case RSS)."
        )
    budget = budget or {}
    if budget.get("scintillation_included") and not budget.get("tropospheric_scintillation_valid", True):
        notes.append(
            f"Tropospheric scintillation (ITU-R P.618) is included, but the method is given for "
            f"elevations of at least 5° and about 4–20 GHz: at {budget.get('frequency_mhz', 0):g} MHz "
            "and this elevation its value is an extrapolation. Consider excluding it; at L/S band "
            "ionospheric scintillation dominates."
        )
    scint = budget.get("ionospheric_scintillation")
    if scint and not scint.get("valid", True):
        notes.append(
            "S4 ≥ 0.4: the f^-1.5 frequency scaling of ITU-R P.531 does not hold for strong "
            "scintillation, so the ionospheric fade is only indicative."
        )
    other = budget.get("other_attenuation_db") or 0.0
    if other >= OTHER_ATTENUATION_WARNING_DB:
        notes.append(
            f"Other Attenuations of {other:g} dB are added on top of the computed atmospheric, "
            "ionospheric and polarisation losses: make sure they are not a total atmospheric "
            "allowance that counts the same effects twice."
        )
    return notes


def pfd_applicable(extra_info: Dict[str, Any]) -> bool:
    """The PFD check applies to space-to-Earth links only (ITU RR Art. 21 limits)."""

    return extra_info.get("pfd_applicable", True)


PFD_EXPLANATION = (
    "The power flux density (PFD) is the power per unit area reaching the receiving site: "
    "PFD = EIRP − 10·log10(4π·d²), with d the slant range. No atmospheric loss is subtracted "
    "(conservative, as required for regulatory checks). To compare it with the ITU limits it is "
    "referred to a 4 kHz bandwidth, assuming the power is spread uniformly over the occupied "
    "bandwidth B: PFD(4 kHz) = PFD − 10·log10(B / 4 kHz). PFD Margin = Limit − PFD(4 kHz) must be "
    "≥ 0 dB. The limit depends on the frequency band and on the elevation angle (ITU Radio "
    "Regulations Art. 21, ECSS-E-ST-50-05C)."
)


def info_rows(budget: Dict[str, Any], extra_info: Dict[str, Any]) -> List[Tuple[str, str]]:
    """Condition-independent figures: data channel, Tx chain and Rx chain details."""

    rows: List[Tuple[str, str]] = []
    if extra_info.get("info_bitrate_kbps") is not None:
        rows.append(("Information (User Data) Bit Rate [kbit/s]", f"{extra_info['info_bitrate_kbps']:.2f}"))
    if budget.get("bitrate_formatted_bps"):
        rows.append(
            ("Transfer Frame Bit Rate (Eb/No ref.) [kbit/s]", f"{budget['bitrate_formatted_bps'] / 1000.0:.2f}")
        )
    if extra_info.get("bitrate_kbps") is not None:
        rows.append(("Coded (Channel) Bit Rate [kbit/s]", f"{extra_info['bitrate_kbps']:.2f}"))
    framing = extra_info.get("ccsds_framing")
    if framing:
        rows.append(
            (
                "CCSDS Transfer Frame [bytes]",
                f"{framing['frame_bytes']} (data field {framing['data_field_bytes']}, user data "
                f"{framing['user_data_bytes']:.1f}, RS I={framing['interleaving_depth']})",
            )
        )
        rows.append(("Formatting / Coding Factor", f"{framing['formatting']:.5f} / {framing['coding']:.5f}"))
    if extra_info.get("code_rate") is not None:
        rows.append(("Overall Rate (information / coded)", f"{extra_info['code_rate']:.3f}"))
    if extra_info.get("occupied_bandwidth_khz") is not None:
        rows.append(("Occupied Bandwidth [kHz]", f"{extra_info['occupied_bandwidth_khz']:.1f}"))
    if "frequency_mhz" in budget:
        rows.append(("Transmitter Frequency [MHz]", f"{budget['frequency_mhz']:.3f}"))
    if "tx_power_dbw" in budget:
        rows.append(("Transmitter Power [dBW]", f"{budget['tx_power_dbw']:.2f}"))
    if "vswr_loss_db" in budget:
        rows.append(("VSWR Losses [dB]", f"{budget['vswr_loss_db']:.2f}"))
    if "effective_gain_dbi" in budget:
        rows.append(("Effective Gain [dBi]", f"{budget['effective_gain_dbi']:.2f}"))
    if "eirp_source" in budget:
        rows.append(
            (
                "EIRP Source",
                "Tx chain (power + effective gain)" if budget["eirp_source"] == "tx_chain" else "Parameters panel",
            )
        )
    if extra_info.get("pointing_loss_db") is not None:
        rows.append(("Tx Antenna Pointing Loss [dB]", f"{extra_info['pointing_loss_db']:.2f}"))
    if extra_info.get("ant_axial_ratio_db") is not None:
        rows.append(("Tx Antenna Axial Ratio [dB]", f"{extra_info['ant_axial_ratio_db']:.2f}"))
    if extra_info.get("rx_axial_ratio_db") is not None:
        rows.append(("Rx Antenna Axial Ratio [dB]", f"{extra_info['rx_axial_ratio_db']:.2f}"))
    if "polarisation_loss_range_db" in budget:
        pol_range = budget["polarisation_loss_range_db"]
        rows.append(
            (
                "Polarisation Loss from Axial Ratios (co-polar gains) [dB]",
                f"{pol_range['co_polar']:.2f}",
            )
        )
    if "rx_depointing_deg" in budget:
        rows.append(("Rx Antenna Depointing Angle [deg]", f"{budget['rx_depointing_deg']:.2f}"))
        rows.append(("Rx Antenna 3dB Beamwidth [deg]", f"{budget['rx_beamwidth_3db_deg']:.2f}"))
    if extra_info.get("station"):
        rows.append(("Station Location", extra_info["station"]))
    if "geomagnetic_latitude_deg" in budget:
        rows.append(
            (
                "Geomagnetic Latitude / Ionospheric Zone",
                f"{budget['geomagnetic_latitude_deg']:.1f}° ({budget['ionospheric_zone']['zone']}, approx.)",
            )
        )
    scint = budget.get("ionospheric_scintillation")
    if scint:
        rows.append(
            (
                "Ionospheric S4 at Link Frequency",
                f"{scint['s4']:.3f} (from {scint['s4_ref']:g} at {scint['ref_freq_mhz']:g} MHz)",
            )
        )
        rows.append(("Ionospheric Peak-to-Peak Fluctuation [dB]", f"{scint['peak_to_peak_db']:.2f}"))
    if extra_info.get("link_availability_pct") is not None:
        rows.append(("Weather Availability [%]", f"{extra_info['link_availability_pct']:g}"))
    return rows


@dataclass
class SummaryCheck:
    """One key result of the summary, checked against its requirement."""

    label: str
    clear: float
    rain: float
    requirement: str
    # "pass", "fail", or None when there is no requirement to check.
    status: Optional[str]


def summary_checks(budget: Dict[str, Any], extra_info: Dict[str, Any]) -> List[SummaryCheck]:
    """Key results (Eb/No, margins, PFD margin) with their pass/fail status.

    A check passes only when it passes in both conditions (clear sky and
    rain faded). The PFD margin does not depend on the weather, so it has
    the same value in both columns.
    """

    required_margin = extra_info.get("required_margin_db")
    by_label = {row.label: row for row in link_budget_rows(budget, required_margin)}
    checks: List[SummaryCheck] = []

    def add(label: str, requirement: str, threshold: Optional[float] = None) -> None:
        row = by_label.get(label)
        if row is None:
            return
        if threshold is not None:
            row = ResultRow(row.group, row.label, row.clear, row.rain, threshold)
        statuses = [row.status(row.clear), row.status(row.rain)]
        overall = "fail" if "fail" in statuses else ("pass" if "pass" in statuses else None)
        checks.append(SummaryCheck(label, row.clear, row.rain, requirement, overall))

    if "required_ebno_db" in budget:
        add("Received Eb/No [dB]", f"≥ {budget['required_ebno_db']:g} dB", budget["required_ebno_db"])
    else:
        add("Received Eb/No [dB]", "–")
    add("Data Recovery Margin [dB]", f"≥ {required_margin:g} dB" if required_margin is not None else "–")
    add("MEAN – 3σ [dB]", "≥ 0 dB")
    add("MARGIN – Worst Case RSS [dB]", "≥ 0 dB")
    if pfd_applicable(extra_info):
        pfd_margin = next((r for r in pfd_rows(budget, extra_info) if r.threshold is not None), None)
        if pfd_margin is not None:
            checks.append(
                SummaryCheck(pfd_margin.label, pfd_margin.value, pfd_margin.value, "≥ 0 dB", pfd_margin.status())
            )
    return checks


def summary_verdict(checks: Sequence[SummaryCheck]) -> Tuple[Optional[str], str]:
    """Overall status ("pass", "fail" or None) and a one-line verdict."""

    failed = [check.label.split(" [")[0] for check in checks if check.status == "fail"]
    if failed:
        return "fail", "Requirement not met: " + ", ".join(failed) + "."
    if any(check.status is not None for check in checks):
        return "pass", "All checked requirements are met in both conditions."
    return None, "No requirement set: enter a Required Eb/No to check the margins."


def rain_percentage(link_availability_pct: float) -> float:
    """Percentage of time p used for the rain-faded condition."""

    return max(0.001, min(50.0, 100.0 - link_availability_pct))


# ---------------------------------------------------------------------------
# PDF report
# ---------------------------------------------------------------------------

@dataclass
class ReportContext:
    """Everything the PDF report shows besides the computed budget."""

    title: str
    direction: str  # e.g. "Downlink (S/C → GS)"
    station: str
    elevation_deg: float
    link_availability_pct: float
    # [(group title, [(label, value text), ...]), ...] -- inputs as entered.
    input_groups: List[Tuple[str, List[Tuple[str, str]]]] = field(default_factory=list)
    generated_at: Optional[datetime] = None
    # Company logo drawn in the page header (light artwork for the navy
    # band); left out when None or unreadable.
    logo_path: Optional[str] = None
    # Project name, shown in the page header and above the title.
    project: str = ""
    # Version of the tool that made the report, shown in the page footer.
    app_version: str = ""


# Colours matching the application's MATLAB-like theme.
_HEADER = "#1f2a36"    # page header band (dark, for the light logo artwork)
_TITLE = "#005a96"     # titles and section headings (MATLAB blue, dark)
_ACCENT = "#0072bd"    # MATLAB blue
_TEXT = "#1a1a1a"
_MUTED = "#5f5f5f"
_BORDER = "#c8c8c8"
_HEAD_BG = "#e6e6e6"   # table column headings, as in MATLAB tables
_ROW_ALT = "#f5f9fd"
_GROUP_BG = "#e3f0fa"
_PASS = "#2e7d32"
_FAIL = "#d64545"

# Characters missing from the standard PDF fonts, drawn with the Symbol font.
_SYMBOL_CHARS = set("σπλΓΔ≥≤−±→")


def _rich(text: str) -> str:
    """Escape ``text`` for a ReportLab Paragraph, using the Symbol font
    for Greek letters and math signs the standard fonts lack."""

    out = []
    for ch in str(text):
        if ch in _SYMBOL_CHARS:
            out.append(f'<font name="Symbol">{ch}</font>')
        elif ch == "&":
            out.append("&amp;")
        elif ch == "<":
            out.append("&lt;")
        elif ch == ">":
            out.append("&gt;")
        elif ch in "✓✗":
            continue
        else:
            out.append(ch)
    return "".join(out)


def _fmt(value: float, digits: int = 2) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "n/a"
    return f"{value:.{digits}f}"


def _load_logo(path: Optional[str]):
    """``(image, width / height)`` for ReportLab, or None when the logo can't be read."""

    if not path:
        return None
    try:
        from reportlab.lib.utils import ImageReader

        image = ImageReader(path)
        width, height = image.getSize()
        return image, width / height
    except Exception:
        return None


def numbered_canvas(header_text: str, right_text: str, stamp: str, app_version: str, logo, margin_x: float):
    """ReportLab canvas class drawing the page header band, footer and "Page x of y"."""

    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.units import mm
    from reportlab.pdfgen import canvas as rl_canvas

    page_w, page_h = A4

    def color(hex_code: str):
        return colors.HexColor(hex_code)

    class NumberedCanvas(rl_canvas.Canvas):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self._pages: list = []

        def showPage(self):
            self._pages.append(dict(self.__dict__))
            self._startPage()

        def save(self):
            total = len(self._pages)
            for state in self._pages:
                self.__dict__.update(state)
                self._decorate(total)
                super().showPage()
            super().save()

        def _decorate(self, total: int):
            self.saveState()
            self.setFillColor(color(_HEADER))
            self.rect(0, page_h - 11 * mm, page_w, 11 * mm, stroke=0, fill=1)
            # MATLAB-blue rule under the band.
            self.setFillColor(color(_ACCENT))
            self.rect(0, page_h - 11.6 * mm, page_w, 0.6 * mm, stroke=0, fill=1)
            text_x = margin_x
            if logo is not None:
                image, aspect = logo
                logo_h = 6.5 * mm
                self.drawImage(
                    image, margin_x, page_h - 8.75 * mm, logo_h * aspect, logo_h, mask="auto",
                )
                text_x = margin_x + logo_h * aspect + 4 * mm
            self.setFillColor(colors.white)
            self.setFont("Helvetica-Bold", 9)
            self.drawString(text_x, page_h - 7 * mm, header_text)
            self.setFont("Helvetica", 8.5)
            self.drawRightString(page_w - margin_x, page_h - 7 * mm, right_text)
            self.setStrokeColor(color(_BORDER))
            self.setLineWidth(0.5)
            self.line(margin_x, 13 * mm, page_w - margin_x, 13 * mm)
            self.setFillColor(color(_MUTED))
            self.setFont("Helvetica", 8)
            tool = f"  ·  Link Budget Tool {app_version}" if app_version else ""
            self.drawString(margin_x, 9 * mm, f"Generated {stamp}{tool}")
            self.drawRightString(page_w - margin_x, 9 * mm, f"Page {self._pageNumber} of {total}")
            self.restoreState()

    return NumberedCanvas


def build_static_report_pdf(
    file_path: str,
    budget: Dict[str, Any],
    extra_info: Dict[str, Any],
    context: ReportContext,
) -> None:
    """Write a compact, paginated A4 PDF report of a fixed-elevation link budget."""

    from reportlab.lib import colors
    from reportlab.lib.enums import TA_RIGHT
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.units import mm
    from reportlab.platypus import (
        CondPageBreak,
        KeepTogether,
        Paragraph,
        SimpleDocTemplate,
        Spacer,
        Table,
        TableStyle,
    )

    generated_at = context.generated_at or datetime.now(timezone.utc)
    page_w, page_h = A4
    margin_x = 18 * mm
    content_w = page_w - 2 * margin_x

    def color(hex_code: str):
        return colors.HexColor(hex_code)

    base = ParagraphStyle("base", fontName="Helvetica", fontSize=9, leading=12, textColor=color(_TEXT))
    styles = {
        "title": ParagraphStyle("title", parent=base, fontName="Helvetica-Bold", fontSize=18, leading=22,
                                textColor=color(_TITLE)),
        "subtitle": ParagraphStyle("subtitle", parent=base, fontSize=10, leading=14, textColor=color(_MUTED)),
        "project": ParagraphStyle("project", parent=base, fontName="Helvetica-Bold", fontSize=9, leading=12,
                                  textColor=color(_ACCENT), spaceAfter=2),
        "h1": ParagraphStyle("h1", parent=base, fontName="Helvetica-Bold", fontSize=12, leading=15,
                             textColor=color(_TITLE), spaceBefore=10, spaceAfter=5),
        "body": base,
        "small": ParagraphStyle("small", parent=base, fontSize=8, leading=10.5, textColor=color(_MUTED)),
        "cell": ParagraphStyle("cell", parent=base, fontSize=8.5, leading=10.5),
        "cell_bold": ParagraphStyle("cell_bold", parent=base, fontName="Helvetica-Bold", fontSize=8.5,
                                    leading=10.5),
        "cell_right": ParagraphStyle("cell_right", parent=base, fontSize=8.5, leading=10.5, alignment=TA_RIGHT),
        "head": ParagraphStyle("head", parent=base, fontName="Helvetica-Bold", fontSize=8.5, leading=10.5,
                               textColor=color(_TEXT)),
        "head_right": ParagraphStyle("head_right", parent=base, fontName="Helvetica-Bold", fontSize=8.5,
                                     leading=10.5, textColor=color(_TEXT), alignment=TA_RIGHT),
        "group": ParagraphStyle("group", parent=base, fontName="Helvetica-Bold", fontSize=8.5, leading=10.5,
                                textColor=color(_TITLE)),
        "verdict": ParagraphStyle("verdict", parent=base, fontName="Helvetica-Bold", fontSize=10, leading=13),
    }

    def heading(text: str) -> list:
        # Start a new page instead of leaving a heading alone at the bottom.
        return [CondPageBreak(35 * mm), Paragraph(_rich(text), styles["h1"])]

    def P(text: str, style: str = "cell") -> Paragraph:
        return Paragraph(_rich(text), styles[style])

    def badge(status: Optional[str]) -> str:
        if status == "pass":
            return f'<font color="{_PASS}"><b>PASS</b></font>'
        if status == "fail":
            return f'<font color="{_FAIL}"><b>FAIL</b></font>'
        return f'<font color="{_MUTED}">–</font>'

    def base_table_style(header_rows: int = 1) -> List[tuple]:
        style = [
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("TOPPADDING", (0, 0), (-1, -1), 3),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
            ("LEFTPADDING", (0, 0), (-1, -1), 5),
            ("RIGHTPADDING", (0, 0), (-1, -1), 5),
            ("LINEBELOW", (0, 0), (-1, -1), 0.4, color(_BORDER)),
            ("BOX", (0, 0), (-1, -1), 0.6, color(_BORDER)),
        ]
        if header_rows:
            style.append(("BACKGROUND", (0, 0), (-1, header_rows - 1), color(_HEAD_BG)))
            style.append(("LINEBELOW", (0, header_rows - 1), (-1, header_rows - 1), 0.8, color(_MUTED)))
        return style

    story: list = []

    # ---- Title block -------------------------------------------------------
    if context.project:
        story.append(Paragraph(f"PROJECT: {_rich(context.project.upper())}", styles["project"]))
    story.append(Paragraph(_rich(context.title), styles["title"]))
    subtitle = (
        f"{context.direction}  ·  {context.station}  ·  "
        f"{budget.get('frequency_mhz', float('nan')):.1f} MHz  ·  elevation {context.elevation_deg:g}°  ·  "
        f"slant range {budget['slant_range_km']:.1f} km"
    )
    story.append(Paragraph(_rich(subtitle), styles["subtitle"]))
    story.append(Spacer(1, 3 * mm))

    # ---- Summary -----------------------------------------------------------
    rain_p = rain_percentage(context.link_availability_pct)

    summary_data = [[P("Key result", "head"), P("Clear Sky", "head_right"),
                     P(f"Rain Faded (p = {rain_p:g}%)", "head_right"), P("Requirement", "head"), P("Status", "head")]]
    checks = summary_checks(budget, extra_info)
    for check in checks:
        summary_data.append(
            [P(check.label, "cell_bold"), P(_fmt(check.clear), "cell_right"), P(_fmt(check.rain), "cell_right"),
             P(check.requirement), Paragraph(badge(check.status), styles["cell"])]
        )
    summary_table = Table(
        summary_data,
        colWidths=[content_w * f for f in (0.34, 0.14, 0.22, 0.17, 0.13)],
        repeatRows=1,
    )
    summary_style = base_table_style()
    for i in range(1, len(summary_data)):
        if i % 2 == 0:
            summary_style.append(("BACKGROUND", (0, i), (-1, i), color(_ROW_ALT)))
    summary_table.setStyle(TableStyle(summary_style))

    overall, verdict_text = summary_verdict(checks)
    if overall == "fail":
        verdict = f'<font color="{_FAIL}">Requirement not met:</font> ' + _rich(verdict_text.split(": ", 1)[1])
    elif overall == "pass":
        verdict = f'<font color="{_PASS}">All checked requirements are met</font> in both conditions.'
    else:
        verdict = _rich(verdict_text)

    story.append(Paragraph("Summary", styles["h1"]))
    story.append(Paragraph(verdict, styles["verdict"]))
    story.append(Spacer(1, 2 * mm))
    if len(summary_data) > 1:
        story.append(summary_table)
    for note in analysis_notes(extra_info, budget):
        story.append(Spacer(1, 1.5 * mm))
        story.append(Paragraph(f'<font color="{_FAIL}"><b>Note:</b></font> {_rich(note)}', styles["small"]))

    # ---- Inputs ------------------------------------------------------------
    input_groups = [(title, items) for title, items in context.input_groups if items]
    # Derived figures that only echo an input (or the availability) are left out.
    entered = {label for _title, items in input_groups for label, _value in items}
    info = [
        (label, value)
        for label, value in info_rows(budget, extra_info)
        if label not in entered and label != "Weather Availability [%]"
    ]
    if info:
        input_groups.append(("Derived Parameters", info))
    if input_groups:
        story.extend(heading("Input and Derived Parameters"))
        data: list = []
        style = base_table_style(header_rows=0)
        for title, items in input_groups:
            style.append(("SPAN", (0, len(data)), (-1, len(data))))
            style.append(("BACKGROUND", (0, len(data)), (-1, len(data)), color(_GROUP_BG)))
            data.append([P(title, "group"), "", "", ""])
            for k in range(0, len(items), 2):
                pair = list(items[k:k + 2]) + [("", "")] * (2 - len(items[k:k + 2]))
                data.append([P(pair[0][0]), P(pair[0][1], "cell_bold"), P(pair[1][0]), P(pair[1][1], "cell_bold")])
        table = Table(data, colWidths=[content_w * f for f in (0.33, 0.17, 0.33, 0.17)])
        table.setStyle(TableStyle(style))
        story.append(table)

    # ---- Link budget -------------------------------------------------------
    story.extend(heading("Link Budget"))
    data = [[P("Parameter", "head"), P("Clear Sky", "head_right"), P(f"Rain Faded (p = {rain_p:g}%)", "head_right")]]
    style = base_table_style()
    current_group = None
    for row in link_budget_rows(budget, extra_info.get("required_margin_db")):
        if row.group != current_group:
            current_group = row.group
            style.append(("SPAN", (0, len(data)), (-1, len(data))))
            style.append(("BACKGROUND", (0, len(data)), (-1, len(data)), color(_GROUP_BG)))
            data.append([P(row.group, "group"), "", ""])
        label_style = "cell_bold" if row.emphasis else "cell"
        cells = [P(row.label, label_style)]
        for value in (row.clear, row.rain):
            status = row.status(value)
            text = _fmt(value)
            if status:
                text = f'{text}  {badge(status)}'
            weight = "<b>%s</b>" % text if row.emphasis else text
            cells.append(Paragraph(weight, styles["cell_right"]))
        data.append(cells)
    table = Table(data, colWidths=[content_w * f for f in (0.52, 0.22, 0.26)], repeatRows=1)
    table.setStyle(TableStyle(style))
    story.append(table)
    story.append(Spacer(1, 1.5 * mm))
    story.append(Paragraph(
        _rich("Clear Sky: gas absorption only. Rain Faded: gas, cloud and rain attenuation (and "
              "scintillation when enabled) exceeded for p% of the time, ITU-R P.618."),
        styles["small"],
    ))

    # ---- PFD ---------------------------------------------------------------
    pfd = pfd_rows(budget, extra_info) if pfd_applicable(extra_info) else []
    pfd_block = [Paragraph("Power Flux Density (regulatory check)", styles["h1"]),
                 Paragraph(_rich(PFD_EXPLANATION), styles["body"]), Spacer(1, 2 * mm)]
    data = [[P("Parameter", "head"), P("Value", "head_right"), P("Status", "head")]]
    for row in pfd:
        data.append([P(row.label, "cell_bold" if row.threshold is not None else "cell"),
                     P(_fmt(row.value, row.digits), "cell_right"),
                     Paragraph(badge(row.status()), styles["cell"])])
    table = Table(data, colWidths=[content_w * f for f in (0.6, 0.22, 0.18)])
    table.setStyle(TableStyle(base_table_style()))
    pfd_block.append(table)
    for note in pfd_notes(budget, extra_info):
        pfd_block.append(Spacer(1, 1.5 * mm))
        pfd_block.append(Paragraph(f'<font color="{_FAIL}">{_rich(note)}</font>', styles["small"]))
    if pfd:
        story.append(KeepTogether(pfd_block))

    # ---- Tolerances --------------------------------------------------------
    tol_rows = tolerance_rows(budget, "rain_faded")
    if tol_rows:
        data = [[P("Parameter", "head")] + [P(h, "head_right") for h in (
            "Nominal", "Favourable", "Adverse", "Mean", "Variance [dB²]")] + [P("PDF", "head")]]
        style = base_table_style()
        for i, row in enumerate(tol_rows, start=1):
            weight = "cell_bold" if row.emphasis else "cell"
            data.append(
                [P(row.label, weight)]
                + [P(f"{v:.2f}", "cell_right") for v in (row.nominal, row.favourable, row.adverse, row.mean)]
                + [P(f"{row.variance:.3f}", "cell_right"), P(row.distribution)]
            )
            if row.emphasis:
                style.append(("BACKGROUND", (0, i), (-1, i), color(_ROW_ALT)))
        table = Table(data, colWidths=[content_w * f for f in (0.31, 0.11, 0.12, 0.11, 0.1, 0.13, 0.12)])
        table.setStyle(TableStyle(style))
        story.append(KeepTogether([
            Paragraph(f"Link Budget with Tolerances (Rain Faded, p = {rain_p:g}%)", styles["h1"]),
            Paragraph(_rich(TOLERANCE_TABLE_EXPLANATION), styles["body"]),
            Spacer(1, 2 * mm),
            table,
        ]))

    # ---- Method ------------------------------------------------------------
    method = [
        "Slant range from the elevation angle and a circular orbit at the given altitude over a "
        "spherical Earth (no TLE involved): a single, typically worst-case, geometry.",
        "Atmospheric attenuation from the ITU-R P.618 slant-path model at the ground station "
        "location (gas P.676, cloud P.840, rain P.618/P.837"
        + (", tropospheric scintillation P.618)." if budget.get("scintillation_included") else
           "; tropospheric scintillation excluded)."),
        "Received S/No = Rx Power + G/T + 228.6 dB (Boltzmann constant); Eb/No referred to the "
        "transfer frame bit rate (coded bit rate without channel coding = information bit rate × "
        "formatting factor), after modulation and receiver degradations.",
        "Statistical margins assume independent contributors (ECSS/CCSDS link budget convention).",
    ]
    if budget.get("ionospheric_scintillation"):
        method.insert(
            2,
            "Ionospheric scintillation fade from the S4 index (ITU-R P.531): S4 scaled to the link "
            "frequency as f^-1.5, peak-to-peak power fluctuation 27.5·S4^1.26 dB, fade taken as half "
            "of it; method and reference S4 (Kiruna, 99.9% of the time) from ESA technical note "
            "AWS-TN-ESA-SYS-0052. The geomagnetic zone of the station is from a dipole model. For the "
            "L-band downlink to Svalbard the note finds 1 to more than 3 dB of ionospheric fading and "
            "1.8 to more than 3 dB of total (tropospheric + ionospheric) attenuation, and recommends "
            "conservative values until site measurements are available.",
        )
    story.append(KeepTogether(
        [Paragraph("Method and Assumptions", styles["h1"])]
        + [Paragraph("• " + _rich(item), styles["body"]) for item in method]
    ))

    # ---- Page decoration (header, footer, "page x of y") -------------------
    header_text = "Link Budget Tool  ·  Static Link Budget Report" + (
        f"  ·  {context.project}" if context.project else ""
    )
    stamp = generated_at.strftime("%Y-%m-%d %H:%M UTC")

    NumberedCanvas = numbered_canvas(
        header_text, context.direction.replace("→", "->"), stamp, context.app_version,
        _load_logo(context.logo_path), margin_x,
    )

    doc = SimpleDocTemplate(
        file_path,
        pagesize=A4,
        leftMargin=margin_x,
        rightMargin=margin_x,
        topMargin=18 * mm,
        bottomMargin=18 * mm,
        title=f"{context.project} – {context.title}" if context.project else context.title,
        author="Exostaff",
        creator=f"Link Budget Tool {context.app_version}".strip(),
        subject=f"{context.direction} static link budget",
    )
    doc.build(story, canvasmaker=NumberedCanvas)
