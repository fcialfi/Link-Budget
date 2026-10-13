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
    per_condition(path, "Total Propagation Loss [dB]", "total_propagation_loss_db", emphasis=True)
    if "gt_dbk" in budget:
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
        both(demod, "Bit Rate incl. Formatting [dBHz]", budget["bitrate_dbhz"])
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
    if "pfd_limit_dbw_m2_4khz" in budget:
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
        rows.append(("Information Bit Rate (excl. coding) [kbit/s]", f"{extra_info['info_bitrate_kbps']:.2f}"))
    if budget.get("bitrate_formatted_bps"):
        rows.append(
            ("Bit Rate incl. Formatting (excl. coding) [kbit/s]", f"{budget['bitrate_formatted_bps'] / 1000.0:.2f}")
        )
    if extra_info.get("bitrate_kbps") is not None:
        rows.append(("Bit Rate incl. Coding [kbit/s]", f"{extra_info['bitrate_kbps']:.2f}"))
    if extra_info.get("code_rate") is not None:
        rows.append(("Code Rate", f"{extra_info['code_rate']:.3f}"))
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
                "Polarisation Loss from Axial Ratios (best/worst) [dB]",
                f"{pol_range['best']:.2f} / {pol_range['worst']:.2f}",
            )
        )
    if "rx_depointing_deg" in budget:
        rows.append(("Rx Antenna Depointing Angle [deg]", f"{budget['rx_depointing_deg']:.2f}"))
        rows.append(("Rx Antenna 3dB Beamwidth [deg]", f"{budget['rx_beamwidth_3db_deg']:.2f}"))
    if extra_info.get("station"):
        rows.append(("Station Location", extra_info["station"]))
    if extra_info.get("link_availability_pct") is not None:
        rows.append(("Weather Availability [%]", f"{extra_info['link_availability_pct']:g}"))
    return rows


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


# Colours matching the application theme.
_NAVY = "#0f2440"
_ACCENT = "#2f6fed"
_TEXT = "#1f2733"
_MUTED = "#5b6472"
_BORDER = "#d7dde6"
_ROW_ALT = "#f3f6fb"
_GROUP_BG = "#e8f0fe"
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
    from reportlab.pdfgen import canvas as rl_canvas
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
                                textColor=color(_NAVY)),
        "subtitle": ParagraphStyle("subtitle", parent=base, fontSize=10, leading=14, textColor=color(_MUTED)),
        "h1": ParagraphStyle("h1", parent=base, fontName="Helvetica-Bold", fontSize=12, leading=15,
                             textColor=color(_NAVY), spaceBefore=10, spaceAfter=5),
        "body": base,
        "small": ParagraphStyle("small", parent=base, fontSize=8, leading=10.5, textColor=color(_MUTED)),
        "cell": ParagraphStyle("cell", parent=base, fontSize=8.5, leading=10.5),
        "cell_bold": ParagraphStyle("cell_bold", parent=base, fontName="Helvetica-Bold", fontSize=8.5,
                                    leading=10.5),
        "cell_right": ParagraphStyle("cell_right", parent=base, fontSize=8.5, leading=10.5, alignment=TA_RIGHT),
        "head": ParagraphStyle("head", parent=base, fontName="Helvetica-Bold", fontSize=8.5, leading=10.5,
                               textColor=colors.white),
        "head_right": ParagraphStyle("head_right", parent=base, fontName="Helvetica-Bold", fontSize=8.5,
                                     leading=10.5, textColor=colors.white, alignment=TA_RIGHT),
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
            style.append(("BACKGROUND", (0, 0), (-1, header_rows - 1), color(_NAVY)))
        return style

    story: list = []

    # ---- Title block -------------------------------------------------------
    story.append(Paragraph(_rich(context.title), styles["title"]))
    subtitle = (
        f"{context.direction}  ·  {context.station}  ·  "
        f"{budget.get('frequency_mhz', float('nan')):.1f} MHz  ·  elevation {context.elevation_deg:g}°  ·  "
        f"slant range {budget['slant_range_km']:.1f} km"
    )
    story.append(Paragraph(_rich(subtitle), styles["subtitle"]))
    story.append(Spacer(1, 3 * mm))

    # ---- Summary -----------------------------------------------------------
    rows = link_budget_rows(budget, extra_info.get("required_margin_db"))
    by_label = {row.label: row for row in rows}
    rain_p = rain_percentage(context.link_availability_pct)
    required_margin = extra_info.get("required_margin_db")

    summary_data = [[P("Key result", "head"), P("Clear Sky", "head_right"),
                     P(f"Rain Faded (p = {rain_p:g}%)", "head_right"), P("Requirement", "head"), P("Status", "head")]]
    checks: List[Tuple[str, Optional[str]]] = []

    def summary_row(label: str, requirement: str, threshold: Optional[float] = None) -> None:
        row = by_label.get(label)
        if row is None:
            return
        if threshold is not None:
            row = ResultRow(row.group, row.label, row.clear, row.rain, threshold)
        statuses = [row.status(row.clear), row.status(row.rain)]
        overall = "fail" if "fail" in statuses else ("pass" if "pass" in statuses else None)
        checks.append((label, overall))
        summary_data.append(
            [P(label, "cell_bold"), P(_fmt(row.clear), "cell_right"), P(_fmt(row.rain), "cell_right"),
             P(requirement), Paragraph(badge(overall), styles["cell"])]
        )

    if "required_ebno_db" in budget:
        summary_row("Received Eb/No [dB]", f"≥ {budget['required_ebno_db']:g} dB", budget["required_ebno_db"])
    else:
        summary_row("Received Eb/No [dB]", "–")
    summary_row("Data Recovery Margin [dB]", f"≥ {required_margin:g} dB" if required_margin is not None else "–")
    summary_row("MEAN – 3σ [dB]", "≥ 0 dB")
    summary_row("MARGIN – Worst Case RSS [dB]", "≥ 0 dB")
    pfd_margin_row = next((r for r in pfd_rows(budget, extra_info) if r.threshold is not None), None)
    if pfd_margin_row is not None:
        status = pfd_margin_row.status()
        checks.append((pfd_margin_row.label, status))
        summary_data.append(
            [P(pfd_margin_row.label, "cell_bold"), P(_fmt(pfd_margin_row.value), "cell_right"),
             P(_fmt(pfd_margin_row.value), "cell_right"), P("≥ 0 dB"), Paragraph(badge(status), styles["cell"])]
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

    failed = [label for label, status in checks if status == "fail"]
    checked = [label for label, status in checks if status is not None]
    if failed:
        verdict = (f'<font color="{_FAIL}">Requirement not met:</font> '
                   + ", ".join(_rich(label.split(" [")[0]) for label in failed) + ".")
    elif checked:
        verdict = f'<font color="{_PASS}">All checked requirements are met</font> in both conditions.'
    else:
        verdict = "No requirement set: enter a Required Eb/No to check the margins."

    story.append(Paragraph("Summary", styles["h1"]))
    story.append(Paragraph(verdict, styles["verdict"]))
    story.append(Spacer(1, 2 * mm))
    if len(summary_data) > 1:
        story.append(summary_table)

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
            data.append([P(title, "cell_bold"), "", "", ""])
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
    for row in rows:
        if row.group != current_group:
            current_group = row.group
            style.append(("SPAN", (0, len(data)), (-1, len(data))))
            style.append(("BACKGROUND", (0, len(data)), (-1, len(data)), color(_GROUP_BG)))
            data.append([P(row.group, "cell_bold"), "", ""])
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
    pfd = pfd_rows(budget, extra_info)
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
    story.append(KeepTogether(pfd_block))

    # ---- Tolerances --------------------------------------------------------
    stats = budget["clear"].get("margin_statistics")
    if stats and stats["contributions"]:
        contributions = stats["contributions"]
        data = [[P(h, "head") for h in ("Contributor", "Favourable [dB]", "Adverse [dB]", "PDF",
                                         "Mean [dB]", "Variance [dB²]")]]
        for c in contributions:
            data.append([P(c["name"]), P(f"+{c['favourable_db']:.2f}", "cell_right"),
                         P(f"−{c['adverse_db']:.2f}", "cell_right"), P(c["distribution"]),
                         P(f"{c['mean_db']:.2f}", "cell_right"), P(f"{c['variance_db2']:.3f}", "cell_right")])
        data.append([P("Total", "cell_bold"),
                     P(f"+{sum(c['favourable_db'] for c in contributions):.2f}", "cell_right"),
                     P(f"−{sum(c['adverse_db'] for c in contributions):.2f}", "cell_right"), "",
                     P(f"{sum(c['mean_db'] for c in contributions):.2f}", "cell_right"),
                     P(f"{stats['variance_db2']:.3f}", "cell_right")])
        table = Table(data, colWidths=[content_w * f for f in (0.32, 0.15, 0.14, 0.09, 0.14, 0.16)])
        style = base_table_style()
        style.append(("LINEABOVE", (0, -1), (-1, -1), 0.8, color(_NAVY)))
        table.setStyle(TableStyle(style))
        story.append(KeepTogether([
            Paragraph("Tolerances", styles["h1"]),
            Paragraph(_rich(
                "Favourable/adverse deviations of each contributor, as effect on the margin. TRI: triangular "
                "peaked at nominal; UNI: uniform; GAU: Gaussian with the limits at ±3σ. Contributors are "
                "independent: means and variances add, σ is the square root of the total variance, and the "
                "Worst Case RSS margin is the nominal margin minus the root-sum-square of the adverse values."
            ), styles["body"]),
            Spacer(1, 2 * mm),
            table,
        ]))

    # ---- Method ------------------------------------------------------------
    method = [
        "Slant range from the elevation angle and a circular orbit at the given altitude over a "
        "spherical Earth (no TLE involved): a single, typically worst-case, geometry.",
        "Atmospheric attenuation from the ITU-R P.618 slant-path model at the ground station "
        "location (gas P.676, cloud P.840, rain P.618/P.837, scintillation P.618).",
        "Received S/No = Rx Power + G/T + 228.6 dB (Boltzmann constant); Eb/No referred to the bit "
        "rate including formatting and excluding coding, after modulation and receiver degradations.",
        "Statistical margins assume independent contributors (ECSS/CCSDS link budget convention).",
    ]
    story.append(KeepTogether(
        [Paragraph("Method and Assumptions", styles["h1"])]
        + [Paragraph("• " + _rich(item), styles["body"]) for item in method]
    ))

    # ---- Page decoration (header, footer, "page x of y") -------------------
    header_text = "Link Budget Tool  ·  Static Link Budget Report"
    stamp = generated_at.strftime("%Y-%m-%d %H:%M UTC")

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
            self.setFillColor(color(_NAVY))
            self.rect(0, page_h - 11 * mm, page_w, 11 * mm, stroke=0, fill=1)
            self.setFillColor(colors.white)
            self.setFont("Helvetica-Bold", 9)
            self.drawString(margin_x, page_h - 7 * mm, header_text)
            self.setFont("Helvetica", 8.5)
            self.drawRightString(page_w - margin_x, page_h - 7 * mm, context.direction.replace("→", "->"))
            self.setStrokeColor(color(_BORDER))
            self.setLineWidth(0.5)
            self.line(margin_x, 13 * mm, page_w - margin_x, 13 * mm)
            self.setFillColor(color(_MUTED))
            self.setFont("Helvetica", 8)
            self.drawString(margin_x, 9 * mm, f"Generated {stamp}")
            self.drawRightString(page_w - margin_x, 9 * mm, f"Page {self._pageNumber} of {total}")
            self.restoreState()

    doc = SimpleDocTemplate(
        file_path,
        pagesize=A4,
        leftMargin=margin_x,
        rightMargin=margin_x,
        topMargin=18 * mm,
        bottomMargin=18 * mm,
        title=context.title,
        author="Link Budget Tool",
        subject=f"{context.direction} static link budget",
    )
    doc.build(story, canvasmaker=NumberedCanvas)
