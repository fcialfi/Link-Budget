"""PDF report of the statistics over several days (contacts, data volume, latency).

Built with ReportLab like the static link budget report, with the same page
header and footer (:func:`static_report.numbered_canvas`).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import List, Optional, Tuple

import static_report
from static_report import _ACCENT, _BORDER, _FAIL, _HEAD_BG, _MUTED, _PASS, _ROW_ALT, _TEXT, _TITLE, _rich


@dataclass
class StatisticsReportContext:
    """Everything the report shows besides the computed result."""

    title: str  # e.g. "EPS-Sterna – Statistics"
    period: str  # e.g. "EPS-Sterna over Svalbard + Punta Arenas, 2026-10-14 to 2026-10-20, min elevation 5°"
    # [(label, value text), ...]: the inputs used.
    settings: List[Tuple[str, str]] = field(default_factory=list)
    # Summary lines as shown in the window.
    summary_lines: List[str] = field(default_factory=list)
    figure_png: Optional[str] = None
    generated_at: Optional[datetime] = None
    logo_path: Optional[str] = None
    project: str = ""
    app_version: str = ""


def _num(value, digits: int = 1) -> str:
    return "–" if value is None else f"{value:.{digits}f}"


def build_statistics_report_pdf(file_path: str, result, context: StatisticsReportContext,
                                max_pass_rows: int = 400) -> None:
    """Write an A4 PDF: verdict, settings, summary, plots and the tables of ``result``."""

    from reportlab.lib import colors
    from reportlab.lib.enums import TA_RIGHT
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.units import mm
    from reportlab.platypus import (
        CondPageBreak, Image, KeepTogether, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle,
    )

    generated_at = context.generated_at or datetime.now(timezone.utc)
    page_w, _page_h = A4
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
        "cell": ParagraphStyle("cell", parent=base, fontSize=8, leading=10),
        "cell_right": ParagraphStyle("cell_right", parent=base, fontSize=8, leading=10, alignment=TA_RIGHT),
        "head": ParagraphStyle("head", parent=base, fontName="Helvetica-Bold", fontSize=8, leading=10),
        "head_right": ParagraphStyle("head_right", parent=base, fontName="Helvetica-Bold", fontSize=8,
                                     leading=10, alignment=TA_RIGHT),
        "verdict": ParagraphStyle("verdict", parent=base, fontName="Helvetica-Bold", fontSize=10, leading=13),
    }

    def heading(text: str) -> list:
        return [CondPageBreak(35 * mm), Paragraph(_rich(text), styles["h1"])]

    def table(head: List[str], rows: List[List[str]], widths: List[float], left_cols: int = 1) -> Table:
        data = [[Paragraph(_rich(h), styles["head"] if k < left_cols else styles["head_right"])
                 for k, h in enumerate(head)]]
        for row in rows:
            data.append([Paragraph(_rich(v), styles["cell"] if k < left_cols else styles["cell_right"])
                         for k, v in enumerate(row)])
        t = Table(data, colWidths=[content_w * w for w in widths], repeatRows=1)
        style = [
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("TOPPADDING", (0, 0), (-1, -1), 2.5),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
            ("LEFTPADDING", (0, 0), (-1, -1), 4),
            ("RIGHTPADDING", (0, 0), (-1, -1), 4),
            ("LINEBELOW", (0, 0), (-1, -1), 0.4, color(_BORDER)),
            ("BOX", (0, 0), (-1, -1), 0.6, color(_BORDER)),
            ("BACKGROUND", (0, 0), (-1, 0), color(_HEAD_BG)),
            ("LINEBELOW", (0, 0), (-1, 0), 0.8, color(_MUTED)),
        ]
        for i in range(2, len(data), 2):
            style.append(("BACKGROUND", (0, i), (-1, i), color(_ROW_ALT)))
        t.setStyle(TableStyle(style))
        return t

    summary = result.summary
    story: list = []
    if context.project:
        story.append(Paragraph(f"PROJECT: {_rich(context.project.upper())}", styles["project"]))
    story.append(Paragraph(_rich(context.title), styles["title"]))
    story.append(Paragraph(_rich(context.period), styles["subtitle"]))
    story.append(Spacer(1, 3 * mm))

    # ---- Verdict on the latency requirement --------------------------------
    story.append(Paragraph("Summary", styles["h1"]))
    if "latency_requirement_min" in summary:
        ok = summary.get("latency_requirement_met")
        text = (
            f"{summary['latency_within_pct']:.1f} % of the data on ground within "
            f"{summary['latency_requirement_min']:g} min (required: {summary['latency_requirement_pct']:g} %); "
            f"{summary['latency_requirement_pct']:g} % of the data within "
            f"{summary['latency_at_required_pct_min']:.1f} min."
        )
        verdict = (f'<font color="{_PASS}">Latency requirement met:</font> ' if ok
                   else f'<font color="{_FAIL}">Latency requirement not met:</font> ') + _rich(text)
        story.append(Paragraph(verdict, styles["verdict"]))
        story.append(Spacer(1, 2 * mm))
    # The window separates the figures of a line with wide spaces: one item per line here.
    for line in context.summary_lines:
        for item in (part.strip() for part in line.split("    ")):
            if item:
                story.append(Paragraph("• " + _rich(item), styles["body"]))

    # ---- Settings ----------------------------------------------------------
    if context.settings:
        story += heading("Inputs")
        story.append(table(["Input", "Value"], [[k, v] for k, v in context.settings], [0.4, 0.6], left_cols=2))

    # ---- Plots -------------------------------------------------------------
    if context.figure_png:
        from reportlab.lib.utils import ImageReader

        width, height = ImageReader(context.figure_png).getSize()
        story.append(KeepTogether([
            Paragraph("Plots", styles["h1"]),
            Image(context.figure_png, width=content_w, height=content_w * height / width),
        ]))

    # ---- Tables ------------------------------------------------------------
    if getattr(result, "stations", None):
        story += heading("Per Station")
        story.append(table(
            ["Station", "Antennas", "Served [min/day]", "Volume [MB/day]", "Busy [%]", "Contention [min/day]"],
            [[s["station"], str(s.get("antennas", 1)), _num(s["served_min_per_day"]),
              _num(s["downloaded_mb_per_day"], 0), _num(s["busy_pct"]), _num(s["contention_min_per_day"])]
             for s in result.stations],
            [0.22, 0.12, 0.17, 0.17, 0.12, 0.2],
        ))
    if getattr(result, "satellites", None):
        story += heading("Per Satellite")
        rows = []
        for s in result.satellites:
            rows.append([
                s["satellite"], _num(s["served_min_per_day"]), _num(s["conflict_min_per_day"]),
                _num(s["downloaded_mb_per_day"], 0), _num(s.get("max_on_board_mb"), 0),
                f"{_num(s.get('mean_latency_min'))} / {_num(s.get('max_latency_min'))}",
                _num(s.get("latency_within_pct")),
            ])
        story.append(table(
            ["Satellite", "Served [min/day]", "Conflicts [min/day]", "Volume [MB/day]", "Max on board [MB]",
             "Latency mean/max [min]", "Within req. [%]"],
            rows, [0.22, 0.12, 0.13, 0.12, 0.13, 0.16, 0.12],
        ))
    story += heading("Per Day")
    story.append(table(
        ["Date", "Contacts", "Contact [min]", "Usable [min]", "Capacity [MB]", "Downloaded [MB]",
         "Max gap [min]", "On board [MB]"],
        [[str(d.date), str(d.contacts), _num(d.contact_min), _num(d.usable_min), _num(d.capacity_mb, 0),
          _num(d.downloaded_mb, 0), _num(d.max_gap_min), _num(d.on_board_mb, 0)] for d in result.days],
        [0.14, 0.1, 0.12, 0.12, 0.12, 0.14, 0.13, 0.13],
    ))
    passes = result.passes[:max_pass_rows]
    story += heading("Per Pass")
    if len(result.passes) > max_pass_rows:
        story.append(Paragraph(
            _rich(f"First {max_pass_rows} of {len(result.passes)} passes; all of them are in the CSV export."),
            styles["small"],
        ))
    story.append(table(
        ["Start (UTC)", "End", "Stations", "Max El [°]", "Contact [min]", "Usable [min]", "Volume [MB]",
         "Min margin [dB]"],
        [[f"{p.start:%Y-%m-%d %H:%M:%S}", f"{p.end:%H:%M:%S}", p.stations, _num(p.max_elevation_deg),
          _num(p.contact_min), _num(p.usable_min), _num(p.downloaded_mb, 1), _num(p.min_margin_db, 2)]
         for p in passes],
        [0.17, 0.09, 0.27, 0.08, 0.09, 0.09, 0.11, 0.1], left_cols=3,
    ))

    header_text = "Link Budget Tool  ·  Statistics Report" + (f"  ·  {context.project}" if context.project else "")
    canvas_class = static_report.numbered_canvas(
        header_text, "Statistics", generated_at.strftime("%Y-%m-%d %H:%M UTC"), context.app_version,
        static_report._load_logo(context.logo_path), margin_x,
    )
    doc = SimpleDocTemplate(
        file_path, pagesize=A4, leftMargin=margin_x, rightMargin=margin_x, topMargin=18 * mm,
        bottomMargin=18 * mm, title=f"{context.project} – {context.title}" if context.project else context.title,
        author="Exostaff", creator=f"Link Budget Tool {context.app_version}".strip(), subject="Statistics",
    )
    doc.build(story, canvasmaker=canvas_class)
