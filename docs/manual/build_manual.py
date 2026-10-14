#!/usr/bin/env python3
"""Build the PDF user manual: docs/Link_Budget_Tool_User_Manual.pdf.

    python docs/manual/take_screenshots.py   # when the GUI layout changed
    python docs/manual/build_manual.py

The text lives in this file, the screenshots in docs/manual/img/ (made by
take_screenshots.py). The version comes from version.py and the
parameters.json key table from the GUI source, so they always match.
"""

import datetime
import os
import re
import sys

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib.utils import ImageReader
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (BaseDocTemplate, CondPageBreak, Frame, Image, KeepTogether, NextPageTemplate,
                                PageBreak, PageTemplate, Paragraph, Spacer, Table, TableStyle)
from reportlab.platypus.tableofcontents import TableOfContents

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
IMG = os.path.join(HERE, "img")
OUT = os.path.join(ROOT, "docs", "Link_Budget_Tool_User_Manual.pdf")
LOGO = os.path.join(ROOT, "assets", "exostaff_logo.png")
sys.path.insert(0, ROOT)
from version import APP_AUTHOR, APP_VERSION  # noqa: E402

TITLE = "Satellite Link Budget Tool"
SUBTITLE = "User Manual"
ISSUE = "Issue 1.0"

# Colours of the application (MATLAB-like theme).
HEADER = colors.HexColor("#1f2a36")
ACCENT = colors.HexColor("#0072bd")
TITLE_BLUE = colors.HexColor("#005a96")
TEXT = colors.HexColor("#1a1a1a")
MUTED = colors.HexColor("#5f5f5f")
HEAD_BG = colors.HexColor("#e6e6e6")
GROUP_BG = colors.HexColor("#e3f0fa")
ROW_ALT = colors.HexColor("#f5f9fd")
BORDER = colors.HexColor("#c8c8c8")


# --------------------------------------------------------------------------
# Fonts and styles
# --------------------------------------------------------------------------
def _register_fonts():
    """DejaVu when available (σ, ≥, →, ° ...), else the PDF base fonts."""
    for folder in ("/usr/share/fonts/truetype/dejavu", "C:/Windows/Fonts", "/Library/Fonts"):
        files = {"": "DejaVuSans.ttf", "-Bold": "DejaVuSans-Bold.ttf", "-Mono": "DejaVuSansMono.ttf"}
        if all(os.path.exists(os.path.join(folder, f)) for f in files.values()):
            for suffix, f in files.items():
                pdfmetrics.registerFont(TTFont("Body" + suffix, os.path.join(folder, f)))
            pdfmetrics.registerFontFamily("Body", normal="Body", bold="Body-Bold", italic="Body",
                                          boldItalic="Body-Bold")
            return "Body", "Body-Bold", "Body-Mono"
    return "Helvetica", "Helvetica-Bold", "Courier"


FONT, BOLD, MONO = _register_fonts()
BASE = ParagraphStyle("base", fontName=FONT, fontSize=9.5, leading=13.5, textColor=TEXT, spaceAfter=5)
S = {
    "h1": ParagraphStyle("h1", parent=BASE, fontName=BOLD, fontSize=17, leading=21, textColor=TITLE_BLUE,
                         spaceBefore=4, spaceAfter=10),
    "h2": ParagraphStyle("h2", parent=BASE, fontName=BOLD, fontSize=12.5, leading=16, textColor=TITLE_BLUE,
                         spaceBefore=10, spaceAfter=5),
    "h3": ParagraphStyle("h3", parent=BASE, fontName=BOLD, fontSize=10.5, leading=14, textColor=TEXT,
                         spaceBefore=6, spaceAfter=3),
    "body": BASE,
    "note": ParagraphStyle("note", parent=BASE, fontSize=9, leading=12.5, textColor=TEXT, backColor=GROUP_BG,
                           borderPadding=(5, 7, 5, 7), spaceBefore=4, spaceAfter=9, leftIndent=7, rightIndent=7),
    "bullet": ParagraphStyle("bullet", parent=BASE, leftIndent=14, bulletIndent=4, spaceAfter=2.5),
    "caption": ParagraphStyle("caption", parent=BASE, fontSize=8.5, leading=11, textColor=MUTED, alignment=1,
                              spaceBefore=3, spaceAfter=10),
    "cell": ParagraphStyle("cell", parent=BASE, fontSize=8.5, leading=11, spaceAfter=0),
    "cell_bold": ParagraphStyle("cell_bold", parent=BASE, fontName=BOLD, fontSize=8.5, leading=11, spaceAfter=0),
    "cell_mono": ParagraphStyle("cell_mono", parent=BASE, fontName=MONO, fontSize=8, leading=11, spaceAfter=0),
    "code": ParagraphStyle("code", parent=BASE, fontName=MONO, fontSize=8.5, leading=11.5, backColor=HEAD_BG,
                           borderPadding=(4, 6, 4, 6), leftIndent=6, rightIndent=6, spaceBefore=3, spaceAfter=8),
    "toc1": ParagraphStyle("toc1", parent=BASE, fontName=BOLD, fontSize=10, leading=13, leftIndent=0,
                           spaceAfter=0, spaceBefore=4),
    "toc2": ParagraphStyle("toc2", parent=BASE, fontSize=9, leading=11.5, leftIndent=16, spaceAfter=0),
}


# --------------------------------------------------------------------------
# Building blocks
# --------------------------------------------------------------------------
_chapter = [0]
_section = [0]


def _esc(text):
    return text.replace("&", "&amp;")


def H1(text):
    text = _esc(text)
    _chapter[0] += 1
    _section[0] = 0
    title = f"{_chapter[0]}  {text}"
    p = Paragraph(title, S["h1"])
    p.toc = (0, title)
    return [PageBreak(), p]


def H2(text):
    text = _esc(text)
    _section[0] += 1
    title = f"{_chapter[0]}.{_section[0]}  {text}"
    p = Paragraph(title, S["h2"])
    p.toc = (1, title)
    return [CondPageBreak(3 * cm), p]


def H3(text):
    return Paragraph(_esc(text), S["h3"])


def P(text):
    return Paragraph(text, S["body"])


def NOTE(text):
    return Paragraph(text, S["note"])


def BUL(items):
    return [Paragraph(item, S["bullet"], bulletText="•") for item in items]


def CODE(text):
    return Paragraph(text.replace("\n", "<br/>"), S["code"])


_FIGURES = [0]


def FIG(name, caption, width_cm=16.5):
    """A screenshot with its caption, numbered in order of appearance."""
    _FIGURES[0] += 1
    caption = f"Figure {_FIGURES[0]} – {caption}"
    path = os.path.join(IMG, f"{name}.png")
    if not os.path.exists(path):
        return [P(f"<i>[missing screenshot {name}.png: run take_screenshots.py]</i>")]
    w, h = ImageReader(path).getSize()
    width = width_cm * cm
    return [KeepTogether([Image(path, width=width, height=width * h / w), Paragraph(caption, S["caption"])])]


def TABLE(rows, widths_cm, header=True, mono_cols=(), bold_cols=(0,)):
    data = []
    for i, row in enumerate(rows):
        cells = []
        for j, cell in enumerate(row):
            if header and i == 0:
                style = "cell_bold"
            elif j in mono_cols:
                style = "cell_mono"
            elif j in bold_cols:
                style = "cell_bold"
            else:
                style = "cell"
            cells.append(Paragraph(str(cell), S[style]))
        data.append(cells)
    table = Table(data, colWidths=[w * cm for w in widths_cm], repeatRows=1 if header else 0)
    style = [
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("BOX", (0, 0), (-1, -1), 0.6, BORDER),
        ("LINEBELOW", (0, 0), (-1, -1), 0.4, BORDER),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
    ]
    if header:
        style += [("BACKGROUND", (0, 0), (-1, 0), HEAD_BG), ("LINEBELOW", (0, 0), (-1, 0), 0.8, MUTED)]
    for i in range(1 if header else 0, len(data)):
        if (i % 2 == 0) == header:
            style.append(("BACKGROUND", (0, i), (-1, i), ROW_ALT))
    table.setStyle(TableStyle(style))
    return [table, Spacer(1, 8)]


def parameter_keys():
    """parameters.json keys of the main panel, read from the GUI source."""
    with open(os.path.join(ROOT, "gui.py"), encoding="utf-8") as f:
        source = f.read()
    block = source[source.index("def _parameter_entries"):]
    block = block[:block.index("}")]
    return re.findall(r'"([a-z_0-9]+)": self\.', block)


# --------------------------------------------------------------------------
# Page decoration
# --------------------------------------------------------------------------
def _cover(canv, doc):
    w, h = A4
    canv.saveState()
    canv.setFillColor(HEADER)
    canv.rect(0, h - 11 * cm, w, 11 * cm, stroke=0, fill=1)
    canv.setFillColor(ACCENT)
    canv.rect(0, h - 11.25 * cm, w, 0.25 * cm, stroke=0, fill=1)
    if os.path.exists(LOGO):
        lw, lh = ImageReader(LOGO).getSize()
        width = 7 * cm
        canv.drawImage(LOGO, 2.2 * cm, h - 4.6 * cm, width=width, height=width * lh / lw, mask="auto")
    canv.setFillColor(colors.white)
    canv.setFont(BOLD, 28)
    canv.drawString(2.2 * cm, h - 7.4 * cm, TITLE)
    canv.setFont(FONT, 16)
    canv.setFillColor(colors.HexColor("#b8c4d0"))
    canv.drawString(2.2 * cm, h - 8.6 * cm, SUBTITLE)
    canv.setFillColor(TEXT)
    canv.setFont(FONT, 11)
    lines = [
        f"Software version {APP_VERSION}",
        ISSUE,
        f"Author: {APP_AUTHOR}",
        datetime.date.today().strftime("%d %B %Y"),
    ]
    for i, line in enumerate(lines):
        canv.drawString(2.2 * cm, h - 13.5 * cm - i * 0.65 * cm, line)
    canv.setFillColor(MUTED)
    canv.setFont(FONT, 9)
    canv.drawString(2.2 * cm, 2.2 * cm, "Exostaff")
    canv.restoreState()


def _page(canv, doc):
    w, h = A4
    canv.saveState()
    canv.setFillColor(HEADER)
    canv.rect(0, h - 1.15 * cm, w, 1.15 * cm, stroke=0, fill=1)
    canv.setFillColor(ACCENT)
    canv.rect(0, h - 1.21 * cm, w, 0.06 * cm, stroke=0, fill=1)
    x = 2 * cm
    if os.path.exists(LOGO):
        lw, lh = ImageReader(LOGO).getSize()
        height = 0.62 * cm
        canv.drawImage(LOGO, x, h - 0.88 * cm, width=height * lw / lh, height=height, mask="auto")
        x += height * lw / lh + 0.4 * cm
    canv.setFillColor(colors.white)
    canv.setFont(BOLD, 9)
    canv.drawString(x, h - 0.72 * cm, f"{TITLE}  ·  {SUBTITLE}")
    canv.setFillColor(MUTED)
    canv.setFont(FONT, 8)
    canv.setStrokeColor(BORDER)
    canv.setLineWidth(0.5)
    canv.line(2 * cm, 1.5 * cm, w - 2 * cm, 1.5 * cm)
    canv.drawString(2 * cm, 1.05 * cm, f"Version {APP_VERSION}  ·  {ISSUE}")
    canv.drawRightString(w - 2 * cm, 1.05 * cm, f"Page {doc.page}")
    canv.restoreState()


class ManualDoc(BaseDocTemplate):
    def afterFlowable(self, flowable):
        toc = getattr(flowable, "toc", None)
        if toc:
            level, text = toc
            key = f"h{id(flowable)}"
            self.canv.bookmarkPage(key)
            self.canv.addOutlineEntry(text, key, level=level, closed=level > 0)
            self.notify("TOCEntry", (level, text, self.page, key))


# --------------------------------------------------------------------------
# Content
# --------------------------------------------------------------------------
def content():
    story = []
    story.append(NextPageTemplate("page"))
    story.append(PageBreak())
    story.append(Paragraph("Contents", S["h1"]))
    toc = TableOfContents()
    toc.levelStyles = [S["toc1"], S["toc2"]]
    # No extra cell padding, so the contents fit on one page.
    toc.tableStyle = TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 0),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
        ("LEFTPADDING", (0, 0), (-1, -1), 0),
        ("RIGHTPADDING", (0, 0), (-1, -1), 0),
    ])
    story.append(toc)

    # 1 -------------------------------------------------------------------
    story += H1("Introduction")
    story.append(P(
        "The Satellite Link Budget Tool computes the radio link budget between a low Earth orbit satellite and "
        "a ground station, for the downlink (satellite to ground) and the uplink (ground to satellite). It "
        "offers two complementary analyses:"))
    story += BUL([
        "<b>Pass analysis</b>: from the satellite TLE, the tool propagates the orbit over a whole day, finds the "
        "passes (contact windows) over the selected ground station and computes, every 30 seconds of each "
        "pass, the elevation, slant range, Doppler shift and the downlink and uplink link budgets.",
        "<b>Static link budget</b>: the classic preliminary, worst-case check at a single fixed elevation "
        "(typically the minimum operational one), independent of any TLE, laid out as an ECSS static link "
        "budget sheet with clear-sky and rain-faded conditions, statistical margins from favourable/adverse "
        "tolerances, the power flux density (PFD) regulatory check and a PDF report.",
    ])
    story.append(P(
        "Both analyses use the same models: ITU-R atmospheric attenuation, ionospheric scintillation from an S4 "
        "index, polarisation mismatch, the satellite antenna pattern and the CCSDS telemetry framing, so that "
        "their results can be compared."))
    story.append(H2("Conventions")[1])
    story += BUL([
        "Names of buttons, fields and tabs are written in <b>bold</b>, as they appear on screen.",
        "Values are in dB, dBW, dBi, dB/K, dBHz unless stated; angles in degrees; bit rates as labelled.",
        "Every input field has an <b>(i)</b> icon next to it: hovering over it shows what the field means and how it is used. "
        "Every row of the result tables shows its definition when hovered.",
    ])

    # 2 -------------------------------------------------------------------
    story += H1("Installation and start")
    story += H2("Program folder")
    story.append(P(
        "The program is delivered as a folder, <b>Link-Budget-Tool</b>, usually as a zip file. Unzip it anywhere "
        "(e.g. C:\\Programs\\Exostaff\\Link-Budget-Tool); no Python installation is needed. The folder contains "
        "<b>Link-Budget-Tool.exe</b>, its libraries (in _internal) and the data files below, which the "
        "program reads at startup and which can be edited or replaced without rebuilding it:"))
    story += TABLE([
        ["File", "Content"],
        ["license.key", "Licence of the installation (see 2.2). Required when the build checks licences."],
        ["parameters.json", "Project name, link budget parameters of the main tabs and the inputs of the two "
                            "static link budget tabs (chapter 9)."],
        ["tle.txt", "Two-line elements of the satellite."],
        ["ground_stations.txt", "Ground station catalogue: name, latitude, longitude, altitude (m)."],
        ["Antenna Pattern Default.txt", "Satellite antenna pattern: angle from nadir (deg), gain (dBi)."],
    ], [5, 11.5], mono_cols=(0,), bold_cols=())
    story.append(P(
        "A missing data file is not an error: the related fields stay empty (or the built-in station list and "
        "antenna pattern are used) and can be filled in, or loaded from another file, in the program."))
    story += H2("Licence")
    story.append(P(
        "The executable can be time-limited. The licence is the file <b>license.key</b>, placed next to "
        "Link-Budget-Tool.exe; it names the customer and the expiry date and is signed, so it cannot be "
        "edited. A valid licence is shown at the top right of the window ('Licensed to … · valid until …', in "
        "orange during its last 30 days). Without a valid licence the program stops at startup with one of "
        "these messages:"))
    story += TABLE([
        ["Message", "What to do"],
        ["No licence found", "Copy the license.key you received into the program folder."],
        ["The licence expired on …", "Ask the supplier for a new licence file."],
        ["The licence file is not valid / damaged", "The file was edited or corrupted: use the original file."],
        ["The system date appears to have been set back", "Correct the computer date and restart."],
    ], [6.5, 10])
    story += H2("Starting the program")
    story.append(P(
        "Double-click <b>Link-Budget-Tool.exe</b>. A window with the Exostaff logo appears immediately and "
        "shows the progress while the program loads (a few seconds). The main window then opens with the data "
        "files already loaded. Close the program with <b>Exit</b> (top right), Ctrl+Q or the window close "
        "button."))

    # 3 -------------------------------------------------------------------
    story += H1("Main window")
    story += FIG("orbit_passes", "Main window, Orbit &amp; Passes tab.")
    story.append(P(
        "The header shows the program version, the licence and the <b>About</b> (version, "
        "author and licence; also F1) and <b>Exit</b> buttons. "
        "Below it, the tabs:"))
    story += TABLE([
        ["Tab", "Content"],
        ["Orbit &amp; Passes", "Project name, TLE, date and ground station, data files, the pass analysis and "
                               "the list of passes (chapter 4)."],
        ["Downlink", "Downlink parameters, and plots and table of the selected pass (chapter 5)."],
        ["Uplink", "Uplink parameters and table of the selected pass (chapter 6)."],
        ["Static LB · Downlink", "Fixed-elevation downlink link budget (chapter 7)."],
        ["Static LB · Uplink", "Fixed-elevation uplink link budget (chapter 7)."],
    ], [4.5, 12])
    story.append(NOTE(
        "<b>Shared fields.</b> Frequency, EIRP, G/T, bit rate, overhead, roll-off, spectral efficiency, "
        "demodulator loss, other attenuations, antenna diameter and link availability belong to the Downlink "
        "and Uplink tabs. The static tabs show them again with a <b>blue label</b>: they are the same fields, so "
        "editing one changes the other."))
    story.append(P(
        "A button that needs to be pressed again because its inputs changed is shown in <b>blue</b> "
        "(Start/Refresh Analysis, Calculate UL Link Budget, Calculate in the static tabs)."))

    # 4 -------------------------------------------------------------------
    story += H1("Orbit & Passes tab")
    story += H2("Satellite and project")
    story.append(P(
        "<b>Satellite</b> lists the satellites of <b>satellites.json</b> (EUMETSAT missions; the file can be "
        "edited). Each satellite has a folder <b>missions/&lt;name&gt;</b> next to the program:"))
    story += BUL([
        "<b>parameters.json</b>: the link budget parameters of the mission (all the fields of the Downlink, "
        "Uplink and Static LB tabs). When it is present they are loaded and the link budget is available.",
        "Without parameters.json only the <b>contact analysis</b> is done: the passes with their times, "
        "duration and maximum elevation. The Downlink, Uplink and Static LB tabs are disabled. To add the "
        "link budget, copy the parameters.json of a similar mission (e.g. AWS) into the folder and edit it.",
        "<b>antenna_pattern.txt</b> (optional): the satellite antenna pattern of the mission (two columns, "
        "angle from nadir in degrees and gain in dBi, separated by commas, semicolons, tabs or spaces; header "
        "lines and # comments are skipped); otherwise Antenna Pattern Default.txt is used. The line under the "
        "pattern buttons and the title of Show Antenna Gain give the file in use; Load Antenna Pattern… only "
        "replaces it for the session.",
        "<b>tle.txt</b>: the last TLE downloaded (or, for a custom orbit, the last synthetic TLE generated). "
        "On selection the TLE is downloaded again when it is older "
        "than one day; without internet the saved one is used.",
    ])
    story.append(P(
        "A mission not launched yet (e.g. EPS-Sterna, of which AWS is the protoflight) has "
        "<b>custom orbit</b>, given by its mean elements: <b>altitude</b> (mean semi-major axis minus the "
        "equatorial radius), <b>eccentricity</b>, argument of perigee <b>ω</b> (90° for a frozen orbit), "
        "<b>inclination</b> (empty: computed for a sun-synchronous orbit) and <b>LTAN</b> (local time of the "
        "ascending node; for a constellation the menu lists the orbital planes). The program builds a "
        "synthetic TLE with its epoch at the analysed date (also with <b>Generate TLE</b>), with the mean "
        "motion set so that the nodal period is that of the J2 theory: the pass times do not drift over the "
        "day. Passes and link budgets are then computed as for a launched satellite. EPS-Sterna uses its "
        "reference orbit: 595.6 km, e = 0.001214, ω = 90°, i = 97.794°, planes at 15:30, 19:30 and 23:30 "
        "(nodal period 96.716 min, altitude over the ellipsoid 600–628 km). SGP4 adds its short-period terms "
        "to the mean elements, so the instantaneous altitude swings a few km more than in the reference "
        "table; the mean altitude and the pass times are those of the reference orbit. The orbit is saved in "
        "the orbit section of the mission's parameters.json."))
    story += FIG("custom_orbit", "EPS-Sterna selected: custom orbit and synthetic TLE.")
    story.append(P(
        "<b>Project Name</b> is shown in the window title, in the PDF report (page header and above "
        "the title), in the proposed PDF file name and in the CSV exports. It is read from and saved to the "
        "parameters file (key project_name)."))
    story += H2("TLE and observation settings")
    story += BUL([
        "<b>TLE Line 1 / Line 2</b>: the satellite elements, filled from tle.txt at startup or with "
        "<b>Load TLE from file</b>. The line below shows the file, the TLE epoch and its age; an "
        "epoch older than 7 days is shown in orange, because the pass times drift.",
        "<b>NORAD ID</b> and <b>Download TLE</b>: download the latest TLE of the satellite and save it in "
        "tle.txt. The sources are tried in turn until one answers: SatNOGS (db.satnogs.org, elements from "
        "Space-Track; it only lists satellites with a known radio transmitter) and N2YO (the TLE on the "
        "satellite's page). When both fail, the status line gives the "
        "reason for each source. "
        "The NORAD ID is read from and saved to the parameters file (key norad_id; AWS is 60543). Without "
        "internet access the TLE already loaded is kept: copy the TLE from a web page into tle.txt instead.",
        "<b>Date</b> (YYYY-MM-DD): the day analysed (00:00–24:00 UTC, 30 s step).",
        "<b>Ground Station</b>: the station of the pass analysis, from the catalogue; <b>Load Ground "
        "Stations…</b> reads another catalogue. The line below shows which file is in use.",
        "<b>Import Link Budget Parameters (.json)…</b>: fills all the parameters, the project name and the "
        "static tabs from a parameters file.",
    ])
    story += H2("Pass analysis")
    story.append(P(
        "<b>Start/Refresh Analysis</b> propagates the orbit for the chosen date and station, finds the passes "
        "above the minimum elevation (5°) and computes the downlink and uplink link budgets along each of them. "
        "The passes are listed under <b>Contact Windows (UTC Time)</b> with their direction ((A)scending or "
        "(D)escending). Select a pass: the Downlink and Uplink tabs show its results. The button turns blue "
        "when the TLE, date or station change."))
    story.append(P(
        "<b>Satellite antenna pattern</b>: the gain of the satellite antenna versus the angle from nadir, used "
        "at every time step (chapter 8.4). It is read from <b>Antenna Pattern Default.txt</b> at startup; "
        "<b>Load Antenna Pattern…</b> loads another file and <b>Show Antenna Gain</b> plots it. The line below "
        "the buttons shows which pattern is in use."))

    # 5 -------------------------------------------------------------------
    story += H1("Downlink tab")
    story += FIG("downlink_params", "Downlink tab: parameters and link losses.")
    story += H2("Parameters")
    story += TABLE([
        ["Field", "Meaning"],
        ["EIRP SAT [dBW]", "Satellite EIRP. Used by the pass analysis only when the Static LB · Downlink tab "
                           "uses an EIRP value (not the Tx chain): it is then taken at the antenna pattern peak."],
        ["Frequency [GHz]", "Downlink carrier frequency."],
        ["C/Io [dBHz]", "Carrier to interference density. Leave empty when there is no interferer."],
        ["G/T GS [dB/K]", "Ground station figure of merit."],
        ["Antenna Diameter GS [m]", "Used by the ITU-R attenuation model and for the 3 dB beamwidth."],
        ["Link Availability [%]", "Rain and clouds are taken for p = 100 − availability [% of time]."],
        ["Other Attenuations [dB]", "Any further fixed loss."],
        ["Tropospheric scintillation", "ITU-R P.618 tropospheric scintillation; off by default (valid from 5° "
                                       "of elevation and about 4–20 GHz)."],
        ["Coded Bit Rate [Mbps]", "Channel bit rate of the transmitter (coding included)."],
        ["Overhead (coded / info)", "Coded bit rate / information bit rate (coding × formatting)."],
        ["Demodulator Loss [dB]", "Receiver implementation loss."],
        ["Roll-off Factor, Spectral Efficiency", "Give the occupied bandwidth: coded rate × (1 + roll-off) / "
                                                 "spectral efficiency (QPSK: 2 b/s/Hz)."],
    ], [5, 11.5])
    story += H2("Link losses from the static link budget")
    story.append(P(
        "With <b>Include the losses entered in the Static LB · Downlink tab</b> (on by default) the pass analysis "
        "uses, from that tab: the Tx chain (transmitter power − circuit and VSWR losses) to which the antenna "
        "pattern gain is added, the ionospheric, polarisation and multipath losses, the modulation degradation, "
        "and the formatting factor, so that Eb/No is referred to the transfer frame bit rate as in the static "
        "budget. The panel shows the values applied. The ground station antenna tracks the satellite and has "
        "no pointing loss."))
    story += H2("Results of the selected pass")
    story += FIG("downlink_results", "Plots and table of a pass.")
    story.append(P(
        "The plots show elevation, Eb/No and C/(No+Io), and Doppler shift along the pass; the table gives every "
        "30 s step. <b>Calculate DL Link Budget</b> recomputes the budget with new parameters without repeating "
        "the orbit propagation; <b>Export Table CSV</b> saves the table (with project and pass on the first "
        "lines). Main columns:"))
    story += TABLE([
        ["Column", "Meaning"],
        ["Off Boresight Angle (°)", "Angle from the satellite nadir to the station."],
        ["S/C Antenna Gain (dBi)", "Satellite antenna gain at that angle, from the pattern."],
        ["EIRP (dBW)", "Tx chain + S/C antenna gain (or EIRP SAT at the pattern peak − roll-off)."],
        ["Atmospheric Att (dB)", "ITU-R attenuation (computed at 5° of elevation, a conservative value)."],
        ["Link Losses (dB)", "Ionospheric + polarisation + multipath, from the static tab."],
        ["C/(No+Io) (dBHz)", "Carrier to noise (and interference) density after the demodulator losses."],
        ["Eb/No (dB)", "C/(No+Io) − 10·log10(transfer frame bit rate)."],
    ], [5, 11.5])

    # 6 -------------------------------------------------------------------
    story += H1("Uplink tab")
    story += FIG("uplink", "Uplink tab.")
    story.append(P(
        "Same principle as the Downlink tab, for the ground station to satellite link: EIRP GS, uplink "
        "frequency, G/T SAT, satellite demodulator loss, uplink attenuations and baseband parameters. G/T SAT is "
        "taken at the satellite antenna boresight; at each time step the pattern roll-off towards the station "
        "(peak gain − gain) is subtracted. <b>Calculate UL Link Budget</b> shows the table of the pass selected "
        "in the Orbit &amp; Passes tab."))

    # 7 -------------------------------------------------------------------
    story += H1("Static link budget tabs")
    story += H2("Layout")
    story += FIG("static_geometry", "Static LB · Downlink: toolstrip, sidebar and an input page.")
    story += BUL([
        "<b>Toolstrip</b>: RUN – <b>Calculate</b> (blue when there are no results yet or the inputs changed); "
        "FILE – <b>Load…</b>, <b>Save…</b>, <b>Reset</b>; EXPORT – <b>Export CSV…</b>, <b>Export PDF…</b> "
        "(enabled once results are up to date).",
        "<b>Sidebar</b>: the input pages (INPUTS) and the result pages (RESULTS, available after Calculate). "
        "At the bottom, the verdict and margin of the last calculation, or OUTDATED when an input changed.",
        "An input error opens the page of the field concerned and puts the cursor in it.",
        "<b>Save…</b> writes the inputs of the tab into the parameters file (section fixed_elevation_downlink "
        "or fixed_elevation_uplink), together with the shared values and the project name, keeping the other "
        "keys of the file. Saving into parameters.json makes them load at startup.",
    ])
    story += H2("Input pages")
    story.append(H3("Station & Geometry"))
    story.append(P(
        "Ground station (its own selection, independent of the pass analysis; its location sets the ITU-R "
        "attenuation and the ionospheric zone), elevation angle (the worst case: the minimum operational "
        "elevation) and satellite altitude, which give the slant range; antenna diameter and link "
        "availability (shared); tropospheric scintillation (off by default)."))
    story += FIG("static_data", "Data Channel page with the CCSDS framing.")
    story.append(H3("Data Channel & Requirements"))
    story.append(P(
        "<b>Data framing</b>: <i>Manual</i> (enter overhead and formatting factor) or <i>CCSDS TM</i> "
        "(downlink): from the coded bit rate and the RS(255,223) interleaving depth I, the tool derives the "
        "transfer frame (223·I bytes: 6-byte primary header, optional OCF 4 B and FECF 2 B), the coding factor "
        "(ASM 4 B per codeblock of 255·I bytes, convolutional 1/2), the formatting factor (frame / user data; "
        "with a Source Packet length also the packet headers, CUC data field header and packet error control "
        "are deducted) and the transfer frame and information bit rates, shown in grey boxes. The (i) icons give "
        "the detailed breakdown. Then the required Eb/No and required nominal margin."))
    story += FIG("static_tx", "Transmitter page: Tx chain and antenna gain from the pattern.")
    story.append(H3("Transmitter"))
    story.append(P(
        "<b>EIRP source</b>: <i>EIRP value</i> (shared field) or <i>Tx chain</i>: EIRP = transmitter power − "
        "antenna circuit loss − VSWR mismatch loss + antenna gain. On the downlink, <b>Antenna gain from the "
        "satellite antenna pattern at this elevation</b> reads the gain from the pattern at the angle from nadir "
        "of the station (e.g. 65.4° at 3° and 628 km); the Tx antenna pointing loss is then not used. The EIRP "
        "actually used is shown under the chain. Tx axial ratio and Tx pointing loss complete the page."))
    story += FIG("static_propagation", "Propagation page: ionospheric loss from S4.")
    story.append(H3("Propagation & PFD Limit"))
    story.append(P(
        "<b>Ionospheric loss</b>: <i>Manual</i> or <i>From S4 scintillation index</i> (ITU-R P.531, chapter "
        "8.3); the page shows the station's geomagnetic latitude and ionospheric zone as guidance for S4. Then "
        "other attenuations (shared; values ≥ 1 dB are flagged, as they may count the atmospheric losses "
        "twice), polarisation mismatch loss (entered, or computed from the Tx and Rx axial ratios when empty), "
        "multipath losses and, on the downlink, the PFD limit."))
    story += FIG("static_ul_rx", "Static LB · Uplink, Receiver page: G/T SAT from the pattern.")
    story.append(H3("Receiver"))
    story.append(P(
        "G/T and demodulator loss (shared), modulation degradation, Rx axial ratio, Rx antenna depointing angle "
        "and 3 dB beamwidth (pointing loss 12·(θ/θ3dB)²; beamwidth 70·λ/D from the GS dish when empty). On the "
        "uplink, <b>G/T SAT at the antenna boresight: subtract the pattern roll-off at this elevation</b> "
        "reduces G/T SAT by the pattern roll-off towards the station, as the pass analysis does at each time "
        "step; the Rx depointing angle is then not used."))
    story.append(H3("Tolerances"))
    story.append(P(
        "Optional favourable/adverse deviations (positive dB) of transmitter power, Tx antenna gain, propagation "
        "loss, Rx G/T and receiver degradation, each with a distribution: TRI (triangular, nominal most likely), "
        "UNI (uniform) or GAU (Gaussian with the limits at ±3σ). They give the statistical margins (chapter "
        "8.7)."))
    story += H2("Result pages")
    story += FIG("static_summary", "Summary page.")
    story += BUL([
        "<b>Summary</b>: overall PASS/FAIL, one tile per key result (Eb/No, Data Recovery Margin, MEAN − 3σ, "
        "Worst Case RSS, PFD margin) with clear-sky / rain-faded values and requirement, notes (e.g. elevation "
        "below 5°) and the derived parameters.",
        "<b>Link Budget</b>: the full table, grouped as Transmitter, Propagation, Receiver, Demodulation, "
        "Margins, for clear sky and rain faded (p = 100 − availability).",
        "<b>PFD Check</b> (downlink): PFD at the ground, in 4 kHz, limit and margin.",
        "<b>With Tolerances</b>: each toleranced quantity down to the margin, with nominal, favourable, adverse, "
        "mean, variance and distribution (rain-faded condition).",
    ])
    story += FIG("static_link_budget", "Link Budget page.")
    story += H2("Exports")
    story.append(P(
        "<b>Export CSV…</b> saves inputs and results as a table. <b>Export PDF…</b> creates an A4 report: "
        "summary with PASS/FAIL, input and derived parameters, link budget table, PFD check, link budget with "
        "tolerances and the method and assumptions, with the project name and the software version."))

    # 8 -------------------------------------------------------------------
    story += H1("Calculation method")
    story += H2("Geometry")
    story.append(P(
        "Spherical Earth (R = 6378.137 km), circular orbit at altitude h. For an elevation el, the slant range "
        "is d = √((R+h)² − ((R+h_GS)·cos el)²) − (R+h_GS)·sin el and the angle from the satellite nadir to the "
        "station is η = asin((R+h_GS)·cos el / (R+h)). In the pass analysis elevation, range and angle come from "
        "the orbit propagation (Skyfield)."))
    story += H2("Propagation")
    story += BUL([
        "Free space loss: FSL = 20·log10(d[km]) + 20·log10(f[GHz]) + 92.45 dB.",
        "Atmospheric attenuation (ITU-R P.618 slant path at the station location): gas (P.676), cloud (P.840), "
        "rain (P.618/P.837) exceeded for p = 100 − availability [%], optional tropospheric scintillation. Clear "
        "sky: gas only; rain faded: gas, cloud and rain.",
        "Polarisation mismatch between co-rotating elliptical antennas from the axial ratios (voltage ratios "
        "r = 10^(AR/20)), co-polar value.",
    ])
    story += H2("Ionospheric scintillation")
    story.append(P(
        "From the S4 index measured at a reference frequency (GPS L1 1575.42 MHz by default), as in ESA "
        "AWS-TN-ESA-SYS-0052: S4 is scaled to the link frequency as f^−1.5, the peak-to-peak power fluctuation "
        "is 27.5·S4^1.26 dB and the fade is half of it. The geomagnetic latitude (dipole) gives the zone: "
        "equatorial (&lt; 20°), mid-latitude (&lt; 58°), auroral oval (&lt; 72°), polar cap."))
    story += H2("Antennas")
    story += BUL([
        "Satellite antenna: the pattern gives the absolute gain (dBi) at the angle from nadir η; outside the "
        "table the gain is 0 dBi. Downlink EIRP = Tx chain + G(η); uplink G/T = G/T SAT − (G_peak − G(η)).",
        "Ground station antenna: tracks the satellite; optional depointing loss 12·(θ/θ3dB)², with "
        "θ3dB = 70·λ/D.",
    ])
    story += H2("Carrier to noise and Eb/No")
    story.append(P(
        "Rx power = EIRP − FSL − atmospheric − other − ionospheric − polarisation − multipath − pointing losses. "
        "S/No = Rx power + G/T + 228.6 dBHz; C/No = S/No − demodulator loss − modulation degradation (combined "
        "with C/Io when given). Eb/No = C/No − 10·log10(R_frame), with the transfer frame bit rate R_frame = "
        "coded bit rate / coding factor = coded bit rate × formatting factor / overhead. Data Recovery Margin = "
        "Eb/No − required Eb/No."))
    story += H2("Power flux density")
    story.append(P(
        "PFD = EIRP − 10·log10(4π·d²) [dBW/m²], without atmospheric loss (conservative). Referred to the 4 kHz "
        "bandwidth of the ITU limits: PFD(4 kHz) = PFD − 10·log10(B / 4 kHz), with the occupied bandwidth B. "
        "PFD margin = limit − PFD(4 kHz) ≥ 0 dB (ITU Radio Regulations Art. 21)."))
    story += H2("Statistical margins")
    story.append(P(
        "Each tolerance with adverse limit a = −adverse and favourable limit b = +favourable contributes a mean "
        "and a variance: TRI mean (a+b)/3, variance (a²+b²−ab)/18; UNI mean (a+b)/2, variance (b−a)²/12; GAU "
        "mean (a+b)/2, variance ((b−a)/6)². Contributions are independent: means and variances add. MEAN − 3σ "
        "= nominal margin + Σ mean − 3·√Σ variance; Worst Case RSS = nominal margin − √Σ adverse²."))

    # 9 -------------------------------------------------------------------
    story += H1("Data files")
    story += H2("parameters.json")
    story.append(P(
        "A JSON file with the values of the main tabs and two optional sections with the inputs of the static "
        "tabs (fixed_elevation_downlink, fixed_elevation_uplink). Missing keys are ignored; a null value leaves "
        "the field empty (e.g. c_io_dbhz: no interferer). Keys of the main tabs:"))
    keys = ["project_name"] + parameter_keys()
    rows = [["Keys"]] + [[", ".join(keys[i:i + 4])] for i in range(0, len(keys), 4)]
    story += TABLE(rows, [16.5], mono_cols=(0,), bold_cols=())
    story.append(P(
        "The static sections are best written with <b>Save…</b> in the static tabs rather than by hand."))
    story += H2("Other files")
    story += TABLE([
        ["File", "Format"],
        ["tle.txt", "The two TLE lines (an optional name line before them is accepted)."],
        ["ground_stations.txt", "One station per line: name, latitude [deg], longitude [deg], altitude [m]; "
                                "lines starting with # are comments."],
        ["Antenna Pattern Default.txt", "One point per line: angle from nadir [deg], gain [dBi], comma "
                                        "separated, from the most negative to the most positive angle. It "
                                        "should cover the lowest elevation used (66° from nadir at 3° and "
                                        "628 km)."],
    ], [5, 11.5], mono_cols=(0,), bold_cols=())

    # 10 ------------------------------------------------------------------
    story += H1("Troubleshooting")
    story += TABLE([
        ["Symptom", "Cause and remedy"],
        ["The program stops with a licence message", "See 2.2."],
        ["Fields empty at startup", "parameters.json or tle.txt are not next to the executable (or not valid): "
                                    "copy them there, or load them from the Orbit &amp; Passes tab."],
        ["'Ground stations: built-in list'", "ground_stations.txt was not found: the built-in stations are "
                                              "used."],
        ["'Pattern: built-in example'", "Antenna Pattern Default.txt was not found: a generic example pattern "
                                         "is used."],
        ["Input Error message", "A field is empty or not a number: the static tabs open its page and select "
                                "it."],
        ["Static results marked OUTDATED", "An input changed after Calculate: press Calculate again (exports "
                                           "are disabled until then)."],
        ["No passes listed", "No pass above 5° on that date, or the TLE is old: update tle.txt."],
    ], [5.5, 11])
    return story


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    doc = ManualDoc(OUT, pagesize=A4, leftMargin=2 * cm, rightMargin=2 * cm, topMargin=2 * cm,
                    bottomMargin=2.2 * cm, title=f"{TITLE} – {SUBTITLE}", author=APP_AUTHOR,
                    subject=f"Version {APP_VERSION}")
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id="body")
    doc.addPageTemplates([
        PageTemplate(id="cover", frames=[frame], onPage=_cover),
        PageTemplate(id="page", frames=[frame], onPage=_page),
    ])
    doc.multiBuild(content())
    print(f"Manual written to {OUT}")


if __name__ == "__main__":
    main()
