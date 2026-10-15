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
ISSUE = "Issue 1.4"

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
        ["ground_stations.csv", "Ground station catalogue (CSV with header, opens in Excel): name, latitude, longitude, "
                                "altitude (m); optional antenna diameter (m), own G/T, its elevation, Tsys and number "
                                "of antennas. Edit Stations... edits it as a table."],
        ["satellites.json", "The satellites of the Satellite menu (chapter 4.1)."],
        ["missions\\<name>\\", "Per mission: parameters.json (project name, link budget parameters, static inputs; chapter "
                                 "9), antenna_pattern.txt and tle.txt (chapter 4.1)."],
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
        "Uplink and Static LB tabs). When it is present they are loaded and the link budget is available. "
        "<b>Mission Parameters…</b>, next to the Satellite menu, shows the whole file as a table grouped by "
        "subject (mission, orbit, downlink, ground station, uplink, statistics, static link budgets), with "
        "units and descriptions: double-click a value to change it (an empty value removes it), Save writes "
        "the file and loads the values in the window.",
        "Without parameters.json only the <b>contact analysis</b> is done: the passes with their times, "
        "duration and maximum elevation. The Downlink, Uplink and Static LB tabs are disabled. To add the "
        "link budget, copy the parameters.json of a similar mission (e.g. AWS) into the folder and edit it.",
        "<b>antenna_pattern.txt</b> (optional; or any .txt/.csv file with \"pattern\" in its name): the "
        "satellite antenna pattern of the mission (two columns, angle from nadir in degrees and gain in dBi, "
        "separated by commas, tabs, spaces or semicolons, then also with a decimal comma as Excel writes it with "
        "Italian settings; header lines and # comments are skipped). There is no common default: a mission "
        "without a readable pattern gets an isotropic 0 dBi antenna, reported in orange. The line under the "
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
        "<b>inclination</b> (empty: computed for a sun-synchronous orbit), on the <b>Custom orbit</b> row. "
        "The <b>Constellation</b> row gives the orbital <b>planes</b> by the local time of their ascending "
        "node (LTAN, hh:mm separated by commas, e.g. 15:30, 19:30, 23:30), the satellites per plane "
        "(<b>Sats/plane</b>, evenly spaced in the plane) and the <b>Phasing</b> (shift of the argument of "
        "latitude from one plane to the next). <b>Analyse: plane / sat #</b> chooses the satellite used by "
        "the pass analysis and the link budgets. The program builds a "
        "synthetic TLE with its epoch at the analysed date (also with <b>Generate TLE</b>), with the mean "
        "motion set so that the nodal period is that of the J2 theory: the pass times do not drift over the "
        "day. Passes and link budgets are then computed as for a launched satellite. EPS-Sterna uses its "
        "reference orbit: 595.6 km, e = 0.001214, ω = 90°, i = 97.794°, planes at 15:30, 19:30 and 23:30 "
        "(nodal period 96.716 min, altitude over the ellipsoid 600–628 km). SGP4 adds its short-period terms "
        "to the mean elements, so the instantaneous altitude swings a few km more than in the reference "
        "table; the mean altitude and the pass times are those of the reference orbit. The orbit is saved in "
        "the orbit section of the mission's parameters.json with <b>Save</b> (Static LB tabs); the "
        "Statistics tab can use the whole constellation. To change the orbits: edit the fields (a new plane "
        "list or Sats/plane updates the plane and satellite menus), check them with <b>Show "
        "Constellation</b>, then Save."))
    story += FIG("custom_orbit", "EPS-Sterna selected: custom orbit, constellation and synthetic TLE.")
    story.append(P(
        "<b>Show Constellation</b> lists every satellite (plane, RAAN, argument of latitude and inclination "
        "at 00:00 UTC of the analysed date) and draws the ground tracks over the first 100 minutes of the day "
        "with the ground stations, and the satellites in their planes (RAAN against argument of latitude). "
        "<b>3D View</b> draws the orbits around the Earth (inertial frame), the satellites and the stations; "
        "drag to rotate, and move the slider to follow the satellites over the first 100 minutes of the day."))
    story += FIG("constellation", "Show Constellation: 3 planes × 2 satellites, phasing 60°.")
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
        "<b>Ground Station</b>: the station of the pass analysis, from the catalogue (ground_stations.csv); "
        "<b>Load Ground Stations…</b> reads another catalogue and <b>Edit Stations…</b> shows the catalogue "
        "as a table to add, change or delete stations (position, diameter, own G/T and Tsys, antennas) and "
        "save it. The line below shows which file is in use.",
        "<b>Min Elevation [deg]</b> (default 5°) and the station's <b>horizon mask</b>: a sample is in "
        "contact only above both. The mask is the file horizon_masks\\&lt;station&gt;.txt next to the program "
        "(two columns: azimuth from North and minimum elevation, linear between the points; see "
        "horizon_masks\\README.txt); <b>Horizon Mask</b> plots it, and the line under the station gives the "
        "file in use. Without a file only the minimum elevation applies.",
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
        "at every time step (chapter 8.4). It is read from the folder of the selected mission; "
        "<b>Load Antenna Pattern…</b> loads another file for the session and <b>Show Antenna Gain</b> plots it. The line below "
        "the buttons shows which pattern is in use."))

    story += H1("Statistics tab")
    story += FIG("statistics", "Statistics over 7 days: summary and data per day.")
    story.append(P(
        "The <b>Statistics</b> tab adds up the pass analysis over several days and ground stations, for the "
        "satellite, minimum elevation and horizon masks of the Orbit &amp; Passes tab and the link budget of the "
        "Downlink tab (same models: antenna pattern, atmospheric attenuation and G/T at each elevation, losses "
        "of the Static LB · Downlink tab). Enter the <b>Start Date</b> and the number of <b>Days</b>, select one "
        "or more <b>Ground Stations</b> (about a second per day and station) and press <b>Compute "
        "Statistics</b>."))
    story += BUL([
        "<b>Station network</b>: the satellite is in contact when any selected station sees it (above the "
        "minimum elevation and that station's horizon mask); at each instant the station with the best margin "
        "is used. Each station uses its own antenna diameter (G/T scaling) and mask; the pass table gives the "
        "stations of each contact.",
        "<b>Usable minutes</b>: samples in contact whose margin, Eb/No minus the <b>Required Eb/No</b>, is at "
        "least the <b>Required Margin</b>; when these fields are empty the values of the Static LB · Downlink "
        "tab are used.",
        "<b>Capacity</b>: the information bit rate (coded bit rate / overhead) over the usable minutes, in MB "
        "of 10^6 bytes: what the link could download.",
        "<b>On-board data</b>: the rates in kbps of what is sent on this downlink band, <b>Instrument Data</b>, "
        "<b>HKTM</b> and <b>Other Data</b> (the uplink, telecommands, is not included). With them the on-board "
        "memory is simulated, starting empty: it fills at the total rate and empties at the downlink rate "
        "during the usable minutes. This gives the volume actually <b>downloaded</b>, the data <b>on board</b> "
        "at the end of each day and its maximum, and the <b>data latency</b> (time from generation to "
        "arrival on ground, first in first out; mean and maximum). The downlink keeps up when the mean "
        "capacity is at least the data generated per day. The rates are read from and saved to the parameters "
        "file (downlink_data_kbps).",
        "<b>Latency Requirement</b> (min) and <b>Within</b> (% of data), with the on-board data rates: e.g. "
        "60 min for 90 % of the data. The summary gives the share of the downloaded data on ground within the "
        "required time, the time within which the required share arrives, and whether the requirement is met. "
        "Saved in the parameters file (latency_requirement).",
        "<b>Real-time TM</b> (kbps): telemetry sent live during the passes and not stored (e.g. real-time "
        "HKTM); it takes that part of the downlink bit rate, leaving less for the stored data.",
        "<b>On-board Memory</b> (GB, with the on-board data rates): memory of each satellite for these data, "
        "e.g. 4 GB for one partition. The summary gives the largest share used, the hours without contacts "
        "until it is full and, if it fills, the data lost (the new data are discarded; the latency of the rest "
        "is then approximate). Saved in the parameters file (memory_gb: 4 GB, the Mission Data partition, "
        "for EPS-Sterna; 64 GB per band in total); lbtool.py --memory-gb.",
        "<b>Latency at end of pass</b>: the data downloaded during a pass count as on ground at the end of that "
        "pass, when the dump is complete and a product can be processed (the convention of product timeliness "
        "requirements); unchecked, each bit counts when it comes down, a few minutes earlier. Saved in the "
        "parameters file (latency_at_pass_end, on for EPS-Sterna); lbtool.py --latency-at-pass-end.",
        "<b>Save to Mission</b> writes these inputs (data rates, latency requirement, memory, minimum pass, "
        "antenna setup, latency at end of pass, ground stations) into the parameters.json of the selected "
        "mission, keeping the rest of the file; they are loaded again whenever the mission is selected. "
        "<b>Reload</b> reads them back; the line below the buttons names the file.",
        "<b>Min Pass Duration</b> (min): contacts (usable minutes of one satellite over one station, across "
        "midnight) shorter than this are not used, being too short to schedule (acquisition, lock, start of the "
        "dump). Empty: 0, every contact is used. Saved in the parameters file (min_pass_min); lbtool.py "
        "--min-pass.",
        "<b>Antenna Setup</b> (min, whole constellation): time to configure an antenna for another satellite. "
        "After serving a satellite the antenna cannot take a different one for this long (the same satellite "
        "resumes at once); a satellite waiting for it counts as a conflict. Saved in the parameters file "
        "(antenna_setup_min; 0 by default: the EPS-Sterna requirement gives at most 7 min to configure a "
        "station for a pass from scratch, a worst case; up to about 3 min it has no effect on SVL+PUN).",
        "<b>Gaps</b>: time from the end of a contact to the start of the next one; a gap counts for the day the "
        "next contact starts.",
        "For a satellite without link budget parameters only contacts, contact minutes and gaps are given.",
        "<b>Whole constellation</b> (custom-orbit missions with several satellites): all the satellites of the "
        "planes (LTAN list) × <b>Sats/plane</b>, evenly spaced in each plane and shifted by <b>Phasing</b> from "
        "one plane to the next (custom orbit fields of the Orbit &amp; Passes tab, saved in the orbit section). "
        "Each satellite has its own memory with the same data rates. Each station tracks one satellite at a "
        "time and keeps it while the link closes; free stations go to the satellites with more data on board. "
        "A satellite ready to download whose stations are all busy waits (<b>conflict</b>). The summary gives "
        "the system figures (served minutes, conflicts, volume, memory, latency and requirement), the Per "
        "Satellite table the figures of each satellite, the Per Station table the use of each station "
        "(served minutes, volume, busy time, minutes with more than one satellite usable).",
        "The <b>Per Station</b> table is also given for one satellite and several stations. Normally one "
        "antenna per site is used, the others being backups; with <b>Use backup antennas</b> (lbtool.py: "
        "--backup-antennas) a station serves as many satellites at a time as its antennas (last column of "
        "ground_stations.csv). Busy is the share of the antenna time in use, Contention the minutes with more "
        "satellites usable than antennas.",
    ])
    story += FIG("statistics_plots", "Statistics plots: AWS over Svalbard, 104 kbps on board, latency "
                 "requirement 90 % within 100 min.")
    story.append(P(
        "Three plots follow the summary. <b>Contacts</b>: a timeline, one row per day (for a constellation, one "
        "row per satellite over the first day), each pass drawn from its start to its end in UTC and coloured by "
        "its station (grey: no usable minute); the gaps are the white spaces. <b>Data on board</b>: the memory "
        "occupation along the whole period (all satellites for a constellation), with its maximum. <b>Data "
        "latency</b>: the share of the data on ground within a given latency, with the requirement and whether "
        "it is met. Without on-board data rates these two plots become the <b>gaps between contacts</b> of each "
        "satellite (its revisit time) and the <b>contact time by station</b> of each satellite (or day), with "
        "the minutes spent waiting for a busy station."))
    story += FIG("statistics_tables", "Per-day and per-pass tables.")
    story.append(P(
        "The tables give each day (contacts, contact and usable minutes, capacity, downloaded volume, longest "
        "gap, data on board) and each pass (times, stations, maximum elevation, minutes, capacity, downloaded "
        "volume, minimum margin; a pass with no usable minute is shown in orange). <b>Export CSV…</b> writes "
        "the summary and both tables. <b>Export PDF…</b> writes an A4 report: the verdict on the latency "
        "requirement, the inputs used, the summary, the plots and the tables per station, satellite, day and "
        "pass (the first 400 passes; all of them are in the CSV)."))

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
        ["G/T GS [dB/K]", "Ground station figure of merit of the reference antenna, at the reference "
                          "elevation in clear sky (Svalbard 3 m: 9.9 dB/K)."],
        ["Antenna Diameter GS [m]", "Used for a station without a diameter in ground_stations.csv (5th "
                                    "column): ITU-R attenuation model, 3 dB beamwidth and G/T scaling."],
        ["G/T Reference Diameter [m]", "Diameter of the antenna the G/T GS refers to. A station with another "
                                       "diameter gets G/T + 20 log10(D / Dref) (same feed and receiver: a "
                                       "2.4 m antenna is 1.94 dB below a 3 m one). Empty: no scaling."],
        ["G/T Reference Elevation [deg], System Noise Temperature [K]",
         "With a system noise temperature (at the reference elevation, clear sky) the G/T follows the sky "
         "noise: T_sky = 275 K × (1 − 10^(−A/10)) + 2.7 K × 10^(−A/10) from the ITU-R attenuation A of the "
         "actual elevation and weather (ITU-R P.618), and G/T = G/T_ref − 10 log10((Tsys − T_sky,ref + "
         "T_sky) / Tsys). At low elevation and in rain the G/T is lower. Ground spill-over is not modelled. "
         "Empty: constant G/T."],
        ["Radome Loss [dB]", "G/T degradation of the radome (e.g. 0.25 dB), subtracted from the G/T."],
        ["Solar Flux [sfu]", "Solar flux density at the downlink frequency, for the Sun noise (1 sfu = "
                             "1e-22 W/m²/Hz). Typically 100–250 sfu at L/S band at solar maximum. Empty: quiet "
                             "Sun, T = 120000 × f^−0.75 K (about 80 000 K, 49 sfu, at 1.7 GHz)."],
        ["Link Availability [%]", "Rain and clouds are taken for p = 100 − availability [% of time]."],
        ["Other Attenuations [dB]", "Any further fixed loss."],
        ["Tropospheric scintillation", "ITU-R P.618 tropospheric scintillation; off by default (valid from 5° "
                                       "of elevation and about 4–20 GHz)."],
        ["Sun noise (Sun outage)", "On by default, with a System Noise Temperature. When the Sun is near the "
                                   "ground antenna beam it adds T_sun × Ω_sun × G(θ) / 4π to the system noise "
                                   "temperature (θ: angle Sun–satellite seen from the station; Gaussian main "
                                   "beam of 70 λ/D, peak gain G/T + 10 log10 Tsys; ITU-R S.465 side lobes "
                                   "32 − 25 log10 θ outside it). G/T and Eb/No drop by 10 log10(1 + ΔT/Tsys). "
                                   "Sun position: Astronomical Almanac low-precision formulae (about 0.01°). "
                                   "A 3 m antenna at 1.7 GHz with the quiet Sun on its axis: +510 K, about "
                                   "−6 dB with Tsys 150 K. The contact list marks the passes with at least "
                                   "0.5 dB (☀ SUN NOISE)."],
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
        "lines). The table uses short headings (hover a heading for its full name and unit) and the "
        "<b>Columns</b> switches above it show or hide the Geometry, Link Budget and PFD groups; samples above "
        "the PFD limit are in red. Main columns (full names, as in the CSV export):"))
    story += TABLE([
        ["Column", "Meaning"],
        ["Off Boresight Angle (°)", "Angle from the satellite nadir to the station."],
        ["S/C Antenna Gain (dBi)", "Satellite antenna gain at that angle, from the pattern."],
        ["EIRP (dBW)", "Tx chain + S/C antenna gain (or EIRP SAT at the pattern peak − roll-off)."],
        ["Atmospheric Att (dB)", "ITU-R attenuation at the elevation of the sample (the model is evaluated "
                                 "on an elevation grid and interpolated)."],
        ["G/T GS (dB/K)", "Ground station G/T at the elevation of the sample, with the Sun noise (see the "
                          "parameters)."],
        ["Sun Separation (°)", "Angle between the Sun and the satellite seen from the station."],
        ["Doppler Rate (Hz/s)", "Rate of change of the downlink Doppler shift, −f r''/c, from the range "
                                "acceleration (relative velocity, satellite gravity, station rotation): exact at "
                                "every sample, peak at the closest approach. Shown dashed in the Doppler plot; "
                                "the contact list gives the largest Doppler shift and rate of each pass."],
        ["Sun Noise (dB)", "G/T degradation from the Sun in the antenna beam."],
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
        "multipath losses and, on the downlink, the PFD limits: <b>PFD Limit 0–5° El</b> and, optionally, "
        "<b>PFD Limit 25–90° El</b>, linear in between as in the ITU Radio Regulations Art. 21 (without the "
        "second one the first applies at every elevation). The static check uses the limit at its elevation; "
        "the pass analysis checks the PFD at every sample against the same limits."))
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
        "PFD margin = limit − PFD(4 kHz) ≥ 0 dB (ITU Radio Regulations Art. 21). The limit is the 0–5° value up "
        "to 5° of elevation, the 25–90° value from 25°, linear in between. Along the passes the PFD is "
        "computed at every sample with the EIRP of the pattern towards the station; the contact list gives the "
        "minimum PFD margin of each pass and flags a pass that exceeds the limit."))
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
        ["ground_stations.csv", "CSV with the header name, latitude_deg, longitude_deg, altitude_m (required) and the "
                                "optional columns (cells may be left empty): antenna_diameter_m; gt_dbk, the "
                                "station's own measured G/T [dB/K, clear sky], which replaces the mission G/T GS "
                                "scaled with the diameter; gt_elevation_deg, the elevation of that G/T (empty: the "
                                "mission's G/T Reference Elevation, usually 5°: G/T values are normally given at 5°, "
                                "clear sky; if it was measured at zenith instead, the low-elevation G/T is about "
                                "0.8 dB optimistic at L band); tsys_k (empty: the mission's); antennas (default 1), "
                                "backups included: used by the statistics only with Use backup antennas. Excel's "
                                "semicolon-separated CSV is accepted; the old ground_stations.csv without header "
                                "is still read. Edit Stations... (Orbit & Passes) shows and edits it as a table."],
        ["antenna_pattern.txt", "One point per line: angle from nadir [deg], gain [dBi], separated by a "
                                "comma, tabs, spaces or a semicolon. It should cover the lowest elevation "
                                "used (66° from nadir at 3° and 628 km)."],
    ], [5, 11.5], mono_cols=(0,), bold_cols=())

    story += H1("Command line")
    story.append(P(
        "<b>lbtool.py</b> runs the pass analysis and the statistics without the window, from the sources "
        "(Python environment of the project), for the missions of satellites.json and with the same "
        "calculations as the window. Examples:"))
    story += BUL([
        "<font face='Courier'>python lbtool.py missions</font>: lists the missions.",
        "<font face='Courier'>python lbtool.py passes --mission AWS --station Svalbard --date 2026-10-14 "
        "-o passes.csv</font>: contacts of the day with maximum elevation, minimum Eb/No and PFD margin; the "
        "samples above the horizon go to the CSV file.",
        "<font face='Courier'>python lbtool.py statistics --mission AWS --stations Svalbard,Kangerlussuaq "
        "--start 2026-10-14 --days 30 --instruments 100 --hktm 4 -o stats.csv</font>: the statistics of the "
        "Statistics tab over a station network, with the on-board data rates in kbps.",
        "Options: <font face='Courier'>--download</font> (download the TLE first), "
        "<font face='Courier'>--plane hh:mm</font> (orbital plane of a custom-orbit mission), "
        "<font face='Courier'>--min-elevation</font>, <font face='Courier'>--required-ebno</font>, "
        "<font face='Courier'>--required-margin</font>, <font face='Courier'>--constellation</font>, <font face='Courier'>--latency</font> / <font face='Courier'>--latency-pct</font>; <font face='Courier'>-h</font> lists them all.",
    ])

    # 10 ------------------------------------------------------------------
    story += H1("Troubleshooting")
    story += TABLE([
        ["Symptom", "Cause and remedy"],
        ["The program stops with a licence message", "See 2.2."],
        ["Fields empty at startup", "parameters.json or tle.txt are not next to the executable (or not valid): "
                                    "copy them there, or load them from the Orbit &amp; Passes tab."],
        ["'Ground stations: built-in list'", "ground_stations.csv was not found: the built-in stations are "
                                              "used."],
        ["'Pattern: isotropic 0 dBi' (orange)", "The mission folder has no readable antenna pattern: put "
                                                 "an antenna_pattern.txt there."],
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
