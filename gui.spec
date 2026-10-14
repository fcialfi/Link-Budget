# -*- mode: python ; coding: utf-8 -*-

from PyInstaller.utils.hooks import collect_data_files

# Include tutti i file di dati della libreria itur
itur_datas = collect_data_files('itur', includes=['data/**/*'])

a = Analysis(
    ['gui.py'],
    # The project folder is also a package (__init__.py), so without it on
    # the search path PyInstaller does not find calculations.py,
    # static_report.py and licensing.py as modules, and misses what they
    # import (reportlab, cryptography).
    pathex=[SPECPATH],
    binaries=[],
    datas=itur_datas + [
        ('calculations.py', '.'),
        ('static_report.py', '.'),
        ('licensing.py', '.'),
        # Logo and icon images; also the licence public key
        # (assets/license_public_key.pem) once license_tool.py keygen has
        # created it, which turns the licence check on (see README).
        ('assets', 'assets'),
    ],
    hiddenimports=[
        'itur',
        # SciPy special functions with compiled extensions that PyInstaller
        # sometimes misses without an explicit hint.
        'scipy.special._ufuncs',
        'scipy.special._ufuncs_cxx',
        # Ensure astropy's optional test runner import can resolve and bring
        # along its dependency on pytest so the bundled executable imports
        # without errors.
        'astropy.tests.runner',
        'pytest',
        # Local modules and the libraries they import only inside functions
        # (PDF report, licence check), listed so a build never leaves them out.
        'calculations',
        'static_report',
        'licensing',
        'version',
        'tle_source',
        'reportlab',
        'cryptography.exceptions',
        'cryptography.hazmat.primitives.serialization',
        'cryptography.hazmat.primitives.asymmetric.ed25519',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        # Avoid optional test-only modules that trigger warnings when they are
        # not present in the runtime environment.
        'scipy.special._cdflib',
    ],
    noarchive=False,
    optimize=0,
)

pyz = PYZ(a.pure)

# Build mode. ONEFILE = False (default) builds a folder, dist/Link-Budget-Tool/,
# with Link-Budget-Tool.exe next to its libraries (in _internal/): it starts
# in seconds because nothing has to be unpacked, and is distributed as a zip
# of the whole folder. ONEFILE = True builds a single .exe that unpacks a few
# hundred MB (scipy, astropy, matplotlib, pandas, the ITU-R data) to a temp
# folder at every start -- minutes with an antivirus scanning the files.
ONEFILE = False
NAME = 'Link-Budget-Tool'

# Splash screen, shown by the bootloader as soon as the .exe is started,
# before any Python code can run (in one-file mode, while the bundle is
# being unpacked). gui.py closes it as soon as its own loading window is up.
splash = Splash(
    'assets/splash.png',
    binaries=a.binaries,
    datas=a.datas,
    text_pos=(40, 262),
    text_size=10,
    text_color='#b8c4d0',
    minify_script=True,
    always_on_top=False,
)

_exe_options = dict(
    name=NAME,
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    # UPX-compressed libraries have to be decompressed at every start, and
    # antivirus programs often scan them as suspicious: off.
    upx=False,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=['assets/exostaff_icon.ico'],
)

if ONEFILE:
    exe = EXE(pyz, a.scripts, splash, splash.binaries, a.binaries, a.datas, [], runtime_tmpdir=None, **_exe_options)
else:
    exe = EXE(pyz, a.scripts, splash, [], exclude_binaries=True, **_exe_options)
    coll = COLLECT(exe, splash.binaries, a.binaries, a.datas, strip=False, upx=False, name=NAME)
