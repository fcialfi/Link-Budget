# -*- mode: python ; coding: utf-8 -*-

from PyInstaller.utils.hooks import collect_data_files

# Include tutti i file di dati della libreria itur
itur_datas = collect_data_files('itur', includes=['data/**/*'])

a = Analysis(
    ['gui.py'],
    pathex=[],
    binaries=[],
    datas=itur_datas + [
        ('calculations.py', '.'),
        ('ground_stations.txt', '.'),
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

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='gui',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,  # Metti a False se vuoi nascondere la console
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=['Satellite.ico'],  # Se non serve, puoi anche rimuoverla
)
