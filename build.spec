# PyInstaller spec for Katana Splicing Tool
# Build example:
#   pyinstaller build.spec

block_cipher = None

datas = [
    ('configs', 'configs'),
    ('assets', 'assets'),
    ('software/Jutils', 'software/Jutils'),
    ('software/SUPPA', 'software/SUPPA'),
    ('software/rmats2sashimiplot/src', 'software/rmats2sashimiplot/src'),
    ('software/rmats2sashimiplot/pyproject.toml', 'software/rmats2sashimiplot'),
    ('software/rmats2sashimiplot/setup.py', 'software/rmats2sashimiplot'),
]

a = Analysis(
    ['app/main.py'],
    pathex=[],
    binaries=[],
    datas=datas,
    hiddenimports=[
        'PySide6.QtCore',
        'PySide6.QtGui',
        'PySide6.QtWidgets',
        'matplotlib.backends.backend_qtagg',
        'seaborn',
        'pandas',
        'numpy',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        'tkinter',
        'PyQt5',
        'PyQt6',
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='KatanaSplicingTool',
    icon='assets/katana_icon.ico',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name='KatanaSplicingTool',
)
