# TI-84 Plus CE Chemistry Apps

This repository contains three CEdev-based C programs for the TI-84 Plus CE:

- `lewis-dot`: Interactive Lewis dot structure generator with support for molecular charge and resonance forms.
- `vsepr-table`: Scrollable VSEPR reference table (geometry, shape, and hybridization).
- `formula-sheet`: Scrollable EM spectrum wavelength chart.

## Requirements

- [CEdev](https://ce-programming.github.io/toolchain/) installed
- `CEDEV` environment variable set to your CEdev install path
- `make` available in your shell

Each project `Makefile` includes:

```make
ifndef CEDEV
$(error CEDEV environment variable is not set.)
endif
```

## Build

Build from each app directory:

```powershell
cd lewis-dot
make
```

```powershell
cd vsepr-table
make
```

```powershell
cd formula-sheet
make
```

Compiled calculator files are written to each app's `bin/` folder.

To clean generated files:

```powershell
make clean
```

## Controls

### `lewis-dot`

- Arrow keys: move periodic-table cursor
- `Enter`: add highlighted element
- `Del`: remove last atom
- `Alpha`: cycle charge (`0`, `+1`, `+2`, `-1`, `-2`)
- `2nd`: generate Lewis structure
- Left/Right: cycle resonance structures
- `Clear`: return to periodic table view
- `Mode`: quit

### `vsepr-table`

- Arrow keys: scroll table
- `Clear`: quit

### `formula-sheet`

- Up/Down: scroll chart
- `Clear`: quit

## Data Generation (VSEPR)

`vsepr-table/src/vsepr_data.h` can be regenerated from CSV:

```powershell
cd vsepr-table
python convert.py vsepr_data.csv src/vsepr_data.h
```

## Project Layout

- `lewis-dot/` - Lewis structure generator source + CEdev build files
- `vsepr-table/` - VSEPR table source/data + CEdev build files
- `formula-sheet/` - formula sheet source + CEdev build files
