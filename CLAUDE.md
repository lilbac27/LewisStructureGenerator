# TI-84 Plus CE Projects

## Project Overview

This is a multi-program workspace for the TI-84 Plus CE graphing calculator, built with C and the CE C/C++ Toolchain. Each program lives in its own subdirectory and compiles to a `.8xp` file that runs on the calculator or in the CEmu emulator.

## Local Environment

- **OS:** Windows 11, shell is Git Bash (MSYS2/MINGW64)
- **Toolchain:** CE C/C++ Toolchain v14.2, installed at `~/CEdev` (`C:\Users\bacon\CEdev`)
- **`CEDEV` and `PATH`** are set in `~/.bashrc` — new terminals pick them up automatically
- **`make`** is installed at `~/bin/make.exe` (Cosmopolitan GNU Make 4.4.1)
- **No Python available** in this shell — use `node` for scripting if needed

## Build System

- **Build a single project:** `make -C <project-dir>` (e.g., `make -C vsepr-table`)
- **Clean a single project:** `make -C <project-dir> clean`
- **Build all projects:** use the "Build All" VSCode task, or loop: `for dir in */; do [ -f "$dir/Makefile" ] && make -C "$dir"; done`
- **If build fails with stale paths:** run `make -C <project-dir> clean` first — dependency files may reference another machine's paths
- **Toolchain docs:** https://ce-programming.github.io/toolchain/

## Hardware Constraints

These constraints are non-negotiable and must be respected in ALL code:

- **CPU:** Zilog eZ80 @ 48 MHz — this is a slow processor, optimize accordingly
- **RAM:** 256 KB total, ~150 KB usable for programs — keep memory usage minimal
- **Screen:** 320x240 pixels, 8-bit indexed color (256-color palette)
- **Input:** Calculator keypad only (no touch, no mouse)
- **No floating point hardware** — never use `float` or `double`; use integer math exclusively

## Coding Rules

### Mandatory

1. **Integer math only.** Never use `float` or `double`. Floating point is software-emulated on the eZ80 and extremely slow. Use fixed-point integer arithmetic if fractional values are needed.
2. **Double buffer all rendering.** Always use `gfx_SetDrawBuffer()` at init and `gfx_SwapDraw()` at the end of each frame to prevent flickering.
3. **Use `kb_Scan()` + `kb_Data[]` for real-time input** in game loops. Only use `os_GetCSC()` for menus or blocking input.
4. **Keep total memory usage well under 150 KB.** Be mindful of large arrays, sprite sheets, and stack depth. Prefer stack allocation for small data, static allocation for persistent data.
5. **Screen bounds are 320x240.** All drawing coordinates must respect this. X: 0-319, Y: 0-239.
6. **Always call `gfx_End()` before returning from `main()`.** Failing to do so leaves the calculator in graphics mode.
7. **Cap frame rate using hardware timers** from `<sys/timers.h>`. Target ~30 FPS for smooth gameplay without burning CPU.

### Sprites and Graphics

- Use `gfx_Sprite()` and `gfx_TransparentSprite()` for drawing sprites
- Sprite format: raw byte arrays — first byte is width, second byte is height, then pixel data row by row
- Use `gfx_SetPalette()` to define custom color palettes
- Color index `0xFF` is white, `0x00` is black in the default palette
- Use `gfx_SetTransparentColor()` to set which palette index is treated as transparent

### Screen Clipping (IMPORTANT)

`gfx_PrintStringXY`, `gfx_HorizLine`, `gfx_VertLine`, and other graphx primitives do **NOT** clip to screen bounds. Drawing at negative coordinates or past 319/239 corrupts the framebuffer:

- **Negative x/y:** Wraps to garbage positions in VRAM — causes flickering and visual corruption
- **Text past x=319:** Characters wrap to the left side of the next pixel row — appears as ghost text on the wrong side of the screen

**Always manually bounds-check coordinates before calling any draw function.** For scrollable content, clip both edges (left/right, top/bottom) before drawing. See `vsepr-table/src/main.c:draw_wrapped_text()` for a working example of per-character clipping with scrolling.

### Save Data

- Use AppVars via `<fileioc.h>` for saving/loading game state or high scores
- Always check return values when reading/writing AppVars

### Code Style

- Use C11 standard
- Keep functions small and focused
- Use descriptive names for game-related variables and functions
- Comment non-obvious game logic
- Prefer `#define` constants over magic numbers for screen positions, speeds, sizes, etc.

## Key Libraries

| Header | Purpose |
|---|---|
| `<graphx.h>` | All 2D graphics: sprites, shapes, text, tilemaps, double buffering |
| `<keypadc.h>` | Fast keyboard scanning (preferred for games) |
| `<ti/getcsc.h>` | Simple blocking key input (menus, not real-time) |
| `<fileioc.h>` | Save/load data to AppVars (save files) |
| `<sys/timers.h>` | Hardware timers for frame rate control |
| `<compression.h>` | Decompress sprites/data at runtime |

## Testing

- **Emulator:** CEmu (https://ce-programming.github.io/CEmu/) — drag and drop the `.8xp` file
- **Physical calculator:** Transfer via TI Connect CE over USB
- After every code change, run `make` to verify it compiles without errors or warnings

## Existing Projects

### vsepr-table
- **Output:** `MYGAME.8xp` — VSEPR molecular geometry reference table
- **Data pipeline:** `vsepr_data.csv` → `convert.py` → `src/vsepr_data.h` (run `python convert.py vsepr_data.csv src/vsepr_data.h` to regenerate)
- **Features:** Scrollable table (arrow keys), text wrapping, screen-clipped rendering
- **Controls:** Arrow keys to scroll, [clear] to exit

## Adding a New Program

1. Create a new subdirectory at the workspace root (e.g., `my-program/`)
2. Inside it, create `src/main.c` and a `Makefile` using this template:
   ```makefile
   NAME        = PROGNAME
   ICON        = icon.png
   DESCRIPTION = "Description here"
   COMPRESSED  = NO
   CFLAGS   = -Wall -Wextra -Oz
   CXXFLAGS = -Wall -Wextra -Oz
   ifndef CEDEV
   $(error CEDEV environment variable is not set.)
   endif
   include $(CEDEV)/meta/makefile.mk
   ```
3. Add the directory name to `.vscode/tasks.json` in the `inputs[0].options` array so VSCode build tasks can find it
4. Run `make -C my-program/` to build

## Project Structure

```
TI-Game/
├── CLAUDE.md                # This file — project rules for Claude
├── .vscode/                 # Shared VSCode config
│   ├── c_cpp_properties.json
│   ├── tasks.json
│   └── settings.json
├── vsepr-table/             # VSEPR geometry reference table
│   ├── src/
│   │   ├── main.c
│   │   └── vsepr_data.h
│   ├── Makefile
│   ├── vsepr_data.csv
│   ├── convert.py
│   └── bin/                 # Build output (generated)
│       └── MYGAME.8xp
└── <new-project>/           # Future programs go here
    ├── src/
    │   └── main.c
    ├── Makefile
    └── bin/
```
