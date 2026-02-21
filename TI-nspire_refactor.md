# TI-Nspire Refactor Guide

## Goal
Port all existing TI-84 Plus CE programs in this repo to TI-Nspire while preserving current behavior and chemistry logic.

This document is a handoff for a future agent. It describes what to change and in what order. It does not perform the refactor.

## Programs in Scope
- `lewis-dot` (`lewis-dot/src/main.c`)
- `formula-sheet` (`formula-sheet/src/main.c`)
- `vsepr-table` (`vsepr-table/src/main.c`, `vsepr-table/src/vsepr_data.h`)

## Current CE-Specific Dependencies
All three apps are tightly coupled to CE libraries and build tooling.

- Graphics: `<graphx.h>` and `gfx_*` calls
- Input: `<keypadc.h>` and `kb_Scan()` / `kb_Data[]`
- Timing: `<sys/timers.h>` and CE timer registers
- Build: CE `Makefile` template via `include $(CEDEV)/meta/makefile.mk`
- Output format: `.8xp`

Reference locations:
- `lewis-dot/src/main.c:19`, `lewis-dot/src/main.c:20`, `lewis-dot/src/main.c:21`
- `lewis-dot/src/main.c:1485` (direct timer register usage)
- `formula-sheet/src/main.c:1`, `formula-sheet/src/main.c:2`, `formula-sheet/src/main.c:3`
- `vsepr-table/src/main.c:1`, `vsepr-table/src/main.c:2`, `vsepr-table/src/main.c:3`
- `lewis-dot/Makefile`, `formula-sheet/Makefile`, `vsepr-table/Makefile`

## Required Decisions Before Coding
A future agent should confirm these first, because they affect every implementation detail:

1. Target device family
- Prefer CX/CX II only (color, 320x240). If non-CX support is required, graphics work increases.

2. Runtime/toolchain
- Choose one stack and stick to it across all programs (for example, Ndless + nSDL).
- Define output artifact type (`.tns` packaging path).

3. Font/text rendering backend
- Current code assumes 8x8 fixed-width text in multiple places.
- Pick a backend that supports fixed metrics or provide a wrapper that emulates fixed-width behavior.

## Recommended Port Strategy
Do not rewrite chemistry/game logic first. Isolate platform calls first.

1. Introduce a platform abstraction layer (new shared module)
- Add a small API for graphics, input, and frame timing.
- Keep app logic in each existing `main.c` mostly intact.

2. Convert CE code to use the abstraction (no behavior changes)
- Replace direct `gfx_*`, `kb_*`, and timer register calls with wrapper functions.
- Verify CE build still works after this step.

3. Implement TI-Nspire backend for that abstraction
- Hook wrappers to Nspire graphics/input/timing APIs.
- Only then adjust any rendering or key mapping quirks.

This reduces risk and allows side-by-side CE and Nspire builds.

## Suggested Platform API (minimal)
A future agent can adjust names, but the capabilities should exist:

```c
/* lifecycle */
void plat_init(void);
void plat_shutdown(void);

/* frame */
void plat_begin_frame(void);
void plat_present_frame(void);
void plat_clear(uint16_t color);
uint32_t plat_ticks_ms(void);
void plat_sleep_until(uint32_t target_ms);

/* drawing */
void plat_set_draw_color(uint16_t color);
void plat_set_text_color(uint16_t fg, uint16_t bg);
void plat_set_text_scale(uint8_t sx, uint8_t sy);
void plat_draw_line(int x1, int y1, int x2, int y2);
void plat_draw_rect(int x, int y, int w, int h);
void plat_fill_rect(int x, int y, int w, int h);
void plat_fill_circle(int x, int y, int r);
void plat_draw_text_xy(const char *s, int x, int y);

/* input */
typedef struct {
    bool up, down, left, right;
    bool enter, clear, del;
    bool alpha, second, mode;
} plat_input_t;

void plat_scan_input(plat_input_t *in);
```

## API Migration Map
Direct mapping work that must happen:

- `gfx_Begin()` / `gfx_End()` -> `plat_init()` / `plat_shutdown()`
- `gfx_SetDrawBuffer()` + `gfx_SwapDraw()` -> `plat_begin_frame()` + `plat_present_frame()`
- `gfx_FillScreen`, `gfx_Line`, `gfx_Rectangle`, `gfx_FillRectangle`, `gfx_FillCircle` -> draw wrappers
- `gfx_PrintStringXY`, `gfx_SetTextFGColor`, `gfx_SetTextBGColor`, `gfx_SetTextScale` -> text wrappers
- `kb_Scan()`, `kb_Data[]` checks -> `plat_scan_input()` fields
- CE timer registers (`timer_Control`, `timer_1_*`, `TIMER1_*`) -> millisecond frame pacing wrapper

## Project-Specific Notes

### `lewis-dot`
Main work:
- Replace CE timer register frame loop with backend-agnostic frame pacing.
- Preserve key repeat behavior and current controls.
- Keep existing chemical generation logic unchanged (`generate_resonance` and related helpers).

Important references:
- `lewis-dot/src/main.c:948` (`generate_resonance`)
- `lewis-dot/src/main.c:1187` (`draw_lewis`)
- `lewis-dot/src/main.c:1463` (`main` loop)

Controls currently used:
- arrows, enter, del, alpha, 2nd, mode, clear

### `formula-sheet`
Main work:
- Replace palette-index assumptions with explicit color values for Nspire backend.
- Keep scrolling and clipped text behavior.

Important references:
- `formula-sheet/src/main.c:53` (`setup_palette`)
- `formula-sheet/src/main.c:68` (`draw_clipped_string`)
- `formula-sheet/src/main.c:94` (`main` loop)

Controls currently used:
- up, down, clear

### `vsepr-table`
Main work:
- Preserve table scrolling and wrapped-cell rendering.
- Keep generated data header flow unchanged unless build system changes require updates.

Important references:
- `vsepr-table/src/main.c:19` (`draw_wrapped_text`)
- `vsepr-table/src/main.c:79` (`main` loop)
- `vsepr-table/src/vsepr_data.h` (static data table)
- `vsepr-table/convert.py` (CSV -> header generator)

Controls currently used:
- up, down, left, right, clear

## Build and Repo Changes Required
A future agent should plan these repo-level edits:

1. Add a shared platform directory
- Example: `common/platform/` with `platform.h`, CE backend, Nspire backend.

2. Split build targets
- Keep existing CE `Makefile` behavior.
- Add Nspire build files (or parallel `Makefile.nspire`) per project.
- Add VSCode tasks for Nspire build/clean in `.vscode/tasks.json`.

3. Keep generated artifacts out of source control where possible
- Current repo includes CE `bin/` and `obj/` outputs.
- Decide whether to continue this practice for Nspire outputs.

## Key Compatibility Risks
- Text width assumptions: many UI calculations assume fixed 8-pixel glyph width.
- Color model differences: CE palette indices versus Nspire color format.
- Input differences: CE key layout constants do not translate directly.
- Clipping behavior: existing code sometimes relies on manual clipping; backend behavior may differ.
- Frame pacing: CE hardware timer logic is direct-register based and must be replaced.

## Validation Checklist for Future Agent
For each program after Nspire port:

1. App starts, renders first screen, and exits cleanly.
2. All listed controls function correctly.
3. Scrolling behavior matches CE version.
4. No off-screen text corruption or drawing artifacts.
5. Frame pacing is stable (no busy-loop burn or visible stutter).
6. No chemistry/result regression in `lewis-dot` output for representative molecules.

## Suggested Execution Order
1. Confirm Nspire toolchain/runtime choice.
2. Create platform abstraction and CE backend.
3. Switch CE apps to abstraction without behavior changes.
4. Add Nspire backend.
5. Port/test `formula-sheet` (smallest surface area).
6. Port/test `vsepr-table` (text + scrolling heavy).
7. Port/test `lewis-dot` (largest feature set).
8. Update docs and tasks with exact Nspire build/run steps.

## Non-Goals for Initial Port
- Do not redesign UI/UX.
- Do not rewrite chemistry algorithms.
- Do not optimize beyond what is needed for parity and stability.

