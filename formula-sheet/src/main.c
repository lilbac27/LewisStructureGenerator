#include <graphx.h>
#include <keypadc.h>
#include <sys/timers.h>
#include <string.h>

/* Layout constants */
#define FONT_WIDTH      8
#define FONT_HEIGHT     8
#define ROW_HEIGHT      20
#define TITLE_HEIGHT    22
#define SWATCH_SIZE     12
#define SWATCH_PAD      4
#define TEXT_X          (SWATCH_SIZE + SWATCH_PAD + 4)
#define WAVE_COL_X      200

/* EM spectrum entry */
typedef struct {
    const char *name;
    const char *wavelength;
    uint8_t color;       /* palette index for swatch */
    uint8_t is_visible;  /* 1 = visible light band */
} SpectrumEntry;

/* Custom palette colors for visible spectrum (indices 0xE0-0xE6) */
#define PAL_RED     0xE0
#define PAL_ORANGE  0xE1
#define PAL_YELLOW  0xE2
#define PAL_GREEN   0xE3
#define PAL_BLUE    0xE4
#define PAL_INDIGO  0xE5
#define PAL_VIOLET  0xE6

/* Spectrum data — longest wavelength to shortest */
static const SpectrumEntry spectrum[] = {
    {"Radio waves",     "> 1 m",            0x00, 0},
    {"Microwaves",      "1 mm - 1 m",       0x00, 0},
    {"Infrared",        "700 nm - 1 mm",    0x00, 0},
    {"  Red",           "620 - 700 nm",     PAL_RED,    1},
    {"  Orange",        "590 - 620 nm",     PAL_ORANGE, 1},
    {"  Yellow",        "570 - 590 nm",     PAL_YELLOW, 1},
    {"  Green",         "495 - 570 nm",     PAL_GREEN,  1},
    {"  Blue",          "450 - 495 nm",     PAL_BLUE,   1},
    {"  Indigo",        "420 - 450 nm",     PAL_INDIGO, 1},
    {"  Violet",        "380 - 420 nm",     PAL_VIOLET, 1},
    {"Ultraviolet",     "10 - 380 nm",      0x00, 0},
    {"X-rays",          "0.01 - 10 nm",     0x00, 0},
    {"Gamma rays",      "< 0.01 nm",        0x00, 0},
};

#define NUM_ENTRIES (sizeof(spectrum) / sizeof(spectrum[0]))

/* Set custom palette entries for visible light colors */
static void setup_palette(void)
{
    /* gfx_palette is a uint16_t array; format is 1555 color (rrrrrggggggbbbbb in 565, but CE uses 1555: gggbbbbb0rrrrrggg) */
    /* Actually CE uses RGB1555 stored as 16-bit: bit layout in memory is gggbbbbb_xrrrrrgg */
    /* Simpler: use gfx_RGBTo1555(r,g,b) macro */
    gfx_palette[PAL_RED]    = gfx_RGBTo1555(220, 40, 40);
    gfx_palette[PAL_ORANGE] = gfx_RGBTo1555(240, 140, 20);
    gfx_palette[PAL_YELLOW] = gfx_RGBTo1555(240, 230, 30);
    gfx_palette[PAL_GREEN]  = gfx_RGBTo1555(40, 180, 50);
    gfx_palette[PAL_BLUE]   = gfx_RGBTo1555(40, 80, 220);
    gfx_palette[PAL_INDIGO] = gfx_RGBTo1555(80, 40, 160);
    gfx_palette[PAL_VIOLET] = gfx_RGBTo1555(140, 30, 180);
}

/* Draw a single clipped string at (x, y). Skips if off-screen. */
static void draw_clipped_string(const char *str, int x, int y)
{
    if (y < 0 || y + FONT_HEIGHT > GFX_LCD_HEIGHT) return;
    if (x >= GFX_LCD_WIDTH) return;

    int len = (int)strlen(str);
    int skip = 0;

    /* Clip left */
    if (x < 0) {
        skip = (-x + FONT_WIDTH - 1) / FONT_WIDTH;
        x += skip * FONT_WIDTH;
    }

    int visible = len - skip;
    int max_fit = (GFX_LCD_WIDTH - x) / FONT_WIDTH;
    if (visible > max_fit) visible = max_fit;
    if (visible <= 0) return;

    char buf[48];
    if (visible >= (int)sizeof(buf)) visible = (int)sizeof(buf) - 1;
    memcpy(buf, str + skip, visible);
    buf[visible] = '\0';
    gfx_PrintStringXY(buf, x, y);
}

int main(void)
{
    gfx_Begin();
    gfx_SetDrawBuffer();
    setup_palette();

    int scroll_y = 0;
    int total_height = TITLE_HEIGHT + (int)(NUM_ENTRIES * ROW_HEIGHT) + ROW_HEIGHT;
    int max_scroll_y = total_height - GFX_LCD_HEIGHT;
    if (max_scroll_y < 0) max_scroll_y = 0;

    bool running = true;
    while (running) {
        /* Input */
        kb_Scan();
        if (kb_Data[6] & kb_Clear)
            running = false;
        if (kb_Data[7] & kb_Up) {
            scroll_y -= 5;
            if (scroll_y < 0) scroll_y = 0;
        }
        if (kb_Data[7] & kb_Down) {
            scroll_y += 5;
            if (scroll_y > max_scroll_y) scroll_y = max_scroll_y;
        }

        /* Draw */
        gfx_FillScreen(0xFF); /* white */

        /* Title */
        {
            int ty = 4 - scroll_y;
            if (ty >= 0 && ty + FONT_HEIGHT <= GFX_LCD_HEIGHT) {
                gfx_SetTextFGColor(0x00);
                gfx_SetTextScale(1, 1);
                gfx_PrintStringXY("EM Spectrum - Wavelength Chart", 20, ty);
            }
            /* Title underline */
            int uy = TITLE_HEIGHT - 2 - scroll_y;
            if (uy >= 0 && uy < GFX_LCD_HEIGHT) {
                gfx_SetColor(0x00);
                gfx_HorizLine(4, uy, 312);
            }
        }

        /* Column headers */
        {
            int hy = TITLE_HEIGHT + 2 - scroll_y;
            if (hy >= 0 && hy + FONT_HEIGHT <= GFX_LCD_HEIGHT) {
                gfx_SetTextFGColor(0x00);
                draw_clipped_string("Type", TEXT_X, hy);
                draw_clipped_string("Wavelength", WAVE_COL_X, hy);
            }
            int div_y = TITLE_HEIGHT + 2 + FONT_HEIGHT + 2 - scroll_y;
            if (div_y >= 0 && div_y < GFX_LCD_HEIGHT) {
                gfx_SetColor(0xB5);
                gfx_HorizLine(4, div_y, 312);
            }
        }

        int header_offset = TITLE_HEIGHT + ROW_HEIGHT;

        /* Spectrum rows */
        for (int i = 0; i < (int)NUM_ENTRIES; i++) {
            int y = header_offset + i * ROW_HEIGHT - scroll_y;

            /* Cull rows off screen */
            if (y + ROW_HEIGHT < 0) continue;
            if (y >= GFX_LCD_HEIGHT) break;

            int text_y = y + (ROW_HEIGHT - FONT_HEIGHT) / 2;

            /* Color swatch for visible light bands */
            if (spectrum[i].is_visible) {
                int sx = 4;
                int sy = y + (ROW_HEIGHT - SWATCH_SIZE) / 2;
                if (sy >= 0 && sy + SWATCH_SIZE < GFX_LCD_HEIGHT &&
                    sx >= 0 && sx + SWATCH_SIZE < GFX_LCD_WIDTH) {
                    gfx_SetColor(spectrum[i].color);
                    gfx_FillRectangle(sx, sy, SWATCH_SIZE, SWATCH_SIZE);
                    /* Border */
                    gfx_SetColor(0x00);
                    gfx_Rectangle(sx, sy, SWATCH_SIZE, SWATCH_SIZE);
                }
            } else {
                /* Gray dash for non-visible bands */
                int dx = 4 + 2;
                int dy = y + ROW_HEIGHT / 2;
                if (dy >= 0 && dy < GFX_LCD_HEIGHT && dx >= 0 && dx + 8 < GFX_LCD_WIDTH) {
                    gfx_SetColor(0xB5);
                    gfx_HorizLine(dx, dy, 8);
                }
            }

            /* Name */
            gfx_SetTextFGColor(0x00);
            if (text_y >= 0 && text_y + FONT_HEIGHT <= GFX_LCD_HEIGHT) {
                draw_clipped_string(spectrum[i].name, TEXT_X, text_y);
            }

            /* Wavelength */
            if (text_y >= 0 && text_y + FONT_HEIGHT <= GFX_LCD_HEIGHT) {
                draw_clipped_string(spectrum[i].wavelength, WAVE_COL_X, text_y);
            }

            /* Row separator */
            int div_y = y + ROW_HEIGHT - 1;
            if (div_y >= 0 && div_y < GFX_LCD_HEIGHT) {
                gfx_SetColor(0xD0);
                gfx_HorizLine(4, div_y, 312);
            }
        }

        /* Footer hint */
        {
            int fy = header_offset + (int)NUM_ENTRIES * ROW_HEIGHT + 4 - scroll_y;
            if (fy >= 0 && fy + FONT_HEIGHT <= GFX_LCD_HEIGHT) {
                gfx_SetTextFGColor(0xB5);
                draw_clipped_string("[clear] to exit", 100, fy);
            }
        }

        gfx_SwapDraw();
    }

    gfx_End();
    return 0;
}
