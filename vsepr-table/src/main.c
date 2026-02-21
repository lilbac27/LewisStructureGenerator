#include <graphx.h>
#include <keypadc.h>
#include <sys/timers.h>
#include <ti/getcsc.h>
#include <string.h>

#include "vsepr_data.h"

#define ROW_HEIGHT 24
#define COL_WIDTH 130
#define NUM_COLS 6
#define HEADER_HEIGHT 26
#define CELL_PADDING 2
#define FONT_WIDTH 8
#define FONT_HEIGHT 8
#define MAX_CHARS_PER_LINE ((COL_WIDTH - 2 * CELL_PADDING) / FONT_WIDTH)

/* Draw text wrapped within a cell, with screen bounds clipping. */
static void draw_wrapped_text(const char *text, int x, int y, int max_width)
{
    int max_chars = max_width / FONT_WIDTH;
    int len = (int)strlen(text);
    int line = 0;
    char buf[32];

    while (len > 0) {
        int draw_y = y + line * FONT_HEIGHT;

        /* Stop if below screen */
        if (draw_y >= GFX_LCD_HEIGHT) break;

        /* Calculate how many chars fit on this line */
        int chunk = len;
        if (chunk > max_chars) {
            chunk = max_chars;
            while (chunk > 0 && text[chunk] != ' ')
                chunk--;
            if (chunk == 0)
                chunk = max_chars;
        }

        /* Only draw if this line is vertically on screen */
        if (draw_y >= 0) {
            int draw_x = x;
            int skip = 0;

            /* Clip characters that are off the left edge */
            if (draw_x < 0) {
                skip = (-draw_x + FONT_WIDTH - 1) / FONT_WIDTH;
                draw_x = x + skip * FONT_WIDTH;
            }

            int visible = chunk - skip;

            /* Clip characters that would extend past the right edge */
            int max_fit = (GFX_LCD_WIDTH - draw_x) / FONT_WIDTH;
            if (visible > max_fit)
                visible = max_fit;

            if (visible > 0 && draw_x < GFX_LCD_WIDTH) {
                if (visible >= (int)sizeof(buf))
                    visible = (int)sizeof(buf) - 1;
                memcpy(buf, text + skip, visible);
                buf[visible] = '\0';
                gfx_PrintStringXY(buf, draw_x, draw_y);
            }
        }

        text += chunk;
        len -= chunk;
        if (len > 0 && *text == ' ') {
            text++;
            len--;
        }
        line++;
    }
}

int main(void)
{
    /* Initialize graphics */
    gfx_Begin();
    gfx_SetDrawBuffer();

    int scroll_x = 0;
    int scroll_y = 0;

    int max_scroll_x = (NUM_COLS * COL_WIDTH) - GFX_LCD_WIDTH;
    if (max_scroll_x < 0) max_scroll_x = 0;
    
    int total_height = HEADER_HEIGHT + (VSEPR_NUM_ROWS * ROW_HEIGHT);
    int max_scroll_y = total_height - GFX_LCD_HEIGHT;
    if (max_scroll_y < 0) max_scroll_y = 0;

    /* Main game loop */
    bool running = true;
    while (running)
    {
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
        if (kb_Data[7] & kb_Left) {
            scroll_x -= 5;
            if (scroll_x < 0) scroll_x = 0;
        }
        if (kb_Data[7] & kb_Right) {
            scroll_x += 5;
            if (scroll_x > max_scroll_x) scroll_x = max_scroll_x;
        }

        /* Draw */
        gfx_FillScreen(0xFF); // white background
        
        // Draw grid and data
        gfx_SetTextFGColor(0x00);

        // Headers
        for (int c = 0; c < NUM_COLS; c++) {
            int x = (c * COL_WIDTH) - scroll_x;
            if (x > GFX_LCD_WIDTH || x + COL_WIDTH < 0) continue; // Culling
            
            draw_wrapped_text(vsepr_headers[c], x + CELL_PADDING, 2 - scroll_y, COL_WIDTH - 2 * CELL_PADDING);
        }

        // Header divider line
        {
            int hy = HEADER_HEIGHT - scroll_y - 2;
            if (hy >= 0 && hy < GFX_LCD_HEIGHT) {
                gfx_SetColor(0x00);
                gfx_HorizLine(0, hy, GFX_LCD_WIDTH);
            }
        }

        // Data rows
        for (int r = 0; r < VSEPR_NUM_ROWS; r++) {
            int y = HEADER_HEIGHT + (r * ROW_HEIGHT) - scroll_y;
            if (y >= GFX_LCD_HEIGHT) break;      // Below screen, done
            if (y + ROW_HEIGHT < 0) continue;    // Above screen, skip

            const char *row_data[NUM_COLS] = {
                vsepr_data[r].valence_pairs,
                vsepr_data[r].ep_geometry,
                vsepr_data[r].bond_pairs,
                vsepr_data[r].lone_pairs,
                vsepr_data[r].shape,
                vsepr_data[r].hybridization
            };

            for (int c = 0; c < NUM_COLS; c++) {
                int x = (c * COL_WIDTH) - scroll_x;
                if (x >= GFX_LCD_WIDTH) break;       // Right of screen, done
                if (x + COL_WIDTH < 0) continue;     // Left of screen, skip

                draw_wrapped_text(row_data[c], x + CELL_PADDING, y + CELL_PADDING, COL_WIDTH - 2 * CELL_PADDING);
            }

            // Row divider — only draw if on screen
            int div_y = y + ROW_HEIGHT - 1;
            if (div_y >= 0 && div_y < GFX_LCD_HEIGHT) {
                gfx_SetColor(0xD0);
                gfx_HorizLine(0, div_y, GFX_LCD_WIDTH);
            }
        }

        // Col dividers — clamp vertical extent to screen
        gfx_SetColor(0x00);
        for (int c = 0; c <= NUM_COLS; c++) {
            int x = (c * COL_WIDTH) - scroll_x;
            if (x >= 0 && x < GFX_LCD_WIDTH) {
                gfx_VertLine(x, 0, GFX_LCD_HEIGHT);
            }
        }

        gfx_SwapDraw();
    }

    gfx_End();
    return 0;
}
