#include "ui_text.h"

#include <graphx.h>
#include <stdbool.h>
#include <string.h>

#include "lewis_model.h"
#include "ui_theme.h"

void safe_print(const char *s, int x, int y)
{
    if (y < 0 || y >= SCR_H - 8) return;
    if (x < 0 || x >= SCR_W) return;

    int max_chars = (SCR_W - x) / 8;
    if (max_chars <= 0) return;

    int len = (int)strlen(s);
    if (len <= max_chars) {
        gfx_PrintStringXY(s, x, y);
    } else {
        char buf[42];
        if (max_chars > 41) max_chars = 41;
        memcpy(buf, s, (size_t)max_chars);
        buf[max_chars] = '\0';
        gfx_PrintStringXY(buf, x, y);
    }
}

uint8_t safe_print_wrapped(const char *text, int x, int y, int max_width, uint8_t max_lines)
{
    if (text == NULL || max_lines == 0 || max_width < 8) {
        return 0;
    }

    int max_chars = max_width / 8;
    if (max_chars <= 0) {
        return 0;
    }

    uint8_t lines_drawn = 0;
    const char *p = text;

    while (*p != '\0' && lines_drawn < max_lines) {
        while (*p == ' ') p++;
        if (*p == '\0') break;

        int remaining = (int)strlen(p);
        int chunk = remaining;

        if (chunk > max_chars) {
            chunk = max_chars;
            while (chunk > 0 && p[chunk] != ' ') {
                chunk--;
            }
            if (chunk == 0) {
                chunk = max_chars;
            }
        }

        while (chunk > 0 && p[chunk - 1] == ' ') {
            chunk--;
        }
        if (chunk <= 0) {
            p++;
            continue;
        }

        int draw_y = y + (int)lines_drawn * 8;
        if (draw_y >= 0 && draw_y <= SCR_H - 8) {
            int draw_x = x;
            int skip_chars = 0;

            if (draw_x < 0) {
                skip_chars = (-draw_x + 7) / 8;
                draw_x += skip_chars * 8;
            }

            int visible_chars = chunk - skip_chars;
            int fit_chars = (SCR_W - draw_x) / 8;
            if (visible_chars > fit_chars) {
                visible_chars = fit_chars;
            }

            if (visible_chars > 0 && draw_x < SCR_W) {
                char buf[48];
                if (visible_chars > (int)sizeof(buf) - 1) {
                    visible_chars = (int)sizeof(buf) - 1;
                }

                memcpy(buf, p + skip_chars, (size_t)visible_chars);
                buf[visible_chars] = '\0';
                gfx_PrintStringXY(buf, draw_x, draw_y);
            }
        }

        lines_drawn++;
        p += chunk;
        while (*p == ' ') p++;
    }

    return lines_drawn;
}

char *int_to_str(int val, char *buf)
{
    if (val == 0) {
        buf[0] = '0';
        buf[1] = '\0';
        return buf;
    }

    char tmp[12];
    int i = 0;
    bool neg = false;
    unsigned int mag;

    if (val < 0) {
        neg = true;
        mag = (unsigned int)(-(val + 1)) + 1u;
    } else {
        mag = (unsigned int)val;
    }

    while (mag > 0) {
        tmp[i++] = (char)('0' + (mag % 10u));
        mag /= 10u;
    }

    int j = 0;
    if (neg) {
        buf[j++] = '-';
    }
    while (i > 0) {
        buf[j++] = tmp[--i];
    }
    buf[j] = '\0';
    return buf;
}

void append_str(char *dst, size_t cap, const char *src)
{
    if (cap == 0) return;
    size_t len = strlen(dst);
    if (len + 1 >= cap) return;
    strncat(dst, src, cap - len - 1);
}

void append_int(char *dst, size_t cap, int val)
{
    char nb[12];
    int_to_str(val, nb);
    append_str(dst, cap, nb);
}

uint8_t text_color_for_bg(uint8_t bg)
{
    return (bg == UI_SELECTED_BG) ? UI_SELECTED_TEXT : UI_TEXT;
}
