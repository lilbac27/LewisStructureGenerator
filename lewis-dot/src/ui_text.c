#include "ui_text.h"

#include <graphx.h>
#include <stdbool.h>
#include <string.h>

#include "lewis_model.h"

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
    switch (bg) {
        case COL_BLACK:
        case COL_DKGRAY:
        case COL_BLUE:
        case COL_RED:
            return COL_WHITE;
        default:
            return COL_BLACK;
    }
}
