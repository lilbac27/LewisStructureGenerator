#ifndef UI_TEXT_H
#define UI_TEXT_H

#include <stddef.h>
#include <stdint.h>

char *int_to_str(int val, char *buf);
void safe_print(const char *s, int x, int y);
void append_str(char *dst, size_t cap, const char *src);
void append_int(char *dst, size_t cap, int val);
uint8_t text_color_for_bg(uint8_t bg);

#endif
