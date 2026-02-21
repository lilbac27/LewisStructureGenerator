# ----------------------------
# Makefile Options
# ----------------------------

NAME        = MYGAME
ICON        = icon.png
DESCRIPTION = "My TI-84 Plus CE Game"
COMPRESSED  = NO

CFLAGS   = -Wall -Wextra -Oz
CXXFLAGS = -Wall -Wextra -Oz

# ----------------------------

ifndef CEDEV
$(error CEDEV environment variable is not set. Set it to your CE toolchain install path.)
endif

include $(CEDEV)/meta/makefile.mk
