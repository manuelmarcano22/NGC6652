#
# macro definitions for the SUN/Solaris ANSI-C Compiler
#
# E.BERTIN 1997
#
SEXMACHINE=linuxpc
CC    = cc				# the C compiler
COPTS = -pg -O2 -malign-double -finline-functions -funroll-loops -DPC_LINUX	# options for the C compiler
CD    = cd
CP    = cp
RM    = rm -f
