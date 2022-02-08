CC  = gcc
COPT  = -m64 -pedantic -Wall -g
# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------
LIBS =-L/usr/local/lib -lmkl_rt -lm
#GUROBILIBDIR   = $(GUROBIDIR)/$(SYSTEM)/lib
#GUROBIINCDIR   = $(GUROBIDIR)/$(SYSTEM)/include

CLNFLAGS  =  -m64 -lm -lpthread 
CFLAGS    := $(COPT)


APPL_NAME := krpca

EXEC = exec/

SOURCE = src/

OBJECT = obj/

APPL_SRCS := \
	main.c \
	read_data.c \
 	krpca.c \
 	pfkpca.c

APPL_OBJS = $(APPL_SRCS:.c=.o)

APPL_NAME := $(addprefix $(EXEC),$(APPL_NAME))

APPL_OBJS := $(addprefix $(OBJECT),$(APPL_OBJS))

APPL_SRCS := $(addprefix $(SOURCE),$(APPL_SRCS))

first: $(APPL_NAME)

$(APPL_NAME): $(APPL_OBJS)
	$(CC) $(CFLAGS) $(APPL_OBJS) $(LIBS) -o $(APPL_NAME) $(CLNFLAGS)

$(APPL_OBJS): $(OBJECT)%.o : $(SOURCE)%.c
	$(CC) -c $(CFLAGS) $< -o $@



depend :
	makedepend -I/usr/local/include/linux -I/usr/lib/gcc/x86_64-linux-gnu/4.4/include $(APPL_SRCS)

clean :
	/bin/rm -f $(APPL_OBJS)
	/bin/rm -f $(APPL_NAME)

src/main.o: src/type.h
src/read_data.o: src/type.h
