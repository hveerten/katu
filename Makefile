CC = clang
CFLAGS = -Wall -Wextra -Igsl -std=gnu99
OPTIMIZATION_FLAGS = -O2 -march=native -flto 
OPTIMIZATION_FLAGS_EXTRA = -DNDEBUG
DEBUG_FLAGS = -ggdb
DEBUG_FLAGS_EXTRA = -fsanitize=undefined
LDFLAGS = -lm -lgsl -lcblas -lpthread


CFLAGS+=$(OPTIMIZATION_FLAGS)
#CFLAGS+=$(OPTIMIZATION_FLAGS_EXTRA)
CFLAGS+=$(DEBUG_FLAGS)
#CFLAGS+=$(DEBUG_FLAGS_EXTRA)
CFLAGS+=-DUSE_THREADS -DUSE_THREAD_POOL -fopenmp
#CFLAGS+=-DUSE_THREADS -DUSE_MY_THREADS -fopenmp
CFLAGS+=-DELECTRON_STEADY_STATE=1
CFLAGS+=-DPROTON_STEADY_STATE=1

EXTERNAL_SOURCES = \
	./src/external_libs/toml.c
EXTERNAL_HEADERS = \
	./src/external_libs/toml.h
EXTERNAL_OBJECTS = \
	./src/external_libs/toml.o

UTILS_SOURCES = \
	./src/utils/tpool.c
UTILS_HEADERS = \
	./src/utils/tpool.h
UTILS_OBJECTS = \
	./src/utils/tpool.o

SOURCES = \
	./src/acceleration.c 		\
	./src/bethe_heitler.c 		\
	./src/config.c		 		\
	./src/distribution.c 		\
	./src/escape.c				\
	./src/decay_and_escape.c	\
	./src/inverse_compton.c 	\
	./src/muon_decay.c 			\
	./src/pair_production.c 	\
	./src/pion_decay.c 			\
	./src/pion_production.c 	\
	./src/population.c			\
	./src/state.c				\
	./src/state_step.c			\
	./src/state_step_common.c	\
	./src/synchrotron.c			\
	./src/test.c

OBJECTS = \
	./src/acceleration.o 		\
	./src/bethe_heitler.o 		\
	./src/config.o		 		\
	./src/distribution.o 		\
	./src/escape.o				\
	./src/decay_and_escape.o	\
	./src/inverse_compton.o 	\
	./src/muon_decay.o 			\
	./src/pair_production.o 	\
	./src/pion_decay.o 			\
	./src/pion_production.o 	\
	./src/population.o 			\
	./src/state.o				\
	./src/state_step.o			\
	./src/state_step_common.o	\
	./src/synchrotron.o			\
	./src/test.o

HEADERS = \
	./src/acceleration.h 		\
	./src/bethe_heitler.h 		\
	./src/config.h		 		\
	./src/constants.h			\
	./src/distribution.h		\
	./src/escape.h				\
	./src/decay_and_escape.h	\
	./src/inverse_compton.h		\
	./src/muon_decay.h 			\
	./src/pair_production.h		\
	./src/pion_decay.h 			\
	./src/pion_production.h 	\
	./src/population.h			\
	./src/state.h				\
	./src/state_step.h			\
	./src/state_step_common.h	\
	./src/synchrotron.h			\

TARGET = model

all: $(TARGET)

$(TARGET): $(OBJECTS) $(HEADERS) $(SOURCES) $(EXTERNAL_OBJECTS) $(UTILS_OBJECTS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJECTS) $(EXTERNAL_OBJECTS) $(UTILS_OBJECTS) $(LDFLAGS)
