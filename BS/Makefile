# Makefile
-include $(DEPS)

CPPFLAGS := `pkg-config opencv --cflags` `pkg-config opencv --libs`

CC := g++

PROG := main
SRCS := main.cpp sm.cpp
OBJS := $(SRCS:.cpp=.o)
DEPS := $(SRCS:.cpp=.d)

all: $(PROG)

$(PROG): $(OBJS) 
	$(CC) -O2 $^ $(CPPFLAGS) -o $@

%.o: %.cpp
	$(CC) -O2 $< -c -MMD -MP -std=c++11

clean:
	rm -rf $(PROG) $(OBJS) $(DEPS)
