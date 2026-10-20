CC = gcc
CFLAGS = -Wall -Wextra -Wpedantic -Werror -O2 -std=c99

SRCS = tropical.c karp.c test_week1.c
OBJS = $(SRCS:.c=.o)
TARGET = test_week1

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.c tropical.h
	$(CC) $(CFLAGS) -c $< -o $@

test: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET) $(TARGET).exe

.PHONY: all test clean
