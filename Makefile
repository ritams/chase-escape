CC = gcc
CFLAGS = -O3 -Wall -std=c99
LDFLAGS = -lm

# Target executable
TARGET = chase_escape_parallel

# Source file
SRC = src/chase_escape_parallel.c

# Default target
all: $(TARGET)

# Compile the program
$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

# Clean up
clean:
	rm -f $(TARGET)

# Make executable
.PHONY: all clean 