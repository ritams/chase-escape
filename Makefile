CC = gcc
CFLAGS = -O3 -Wall -std=c99
LDFLAGS = -lm

# Target executables
TARGET = chase_escape_parallel
TARGET_LAMBDA = chase_escape_lambda_scan

# Source files
SRC = src/chase_escape_parallel.c
SRC_LAMBDA = src/chase_escape_lambda_scan.c

# Default target
all: $(TARGET) $(TARGET_LAMBDA)

# Compile the original program
$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

# Compile the lambda scan program
$(TARGET_LAMBDA): $(SRC_LAMBDA)
	$(CC) $(CFLAGS) -o $(TARGET_LAMBDA) $(SRC_LAMBDA) $(LDFLAGS)

# Clean up
clean:
	rm -f $(TARGET) $(TARGET_LAMBDA)

# Make executable
.PHONY: all clean 