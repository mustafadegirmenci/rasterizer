CC = g++
CFLAGS = -std=c++17 -O3

SRC_DIR = src
INCLUDE_DIR = include

INCLUDE_HEADERS = $(wildcard $(INCLUDE_DIR)/*.h)
INCLUDE_SOURCES = $(wildcard $(SRC_DIR)/*.cpp)

SUBDIR_HEADERS = $(shell find $(INCLUDE_DIR) -type f -name "*.h")
SUBDIR_SOURCES = $(shell find $(SRC_DIR) -type f -name "*.cpp")

HEADERS = $(INCLUDE_HEADERS) $(SUBDIR_HEADERS)
SOURCES = $(INCLUDE_SOURCES) $(SUBDIR_SOURCES)

rasterizer: $(HEADERS) $(SOURCES)
	$(CC) $^ -o$@ $(CFLAGS)

clean:
	rm -f rasterizer
