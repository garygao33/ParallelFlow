CXX = g++
OMPFLAGS  = -fopenmp
CXXFLAGS = -std=c++17 -Wall -Wextra $(OMPFLAGS)
LDFLAGS   = $(OMPFLAGS)


SRC_DIR = src
BUILD_DIR = build

SRCS = $(SRC_DIR)/parallel.cpp

OBJS = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRCS))

TARGET = flow

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) $(LDFLAGS) -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR) $(TARGET)

.PHONY: all clean
