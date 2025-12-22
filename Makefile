CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall -Wextra -Iinclude -fopenmp
LDFLAGS = -lgsl -lgslcblas -lm -fopenmp

NAME = brownian_dynamics
TARGET = bin/$(NAME)

SRC_DIR = src
OBJ_DIR = build
BIN_DIR = bin

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

all: directories $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

directories:
	mkdir -p $(OBJ_DIR) $(BIN_DIR)

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

run: all
	./$(TARGET) 2D

.PHONY: all clean directories run
