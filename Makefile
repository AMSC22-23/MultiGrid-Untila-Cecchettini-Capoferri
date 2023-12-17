CXX := g++
CXXFLAGS := -Wall -std=c++20 -O3

# There are all var for directories declaration
SRCDIR := src
INCDIR := include
BINDIR := bin
TARGET := Multigrid  # the name of executable

SRC := $(wildcard $(SRCDIR)/*.cpp)
OBJ := $(patsubst $(SRCDIR)/%.cpp, $(BINDIR)/%.o, $(SRC))

.PHONY: all clean

all: $(TARGET)

# Rules for linking
$(TARGET): $(OBJ)
	@mkdir -p $(BINDIR)
	@$(CXX) $(CXXFLAGS) $^ -o $@

# Rules for compiling
$(BINDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(BINDIR)
	@$(CXX) $(CXXFLAGS) -I $(INCDIR) -c $< -o $@

clean:
	@rm -rf $(BINDIR)
	@rm $(TARGET)