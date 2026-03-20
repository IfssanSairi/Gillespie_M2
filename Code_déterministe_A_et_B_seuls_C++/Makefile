# Compilateur
CXX = g++

# Options de compilation
CXXFLAGS = -Wall -Wextra -std=c++14

# include paths
INCLUDE = "-I/opt/homebrew/Cellar/boost/1.90.0_1/include"

# Fichiers sources
SRCS = main.cpp Entite.cpp

# Fichiers objets
OBJS = $(SRCS:.cpp=.o)

# Exécutable final
TARGET = deterministic

# Règle par défaut
all: $(TARGET)
	rm -f $(OBJS)

# Règle pour créer l'exécutable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^

# Règle pour compiler les fichiers .cpp en .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

# Règle pour nettoyer les fichiers objets et l'exécutable
clean:
	rm -f $(OBJS) $(TARGET)

# Règle pour tout nettoyer et recompiler
rebuild: clean all
