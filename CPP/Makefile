# compilateur c++
CC = g++

# options de compilation
CFLAGS = -Wall -Wextra -Werror -std=c++11 -g

# options de linkage

# nom de l'executable
EXEC = Simulation

# fichiers sources
SRC = $(wildcard *.cpp)

# fichiers objets
OBJ = $(SRC:.cpp=.o)

# regle de compilation
all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) -o $@ $^

%.o: %.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

# regle de nettoyage
clean:
	rm -rf *.o

# regle de nettoyage totale
mrproper: clean
	rm -rf $(EXEC)

# regle de compilation automatique
.PHONY: clean mrproper

# fin du fichier Makefile
