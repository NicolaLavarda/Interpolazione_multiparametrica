CC = c++
CFLAGS = -Iinclude
SRC = $(wildcard src/*.cpp) main.cpp  # Espande correttamente tutti i file .cpp in src/
OBJ = $(SRC:.cpp=.o)
EXEC = interpolate

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) -o $@ $^

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJ) $(EXEC)