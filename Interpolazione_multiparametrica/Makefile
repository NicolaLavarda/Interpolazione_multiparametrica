# Compilatore e opzioni
CXX = c++
CXXFLAGS = -O2 -Iinclude		#-Wall -O2 -Iinclude

# Obiettivo principale
app: main.o src/file.o src/matrix.o src/interpolation.o
	$(CXX) $(CXXFLAGS) -o app main.o src/file.o src/matrix.o src/interpolation.o

# Regole per la compilazione dei file oggetto
main.o: main.cpp include/interpolating_function.h
	$(CXX) $(CXXFLAGS) -c main.cpp

src/file.o: src/file.cpp include/file.h
	$(CXX) $(CXXFLAGS) -c src/file.cpp -o src/file.o

src/matrix.o: src/matrix.cpp include/matrix.h
	$(CXX) $(CXXFLAGS) -c src/matrix.cpp -o src/matrix.o

src/interpolation.o: src/interpolation.cpp include/interpolation.h
	$(CXX) $(CXXFLAGS) -c src/interpolation.cpp -o src/interpolation.o

# Pulizia dei file generati
clean:
	rm -f *.o src/*.o app
