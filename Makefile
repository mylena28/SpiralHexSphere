FC     = gfortran
FFLAGS = -Wall -Wextra -O2
TARGET = spiral
DATA   = spiral.dat

OBJS = params.o geometry.o io.o main.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

params.o: params.f90
	$(FC) $(FFLAGS) -c $< -o $@

geometry.o: geometry.f90 params.o
	$(FC) $(FFLAGS) -c $< -o $@

io.o: io.f90 params.o
	$(FC) $(FFLAGS) -c $< -o $@

main.o: main.f90 geometry.o io.o params.o
	$(FC) $(FFLAGS) -c $< -o $@

run: $(TARGET)
	./$(TARGET)

plot: $(DATA)
	python3 plot_mesh.py

all_run: $(TARGET)
	./$(TARGET)
	python3 plot_mesh.py

clean:
	rm -f *.o *.mod $(TARGET) $(DATA)

.PHONY: all run plot all_run clean
