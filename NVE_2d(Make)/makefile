compile = gcc -O1
link = -lm

info:
	@echo "Usage:"
	@echo "  make run       # Clean, compile, and run the program"
	@echo "  make execute   # Just run the compiled program"
	@echo "  make clean     # Remove all object files and executable"

# Object files
global.o: global.c global.h subroutine.h
	$(compile) -c global.c -o global.o

initialize.o: initialize.c global.h subroutine.h
	$(compile) -c initialize.c -o initialize.o

cell_list.o: cell_list.c global.h subroutine.h
	$(compile) -c cell_list.c -o cell_list.o

time_evolve.o: time_evolve.c global.h subroutine.h
	$(compile) -c time_evolve.c -o time_evolve.o

main.o: main.c global.h subroutine.h
	$(compile) -c main.c -o main.o

# List of object files
object = global.o initialize.o cell_list.o time_evolve.o main.o

# Run target: Cleans, compiles, and then executes
run: clean $(object) 
	$(compile) $(object) $(link) -o run
	@echo "Compilation complete...."

# Execute target: Runs the program without recompiling
execute: 
	@echo "Running existing compiled program..."
	@./run

# Clean target: removes object files and executable
clean:
	@rm -fv $(object) run
	@echo "All object files and executable have been removed."
