HEADER_FILES = short_read_mapper.h layer.h bml_selector.h
CPP_FILES = main.cpp short_read_mapper.cpp layer.cpp bml_selector.cpp
EXECUTABLE = short_read_mapper

all: main run

main: $(HEADER_FILES) $(CPP_FILES)
	g++ -std=c++11 -O3  -o $(EXECUTABLE) $^

.PHONY: run
run:
	./$(EXECUTABLE)

.PHONY: clean
clean:
	@rm -f *.hex *.dat
	@rm -f $(EXECUTABLE)
