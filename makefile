HEADER_FILES = short_seq_mapper.h layer.h
CPP_FILES = main.cpp short_seq_mapper.cpp layer.cpp
EXECUTABLE = short_seq_mapper

all: main run

main: $(HEADER_FILES) $(CPP_FILES)
	@g++ -std=c++11 -O3  -o $(EXECUTABLE) $^

.PHONY: run
run:
	@$(EXECUTABLE)
