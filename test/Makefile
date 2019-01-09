cc_options=-std=c99 -g -Wall -Wextra -Werror
test_exec=ulapack_test

src=../src

all : build
	./$(test_exec)

build: unit_tests.o ulapack.o
	cc $(cc_options) unit_tests.o ulapack.o -o $(test_exec)

unit_tests.o : unit_tests.c ulapack.o
	cc $(cc_options) -I$(src) -c unit_tests.c

ulapack.o : $(src)/ulapack.c $(src)/ulapack.h $(src)/ulapack_options.h $(src)/ulapack_type.h
	cc $(cc_options) -c $(src)/ulapack.c

clean :
	rm *.o $(test_exec)