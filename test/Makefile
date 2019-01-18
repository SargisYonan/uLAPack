cc_options=-std=c99 -Wall -Wextra -Werror

static_test_exec=ulapack_test_static
static_test_options=-DULAPACK_USE_STATIC_ALLOC=1

dynamic_test_exec=ulapack_test_dynamic
dynamic_test_options=-DULAPACK_USE_DYNAMIC_ALLOC=1

src=../src

test : $(static_test_exec) $(dynamic_test_exec)
	./$(static_test_exec)
	./$(dynamic_test_exec)

all : build_static build_dynamic

$(static_test_exec) : static_unit_tests.o ulapack_static.o
	cc $(cc_options) unit_tests.o ulapack.o -o $(static_test_exec)

$(dynamic_test_exec) : dynamic_unit_tests.o ulapack_dynamic.o
	cc $(cc_options) unit_tests.o ulapack.o -o $(dynamic_test_exec)

dynamic_unit_tests.o : unit_tests.c ulapack_dynamic.o
	cc $(cc_options) -I$(src) -c unit_tests.c $(dynamic_test_options)

static_unit_tests.o : unit_tests.c ulapack_static.o
	cc $(cc_options) -I$(src) -c unit_tests.c $(static_test_options)

ulapack_static.o : $(src)/ulapack.c $(src)/ulapack.h $(src)/ulapack_options.h $(src)/ulapack_type.h
	cc $(cc_options) -c $(src)/ulapack.c $(static_test_options)

ulapack_dynamic.o : $(src)/ulapack.c $(src)/ulapack.h $(src)/ulapack_options.h $(src)/ulapack_type.h
	cc $(cc_options) -c $(src)/ulapack.c $(dynamic_test_options)

clean :
	rm *.o $(dynamic_test_exec) $(static_test_exec)