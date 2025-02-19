CC=g++
CFLAGS=-Wall -Wextra -Werror -std=c++17 -lstdc++
LDFLAGS=-O2 -Wall -Wextra -Werror -std=c++17 -lstdc++
DFLAGS= -fsanitize=address -g -Wall -Wextra -Werror -std=c++17 -lstdc++
SOURCE=$(wildcard s21_*.cc)
OBJECTS=$(patsubst %.cc,%.a, ${SOURCE}) # заменяет в именах всех исходников расширение
UNAME:=$(shell uname -s) # определяем ОС
TEST_FLAGS:=-lgtest -lpthread
GCOV_FLAGS:=-ftest-coverage -fprofile-arcs

ifeq ($(UNAME),Linux ) # Linux
	OPEN_CMD=xdg-open # команда для открытия отчета
	ADD_LIB=-lm -lrt -lpthread -lsubunit # дополнительные библиотеки
	BREW=/home/darkwolf/homebrew
endif
ifeq ($(UNAME),Darwin ) # MacOS
	OPEN_CMD=open
	ADD_LIB= 
	BREW=${HOME}/homebrew
endif

PROJECTNAME=s21_matrix_oop
TESTC=test_s21_matrix_oop.cc

all: clean ${PROJECTNAME}.a test

$(PROJECTNAME).a: ${SOURCE}
	${CC} -c $(CFLAGS) ${ADD_LIB} ${SOURCE}
	ar rc $(PROJECTNAME).a $(PROJECTNAME).o
	ranlib $(PROJECTNAME).a

test: $(TESTC) ${PROJECTNAME}.a
	g++ $(CFLAGS) $(TEST_FLAGS) $(TESTC) ${PROJECTNAME}.a -o $@ 
	./$@

gcov_report:
	g++ $(CFLAGS) -c $(TESTC)
	g++ $(CFLAGS) $(GCOV_FLAGS) -c $(SOURCE)
	g++ $(CFLAGS) $(GCOV_FLAGS) *.o $(TEST_FLAGS) -o test
	./test
	lcov -t ${PROJECTNAME} -o ${PROJECTNAME}.info -c -d .;\
	genhtml -o report ${PROJECTNAME}.info;
	open report/index.html
	-$(OPEN_CMD) ./report/index.html

rebuild: clean all

# main: main.c s21_matrix.c
# 	gcc main.c s21_matrix.c -lm -g

clean:
	-rm -rf *.o && rm -rf *.gcno
	-rm -rf *.a && rm -rf *.gcda
	-rm -rf *.info && rm -rf *.gcov
	-rm -rf ./test && rm -rf ./gcov_report
	-rm -rf ./report/

valgrind: test
	valgrind --tool=memcheck --leak-check=yes -s ./test
	
leaks: test
	leaks -atExit -- ./a.out 

clang:
	cp ../materials/linters/.clang-format .clang-format
	clang-format -i *.cc *.h
	rm .clang-format	

clang_err:
	cp ../materials/linters/.clang-format .clang-format
	clang-format -n *.cc *.h
	rm .clang-format

.PHONY: all clean test