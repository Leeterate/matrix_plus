CC=g++
OS=$(shell uname)

ifeq ($(OS), Linux)
	LIBS=-lcheck -lgtest -lpthread -lrt -lm -lsubunit
	LEAKS=CK_FORK=no valgrind -s --tool=memcheck --leak-check=full --show-leak-kinds=all --show-reachable=yes ./test
else
	LIBS=-lcheck -lgtest -pthread
	LEAKS=leaks --atExit -- ./test
endif

all: clean matrix_oop.a test

test: clean matrix_oop.a
	$(CC) -g test.cpp matrix_oop.a $(LIBS) -o test
	./test

gcov_report: matrix_oop.a
	$(CC) --coverage -c matrix_oop.cpp
	ar rc matrix_oop.a *.o
	$(CC) test.cpp -L. matrix_oop.a -lcheck -lgtest -pthread -lpthread -lm -lgcov -o test_cov
	rm -rf *.o
	./test_cov
	gcov *.gcda
	lcov -t "gcov_report" -o gcov_report.info -c -d .
	mkdir -p ./report
	genhtml -o report gcov_report.info
	open ./report/index.html
	
matrix_oop.a:
	$(CC) -Wall -Wextra -Werror -std=c++17 -c matrix_oop.cpp -o matrix_oop.o
	ar rcs matrix_oop.a matrix_oop.o
	ranlib matrix_oop.a

leaks: clean test
	$(LEAKS)

CPPcheck:
	cppcheck --enable=all --std=c++17  --language=c++ --suppress=missingIncludeSystem *.cpp *.h

style:
	cp ../materials/linters/.clang-format ./
	clang-format -i *.cpp *.h
	clang-format -n *.h *.cpp
	rm -rf .clang-format

clean:
	rm -rf *.0 *.a *.gcno *.gcda *.info *.out report *.dSYM *.gch *.o test2.out .vscode .DS_Store test
	rm -rf *.gcov test_cov CPPLINT.cfg matrix_oop ../.DS_Store ../.vscode a.out