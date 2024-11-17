tests:
	gcc -Wall -Wextra -std=gnu23 -fdiagnostics-color \
		-g -Og -fsanitize=address,undefined \
		-I include \
		-o test \
		src/quat.c src/test.c \
		-lm
