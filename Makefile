# Specify the target executable and the source files needed to build it
spkeamns:
	gcc -o spkmeans spkmeans.c -ansi -Wall -Wextra -Werror -pedantic-errors -lm  
