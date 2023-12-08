# Makefile

TARGET = prova
SOURCE = main.cpp

${TARGET}: ${SOURCE}
	g++ -Wall -O3 -I${mkEigenInc} -o ${TARGET} ${SOURCE}

clean:
	@rm -f ${TARGET}

remake:
	@make clean && make