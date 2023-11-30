# Makefile

TARGET = prova
SOURCE = main.cpp

${TARGET}: ${SOURCE}
	g++ -Wall -I${mkEigenInc} -o ${TARGET} ${SOURCE}

clean:
	@rm -f ${TARGET}

remake:
	@make clean && make