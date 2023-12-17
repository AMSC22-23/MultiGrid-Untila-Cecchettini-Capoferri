# Makefile

TARGET = prova
SOURCE = ./src/*

${TARGET}: ${SOURCE}
	g++ -Wall -O3 -I ./include -o ${TARGET} ${SOURCE}

clean:
	@rm -f ${TARGET}

remake:
	@rm -f ${TARGET} && make