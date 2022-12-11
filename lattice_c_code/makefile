CC=gcc
CFLAGS=-Wall -pedantic -std=c99 -g3

all: Main

Main: main.o helper.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

main.o: helper.h
#//helper.o: helper.h
