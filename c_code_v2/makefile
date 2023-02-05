CC=gcc
CFLAGS=-O2 -Wall -pedantic -std=c99 -g3

all: Main

Main: main.o helper.o vectors.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

main.o: helper.h vectors.h
#//helper.o: helper.h