CC = gcc
CFLAGS= -O2
all:  makegp gprun conrun grrun
.c.o:
	${CC} -c $(CFLAGS) $*.c 
makegp: mg.o
	${CC} $(CFLAGS) -o $@ $<
gprun: gpd.o gpp.o permfns.o
	${CC} $(CFLAGS) -o $@ $<
conrun: conpdcd.o conpdc1.o conpdc2.o
	${CC} $(CFLAGS) -o $@ $<
grrun: grd.o grp.o permfns.o
	${CC} $(CFLAGS) -o $@ $<

.PHONY: all
