CC = gcc
CFLAGS= -O -fno-builtin
BIN= ../exec.d
gap:    ${BIN}/gprun ${BIN}/egrun ${BIN}/sylrun ${BIN}/normrun ${BIN}/conrun\
	${BIN}/pcrun ${BIN}/scrun ${BIN}/selgen ${BIN}/matcalc ${BIN}/nqmrun\
	${BIN}/nqrun ${BIN}/grrun ${BIN}/extprun ${BIN}/crrun ${BIN}/readrels\
	${BIN}/cohomology.gap ${BIN}/execcmd.gap ${BIN}/calcpres.gap
all:  ${BIN}/makegp ${BIN}/gprun ${BIN}/egrun ${BIN}/sylrun ${BIN}/normrun\
	${BIN}/sylnorm ${BIN}/conrun ${BIN}/pcrun ${BIN}/scrun\
	${BIN}/selgen ${BIN}/chbrun ${BIN}/optrun ${BIN}/matcalc\
	${BIN}/matperm ${BIN}/nqmrun ${BIN}/nqrun ${BIN}/grrun\
	${BIN}/extprun ${BIN}/crrun ${BIN}/readrels ${BIN}/testmf\
	${BIN}/testpf ${BIN}/subdir ${BIN}/wreath ${BIN}/permmat\
	${BIN}/tcrun ${BIN}/tcrunb ${BIN}/pkrun ${BIN}/tens ${BIN}/matsum\
	${BIN}/cohomology ${BIN}/execcmd ${BIN}/clrgrp\
	${BIN}/cohomology.gap ${BIN}/execcmd.gap ${BIN}/calcpres.gap
.c.o:
	${CC} -c $(CFLAGS) $*.c

${BIN}/%:
	${CC} $(LDFLAGS) -o $@ $^

${BIN}/makegp: mg.o

${BIN}/gprun: gpd.o gpp.o permfns.o

${BIN}/egrun: egd.o egp.o permfns.o

${BIN}/sylrun: syld.o sylp.o morepermfns.o permfns.o

${BIN}/normrun: normd.o normp1.o normp2.o chb.o permfns.o morepermfns.o

${BIN}/sylnorm: snd.o sylp2.o normp1.o normp2.o permfns.o morepermfns.o chb.o

${BIN}/conrun: conpdcd.o conpdc1.o conpdc2.o

${BIN}/pcrun: pcd.o pcp.o permfns.o pcscfns.o

${BIN}/scrun: scd.o scp.o permfns.o pcscfns.o

${BIN}/selgen: slg.o

${BIN}/chbrun: testchb.o chb.o permfns.o

${BIN}/optrun: optd.o optp.o permfns.o

${BIN}/matcalc: mcd.o mcp.o matfns.o

${BIN}/matperm: matperm.o matfns.o

${BIN}/nqmrun: nqmd.o nqmp.o nqmfns.o

${BIN}/nqrun: nqd.o nqp1.o nqp2.o nqp3.o nqfns.o matfns.o

${BIN}/nqrunb: nqdb.o nqp1b.o nqp2b.o nqp3b.o nqfnsb.o matfnsb.o

${BIN}/grrun: grd.o grp.o permfns.o

${BIN}/extprun: extpd.o exa.o exb.o exc.o permfns.o matfns.o

${BIN}/crrun: crd.o crp1.o crp2.o matfns.o

${BIN}/readrels: readrels.o

${BIN}/testmf: testmf.o matfns.o

${BIN}/testpf: testpf.o permfns.o

${BIN}/subdir: subdir.o

${BIN}/wreath: wreath.o

${BIN}/permmat: permmat.o

${BIN}/tcrun: tcd.o tcp.o

${BIN}/tcrunb: tcdb.o tcpb.o

${BIN}/pkrun: pkd.o pkp.o pkmatfns.o

${BIN}/tens: tens.o

${BIN}/matsum: matsum.o


${BIN}/cohomology: cohomology
	cp cohomology $(BIN)/cohomology
	chmod +x $(BIN)/cohomology
${BIN}/execcmd: execcmd
	cp execcmd $(BIN)/execcmd
	chmod +x $(BIN)/execcmd
${BIN}/clrgrp: clrgrp
	cp clrgrp $(BIN)/clrgrp
	chmod +x $(BIN)/clrgrp
${BIN}/cohomology.gap: cohomology.gap
	cp cohomology.gap $(BIN)/cohomology.gap
	chmod +x $(BIN)/cohomology.gap
${BIN}/calcpres.gap: calcpres.gap
	cp calcpres.gap $(BIN)/calcpres.gap
	chmod +x $(BIN)/calcpres.gap
${BIN}/execcmd.gap: execcmd.gap
	cp execcmd.gap $(BIN)/execcmd.gap
	chmod +x $(BIN)/execcmd.gap

clean:
	/bin/rm -f *.o
distclean:
	/bin/rm -f *.o
	cd ../bin; /bin/rm -f *

.PHONY: clean distclean
