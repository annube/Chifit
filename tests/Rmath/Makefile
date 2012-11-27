
R_HOME=/afs/ifh.de/group/nic/scratch/poole/projects_andreas/tools/R-2.15.1

LDFLAGS=-L${R_HOME}/src/nmath/standalone -lRmath
CFLAGS=-I${R_HOME}/include


rmex1: rmex1.cc Makefile
	${CXX} ${CFLAGS} $< ${LDFLAGS} -o $@
