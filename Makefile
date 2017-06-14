CC=g++
CFLAGS=-O3
LDFLAGS=-static
PROGRAM=icosphere FibonacciSphere RandomSphereSampling UniformSphereSampling

all: ${PROGRAM}

icosphere: icosphere.cpp icosphere.hpp MathTools.hpp
	${CC} ${CFLAGS} icosphere.cpp -o icosphere ${LDFLAGS}

FibonacciSphere: FibonacciSphere.cpp FibonacciSphere.hpp RandomSphereSampling.hpp MathTools.hpp
	${CC} ${CFLAGS} FibonacciSphere.cpp -o FibonacciSphere ${LDFLAGS}

RandomSphereSampling: RandomSphereSampling.cpp RandomSphereSampling.hpp MathTools.hpp
	${CC} ${CFLAGS} RandomSphereSampling.cpp -o RandomSphereSampling ${LDFLAGS}

UniformSphereSampling: UniformSphereSampling.cpp UniformSphereSampling2.hpp RandomSphereSampling.hpp MathTools.hpp
	${CC} ${CFLAGS} UniformSphereSampling.cpp -o UniformSphereSampling ${LDFLAGS}

clean:
	rm ${PROGRAM}
