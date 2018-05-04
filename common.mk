FC = gfortran
CC = gcc
LD = gfortran

CPPFLAGS = -I/usr/include/suitesparse
LDFLAGS = -L/usr/include/suitesparse -lumfpack -llapack -lblas -lm
