
# compiile option
FC = gfortran
FFLAG1 = -c -g
FFLAG2 = -lumfpack -lamd -lblas -lm

# variables
obj0 = constants.o vectors.o
obj1 = toggles.o domainmesh.o getlgmap.o prob.o\
openbasis.o localid.o eleinfo.o\
integfunction.o assembstiff.o solveequation.o geterror.o
obj2 = umfdglr4.o
obj_other = obj/umf4_f77wrapper.o

main.out: main.f90 $(obj0) $(obj1) $(obj2) $(obj_other)
	$(FC)  -o main.out  main.f90 $(obj0) $(obj1) $(obj2) $(obj_other) $(FFLAG2) 

$(obj0): %.o : %.f90
	$(FC) $(FFLAG1) $< -o $@

$(obj1): %.o : %.f90
	$(FC) $(FFLAG1) $< -o $@

$(obj2): %.o : %.f $(obj_other)
	$(FC) $(FFLAG1) $< -o $@ $(FFLAG2) 

# .PHONY
.PHONY: run clean

run: main.out
	./main.out

clean:
	- rm -f main.out *.o *.mod \
	./data_out/* ./err/* 
	clear
