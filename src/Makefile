TARGET = shinchan.exe 
OBJECTS = calc.o  rmsd_mod.o \
          PDF_Mod.o VCV_Mod.o \
          periodic.o reading_mod.o \
          main.o 
MOD_FILES  = calculation.mod      \
						 read_assign_pdf.mod  \
						 rmsd_module.mod      \
						 periodic.mod         \
						 reading.mod          \
						 vcv_mod.mod
FC = gfortran

#for gfotran 
ifeq (${FC},gfortran)
  FFLAGS += -fimplicit-none -fbounds-check
  #LDFLAGS += -fopenmp -llapack -lblas
endif

.SUFFIXES : .o .f90
.f90.o:
	${FC} -c $<

${TARGET}:${OBJECTS}
	${FC} -o $@ ${OBJECTS} ${LDFLAGS} ${FFLAGS}


.PHONY: clean
clean:
	${RM} ${TARGET} ${OBJECTS} ${MOD_FILES}
