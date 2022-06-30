# FICHIER DE CREATION D'OBJETS ET D'UN EXECUTABLE
#     RANGEMENT PAR ORDRE ALPHABETIQUE
#
# VERSION DEBBUGGEUR -db
# VERSION OPTIMISEUR   
F90= gfortran #ifort
OBJDIR = obj1D
MODDIR = mod1D
BINDIR = bin1D
SRCDIR = Src
#cons
EXE    = main


FFLAGS = -DLINUX -J$(MODDIR)  -cpp -c $(OPT) -ffree-line-length-none -fcheck=all -Wall
LDFLAGS= -J$(MODDIR)  -cpp $(OPT)  -ffree-line-length-none -fcheck=all -Wall




OBJS = $(addprefix $(OBJDIR)/, global.o  precisions.o interface.o io.o scheme.o)


code: $(MODDIR) $(OBJDIR) $(BINDIR) $(OBJS) $(SRCDIR)/RDeuler1DHighOrder.f90 
	$(F90) $(LDFLAGS) -o $(BINDIR)/$(EXE) $(SRCDIR)/RDeuler1DHighOrder.f90 $(OBJS)



$(MODDIR): 
	mkdir -p $(MODDIR)

$(OBJDIR):	
	mkdir -p $(OBJDIR)

$(BINDIR):	
	mkdir -p $(BINDIR)	

$(OBJDIR)/global.o: $(SRCDIR)/global.f90 $(OBJDIR)/precisions.o
	$(F90) $(FFLAGS) -o $(OBJDIR)/global.o $(SRCDIR)/global.f90

$(OBJDIR)/precisions.o: $(SRCDIR)/precisions.f90
		       $(F90) $(FFLAGS) -o $(OBJDIR)/precisions.o $(SRCDIR)/precisions.f90


$(OBJDIR)/interface.o: $(SRCDIR)/interface.f90 $(OBJDIR)/precisions.o
	$(F90) $(FFLAGS) -o $(OBJDIR)/interface.o $(SRCDIR)/interface.f90

$(OBJDIR)/io.o: $(SRCDIR)/io.f90 $(OBJDIR)/precisions.o
	$(F90) $(FFLAGS) -o $(OBJDIR)/io.o $(SRCDIR)/io.f90


$(OBJDIR)/scheme.o: $(SRCDIR)/scheme.f90 $(OBJDIR)/precisions.o
	$(F90) $(FFLAGS) -o $(OBJDIR)/scheme.o $(SRCDIR)/scheme.f90


#$(OBJDIR)/elements_1D.o: $(SRCDIR)/elements_1D.f90 $(OBJDIR)/precisions.o
#	$(F90) $(FFLAGS) -o $(OBJDIR)/elements_1D.o $(SRCDIR)/elements_1D.f90




clean:	
	rm -rf $(OBJDIR)
	rm -rf $(MODDIR)
	rm $(SRCDIR)/*.f90~ 
	rm *.mod



