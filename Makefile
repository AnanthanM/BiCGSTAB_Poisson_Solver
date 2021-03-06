NAME 	:= EXECUTABLE
OBDIR	:= OBJ/
RUNDIR	:= RUN

OBJECTS =  $(OBDIR)Helper_Functions.o $(OBDIR)VTK_Output.o $(OBDIR)Set_Ghost_Cells.o $(OBDIR)FVM.o  $(OBDIR)BiCGSTAB.o $(OBDIR)Pressure.o

HEADER 	:=	fvm.h

INCPATH	=	-I ~/ 

LIBPATH	= 	-L ~/

LIBS	=	-lm 

CC	:=	gcc

FLAGS	:=	  $(INCPATH) -g  -Wall -O3 -march=native -ftree-vectorize -ftree-vectorizer-verbose=1 -fopenmp

LFLAGS	:=	$(LIBPATH) $(LIBS) 

####################################################################

$(NAME): $(OBJECTS) $(RUNDIR)
	$(CC) -o $(NAME) $(OBJECTS) $(FLAGS) $(LFLAGS)
	@mv $(NAME) $(RUNDIR) 

$(RUNDIR): 
	@test -d $(RUNDIR) || mkdir $(RUNDIR)

$(OBDIR)FVM.o: FVM.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)FVM.o FVM.c

$(OBDIR)Set_Ghost_Cells.o: Set_Ghost_Cells.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)Set_Ghost_Cells.o Set_Ghost_Cells.c

$(OBDIR)VTK_Output.o: VTK_Output.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)VTK_Output.o VTK_Output.c

$(OBDIR)Helper_Functions.o: Helper_Functions.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)Helper_Functions.o Helper_Functions.c

$(OBDIR)BiCGSTAB.o: BiCGSTAB.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)BiCGSTAB.o BiCGSTAB.c

$(OBDIR)Pressure.o: Pressure.c $(HEADER)
	@test -d OBJ || mkdir OBJ        
	$(CC) -c $(FLAGS) -o $(OBDIR)Pressure.o Pressure.c

clean: 		
	rm -f *~
	rm $(OBDIR)*.o

Delete_Results: 		
	rm -f *~
	rm -r  $(RUNDIR)
	rm -r  $(OBDIR)
