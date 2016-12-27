# MpiLineDetection

CMPE300 PROGRAMMING PROJECT
2016400372 MEHMET SEFA BALIK

Smoothing and line detecting algorithm using
C programming language with the help of MPI (Message Passing Interface) libraries.

27.12.2016






Introduction

The goal of this project is to detect lines in a grayscale image. The program :
•	Takes 200x200 grayscale image array(should be given as 1. Command line argument)
•	Master processor reads and scatters divided arrays to the slaves 
•	Slaves store divided arrays
•	Slaves communicate  with each other, sends and receives the boundary rows
•	When communication is done, Slaves smooth the arrays in parallel
•	When smoothing is done, slaves again communicate with each other, sends and receives the boundary smoothed rows
•	After communication is done, slaves do the convolution in parallel to detect horizontal, vertical and oblique lines, and treshold with the given number(as 3. Command line argument)
•	After convolution is done, slaves gather convolved arrays in Master processor.
•	After gathering, Master writes the convolved array in output file(should be given as 2. Command line argument)








Program Interface
Compile:
•	mpicc <project c file> -o <projectexec>
e.g. 	mpicc project.c –o project
Run:
•	mpiexec –n  <number of processors>  <projectexec> <inputfile> <outputfile> <treshold>
e.g.	mpiexec –n 101 project input.txt output10.txt 10

Program Execution:
Compile:
 <project c file>  is the name of the .c project file . <projectexec> is the name of the executable.
Run:
 <number of processors> is the number of processors that you want to divide plus master.  <projectexec> is the name of the executable that is compiled beforehand. <inputfile>  should contain the 200x200 grayscale image array. <outputfile>  is the name of the file that you want to write the output. <treshold> is the integer number that would be the treshold for the output image.





Program Structure
C Programming language is used to implement the project. For parallel programming, Mpi library is used. Since C is a procedurel language, project is implemented highly procedurel. Apart from Mpi library functions, there are no functions, no objects. To me, procedurel style is easier to implement when using  mpi library. 

•	Mpi Functions That Are Used:

MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,  MPI_Comm comm)
//Sends data from one process to all other processes 
MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
	//Gathers together values from a group of processes
MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,MPI_Comm comm);
             //Performs a blocking send
  MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,MPI_Comm comm,MPI_Status_IGNORE)
             //Blocking receive for a message
  MPI_Init(&argc, &argv)
 //Initialize the MPI execution environment 

  MPI_Comm_rank(MPI_COMM_WORLD, &rank)
//Determines the rank of the calling process in the communicator 
  MPI_Comm_size(MPI_COMM_WORLD, &size)
//Determines the size of the group associated with a communicator 
  MPI_Barrier(MPI_COMM_WORLD
 //Blocks until all processes in the communicator have reached this routine. 
  MPI_Finalize()
//Terminates MPI execution environment

•	Data Structures That Are Used:

Only 1D  integer arrays are used as data structure. I used 1D arrays instead of 2D for simplicity. I used 2D array logic in 1D integer. For Example, instead of array[i][j] , I used array[i +N*j]  ,(N is the number of coloumns).

 

Conclution:
	To conclude,  the goal of the program is to detect lines in a grayscale image that is given to the program as parameter. The program uses mpi library for parallel programming. I tried the program with 200x200 gray scale image input and 101 processors.  Program is implemented in C and with highly procedurel approach. I also highly commented the code.

