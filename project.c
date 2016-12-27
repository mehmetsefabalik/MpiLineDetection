#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"


/*
*Master reads the 200x200 array, divides the array size-1 rows, sends divided arrays to the slaves
*Slaves get boundary rows, smooth arrays
*Slaves share smoothed boundary arrays, apply convolving and tresholding
*Slaves send convolved arrays to Master
*Master concatanetes arrays and writes final array to the output file
*
*/
int main(int argc, char* argv[])
{
  if(argc<=3) {
      printf("BAD USAGE! Use this way: mpiexec -n 101 project input.txt 25");
      exit(1);
   }

    int rank; // rank of the current processor
    int size; // total number of processors

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // gets the rank of the current processor
    MPI_Comm_size(MPI_COMM_WORLD, &size); // gets the total number of processors

    int N = 200;
    FILE *cin = fopen(argv[1], "r"); // reads from argument
    FILE *f = fopen(argv[2], "w");
    int arr[N*N]; // stores preferences of all
    int all[(N-4)*(N-4)];
    int a = N/(size-1);
    int array[N*N+(N*a)];
    int pref[N*a+2*N]; // stores preferences of each
    int smoothEach[(N-2)*a + (2*(N-2))];
    int convolvedEach[(N-4)*a + (2*(N-4))];
    int sendTop[N];
    int receiveBottom[N];
    int sendBottom[N];
    int receiveTop[N];
    int topBottom = 0;
    int sendTopC[N-2];
    int receiveBottomC[N-2];
    int sendBottomC[N-2];
    int receiveTopC[N-2];
    int topBottomC = 0;
    int numofel;
    int snumofel;
    int cnumofel;

    // If it's master processor, reads from input file
    if(rank==0){
    	int num=0,i=0;
    	while(fscanf(cin, "%d", &num)==1){
    		arr[i]=num;
    		i++;
    	}
      for (int i = 0; i <  (N*a); i++) {
        array[i] = 0;
      }
      int w = 0;
      for (int i = N*a; i < N*N + (N*a); i++) {
        array[i] = arr[w];
        w++;
      }
    	fclose(cin);
    }

    // sends data from root array arr to pref array on each processor

    MPI_Scatter(array,N*a,MPI_INT,pref,N*a,MPI_INT,0,MPI_COMM_WORLD);

    if (rank%2==1 && rank!=size-1) {

      int b=0;
      for (int i = N*a-N; i < N*a; i++) {
          sendBottom[b]=pref[i];
          b++;
       }

        MPI_Send(&sendBottom, N, MPI_INT, rank+1, 0, MPI_COMM_WORLD);

        MPI_Recv(&receiveBottom, N, MPI_INT, rank+1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        int r=0;
        for (int i = N*a; i < N*a+N; i++) {
            pref[i]=receiveBottom[r];
            r++;
        }
         topBottom = 1;
         if (rank==1) {
           int temp2[N*a+N];
           for (int i = 0; i < N*a+N; i++) {
               temp2[i]=pref[i];
           }
           int r2=0;
           int pref[N*a+N];
           for (int i = 0; i < N*a+N; i++) {
               pref[i]=temp2[r2];
               r2++;
           }
      /*     int i = 0;
               for(; i<N*a+N;i++){
                   printf("Process Numb %d and %d th element of my list is %d\n",rank,i+1,pref[i] );
               }*/
         }
//////////////////////////////////////////////////////////////////

    }
    else if (rank%2==0 && rank!=0) {
      int b=0;
      for (int i = 0; i < N; i++) {
          sendTop[b]=pref[i];
          b++;
      }
        MPI_Recv(&receiveTop, N, MPI_INT, rank-1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        int temp[N*a+N];
        for (int i = 0; i < N*a+N; i++) {
            temp[i]=pref[i];
        }
        int r=0;
        for (int i = 0; i < N; i++) {
            pref[i]=receiveTop[r];
            r++;
        }
        r=0;
        for (int i = N; i < N*a+N; i++) {
            pref[i]=temp[r];
            r++;
        }
        MPI_Send(&sendTop, N, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
        topBottom = 1;
        if (rank==size-1) {
          int temp2[N*a+N];
          for (int i = 0; i < N*a+N; i++) {
              temp2[i]=pref[i];
          }
          int r2=0;
          int pref[N*a+N];
          for (int i = 0; i < N*a+N; i++) {
              pref[i]=temp2[r2];
              r2++;
          }
        /*  int i = 0;
              for(; i<N*a+N;i++){
                  printf("Process Numb %d and %d th element of my list is %d\n",rank,i+1,pref[i] );
              }*/
        }
////////////////////////////////////////////////////////////////////////////////////////

    }
    MPI_Barrier(MPI_COMM_WORLD);
//***********************************************************************************
    if(rank%2==1 && topBottom == 1 && rank!=1){
      int b1=0;
      for (int i = 0; i < N; i++) {
          sendTop[b1]=pref[i];
          b1++;
       }

        MPI_Send(&sendTop, N, MPI_INT, rank-1, 0, MPI_COMM_WORLD);

        MPI_Recv(&receiveTop, N, MPI_INT, rank-1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

      int temp1[N*a+N];
      for (int i = 0; i < N*a+N; i++) {
          temp1[i]=pref[i];
      }
      int r1=0;
      for (int i = 0; i < N; i++) {
          pref[i]=receiveTop[r1];
          r1++;
      }
      r1=0;
      for (int i = N; i < N*a+2*N; i++) {
          pref[i]=temp1[r1];
          r1++;
      }


    /*  int i = 0;
          for(; i<N*a+2*N;i++){
              printf("Process Numb %d and %d th element of my list is %d\n",rank,i+1,pref[i] );
          }*/
    }

    else if(rank%2==0 && rank!=0 && topBottom == 1 && rank!=size-1){
      int b1=0;
      for (int i = N*a; i < N*a+N; i++) {
          sendBottom[b1]=pref[i];
          b1++;
      }

        MPI_Recv(&receiveBottom, N, MPI_INT, rank+1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        int r1=0;
        for (int i = N*a+N; i < N*a+2*N; i++) {
            pref[i]=receiveBottom[r1];
            r1++;
        }
        MPI_Send(&sendBottom, N, MPI_INT, rank+1, 0, MPI_COMM_WORLD);

    /*  int i = 0;
          for(; i<N*a+2*N;i++){
              printf("Process Numb %d and %d th element of my list is %d\n",rank,i+1,pref[i] );
          }*/
    }
//////////////////////////////////////////SMOOTHING///////////////////////////////////
MPI_Barrier(MPI_COMM_WORLD);
    if(rank!=0){
    if(rank!=1 && rank!=size-1){
      float sum=0.0;
      for (int j = 2; j < ((int)(sizeof(pref)/sizeof(pref[0])/N)) ; j++) {
        for (int i = 0; i < N-2; i++) {

          sum = (pref[(j-2)*N+i]+pref[(j-2)*N+i+1]+pref[(j-2)*N+i+2]+pref[(j-2)*N+N+i]+pref[(j-2)*N+N+i+1]+pref[(j-2)*N+N+i+2]+pref[(j-2)*N+2*N+i]+pref[(j-2)*N+2*N+i+1]+pref[(j-2)*N+2*N+i+2])/9;
          smoothEach[i+(j-2)*(N-2)]=(int)sum;
        }
      }
    }else{
      float sum=0.0;
      for (int j = 2; j < 3 ; j++) {
        for (int i = 0; i < N-2; i++) {

          sum = (pref[(j-2)*N+i]+pref[(j-2)*N+i+1]+pref[(j-2)*N+i+2]+pref[(j-2)*N+N+i]+pref[(j-2)*N+N+i+1]+pref[(j-2)*N+N+i+2]+pref[(j-2)*N+2*N+i]+pref[(j-2)*N+2*N+i+1]+pref[(j-2)*N+2*N+i+2])/9;
          smoothEach[i+(j-2)*(N-2)]=(int)sum;
        }
      }

    }

    }

    /////////////////////////////////////////////////////////////////////////////
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank!=0){
      if (rank==size-1 || rank == 1) {
        int temp3[(N-2)*a-(N-2)];
        for (int i = 0; i < (N-2)*a-(N-2); i++) {
            temp3[i]=smoothEach[i];
        }
        int r3=0;
        int smoothEach[(N-2)*a-(N-2)];
        for (int i = 0; i < (N-2)*a-(N-2); i++) {
            smoothEach[i]=temp3[r3];
            r3++;
        }
        snumofel=((int)(sizeof(smoothEach)/sizeof(smoothEach[0]))/(N-2));
      }else if(rank!=size-1 && rank != 1){
        snumofel=((int)(sizeof(smoothEach)/sizeof(smoothEach[0]))/(N-2));
      }

    }

    ////////////////////////////////////////////////CONVOLUTION/////////////////////////////////
MPI_Barrier(MPI_COMM_WORLD);

    if (rank%2==1 && rank!=size-1) { // sendBottom - receiveBottom
      if(rank==1){
        int b5=0;
        for (int i = 0; i < (N-2); i++) {
            sendBottomC[b5]=smoothEach[i];
            b5++;
         }

      }else{
        int b5=0;
        for (int i = (N-2)*a-(N-2); i < (N-2)*a; i++) {
            sendBottomC[b5]=smoothEach[i];
            b5++;
         }
     }

        MPI_Send(&sendBottomC, N-2, MPI_INT, rank+1, 0, MPI_COMM_WORLD);

        MPI_Recv(&receiveBottomC, N-2, MPI_INT, rank+1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

        int r5=0;
        for (int i = (N-2)*a; i < (N-2)*a+(N-2); i++) {
            smoothEach[i]=receiveBottomC[r5];
            r5++;
        }
        if(rank!=1){
         topBottomC = 1;
        }

//////////////////////////////////////////////////////////////////

    }
    else if (rank%2==0 && rank!=0) { //sendTop - receiveTop
      int b6=0;
      for (int i = 0; i < N-2; i++) {
          sendTopC[b6]=smoothEach[i];
          b6++;
      }

        MPI_Recv(&receiveTopC, N-2, MPI_INT, rank-1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        int tempC5[(N-2)*a+(N-2)];
        for (int i = 0; i < (N-2)*a+(N-2); i++) {
            tempC5[i]=smoothEach[i];
        }
        int r6=0;
        for (int i = 0; i < N-2; i++) {
            smoothEach[i]=receiveTopC[r6];
            r6++;
        }
        r6=0;
        for (int i = N-2; i < (N-2)*a+(2*(N-2)); i++) {
            smoothEach[i]=tempC5[r6];
            r6++;
        }


        MPI_Send(&sendTopC, N-2, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
      if(rank!=size-1){
        topBottomC = 1;

        }

////////////////////////////////////////////////////////////////////////////////////////
    }
//***********************************************************************************
MPI_Barrier(MPI_COMM_WORLD);
    if(rank%2==1 && topBottomC == 1 && rank!=1){ // sendTop - receiveTop
      int b7=0;
      for (int i = 0; i < N-2; i++) {
          sendTopC[b7]=smoothEach[i];
          b7++;
       }

        MPI_Send(&sendTopC, N-2, MPI_INT, rank-1, 0, MPI_COMM_WORLD);

        MPI_Recv(&receiveTopC, N-2, MPI_INT, rank-1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

      int tempC6[(N-2)*a+(N-2)];
      for (int i = 0; i < (N-2)*a+(N-2); i++) {
          tempC6[i]=smoothEach[i];
      }
      int r7=0;
      for (int i = 0; i < N-2; i++) {
          smoothEach[i]=receiveTopC[r7];
          r7++;
      }
      r7=0;
      for (int i = N-2; i < (N-2)*a+(2*(N-2)); i++) {
          smoothEach[i]=tempC6[r7];
          r7++;
      }


    /*  int i = 0;
          for(; i<(N-2)*a+(2*(N-2));i++){
              printf("Process Numb %d and %d th element of my list is %d\n",rank,i+1,smoothEach[i] );
          }*/
    }

    else if(rank%2==0 && rank!=0 && topBottomC == 1 && rank!=size-1){ // sendBottom - receiveBottom
      int b8=0;
      for (int i = (N-2)*a; i < (N-2)*a+(N-2); i++) {
          sendBottomC[b8]=smoothEach[i];
          b8++;
      }

        MPI_Recv(&receiveBottomC, N-2, MPI_INT, rank+1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        int r8=0;
        for (int i = (N-2)*a+(N-2); i < (N-2)*a+(2*(N-2)); i++) {
            smoothEach[i]=receiveBottomC[r8];
            r8++;
        }
        MPI_Send(&sendBottomC, N-2, MPI_INT, rank+1, 0, MPI_COMM_WORLD);

    /*  int i = 0;
          for(; i<(N-2)*a+(2*(N-2));i++){
              printf("Process Numb %d and %d th element of my list is %d\n",rank,i+1,smoothEach[i] );
          }*/
    }

MPI_Barrier(MPI_COMM_WORLD);

    if(rank!=0 && rank!=1 && rank!=size-1){
      int t1, t2, t3, t4;
      int treshold = atoi(argv[3]);
      printf("treshold is %d\n",treshold);
      for (int j = 2; j < ((int)(sizeof(smoothEach)/sizeof(smoothEach[0]))/(N-2)) ; j++) {
        for (int i = 0; i < N-4; i++) {

        t1 = (smoothEach[(j-2)*(N-2)+i]*(-1)) + (smoothEach[(j-2)*(N-2)+i+1]*(-1)) + (smoothEach[(j-2)*(N-2)+i+2]*(-1))
             + (smoothEach[(j-2)*(N-2)+(N-2)+i]*2) + (smoothEach[(j-2)*(N-2)+(N-2)+i+1]*2) + (smoothEach[(j-2)*(N-2)+(N-2)+i+2]*2)
               + (smoothEach[(j-2)*(N-2)+2*(N-2)+i]*(-1)) + (smoothEach[(j-2)*(N-2)+2*(N-2)+i+1]*(-1)) + (smoothEach[(j-2)*(N-2)+2*(N-2)+i+2]*(-1));

        t2 = (smoothEach[(j-2)*(N-2)+i]*(-1)) + (smoothEach[(j-2)*(N-2)+i+1]*(2)) + (smoothEach[(j-2)*(N-2)+i+2]*(-1))
             + (smoothEach[(j-2)*(N-2)+(N-2)+i]*(-1)) + (smoothEach[(j-2)*(N-2)+(N-2)+i+1]*2) + (smoothEach[(j-2)*(N-2)+(N-2)+i+2]*(-1))
               + (smoothEach[(j-2)*(N-2)+2*(N-2)+i]*(-1)) + (smoothEach[(j-2)*(N-2)+2*(N-2)+i+1]*(2)) + (smoothEach[(j-2)*(N-2)+2*(N-2)+i+2]*(-1));

        t3 = (smoothEach[(j-2)*(N-2)+i]*(-1)) + (smoothEach[(j-2)*(N-2)+i+1]*(-1)) + (smoothEach[(j-2)*(N-2)+i+2]*(2))
             + (smoothEach[(j-2)*(N-2)+(N-2)+i]*(-1)) + (smoothEach[(j-2)*(N-2)+(N-2)+i+1]*2) + (smoothEach[(j-2)*(N-2)+(N-2)+i+2]*(-1))
               + (smoothEach[(j-2)*(N-2)+2*(N-2)+i]*(2)) + (smoothEach[(j-2)*(N-2)+2*(N-2)+i+1]*(-1)) + (smoothEach[(j-2)*(N-2)+2*(N-2)+i+2]*(-1));

        t4 = (smoothEach[(j-2)*(N-2)+i]*(2)) + (smoothEach[(j-2)*(N-2)+i+1]*(-1)) + (smoothEach[(j-2)*(N-2)+i+2]*(-1))
            + (smoothEach[(j-2)*(N-2)+(N-2)+i]*(-1)) + (smoothEach[(j-2)*(N-2)+(N-2)+i+1]*2) + (smoothEach[(j-2)*(N-2)+(N-2)+i+2]*(-1))
              + (smoothEach[(j-2)*(N-2)+2*(N-2)+i]*(-1)) + (smoothEach[(j-2)*(N-2)+2*(N-2)+i+1]*(-1)) + (smoothEach[(j-2)*(N-2)+2*(N-2)+i+2]*(2));

        if(t1<treshold && t2<treshold && t3<treshold && t4<treshold){
          convolvedEach[i+(j-2)*(N-4)]=0;
        }else{
          convolvedEach[i+(j-2)*(N-4)]=255;
        }
      }
    }
    /*int i = 0;
        for(; i<(N-4)*a;i++){
            printf("CONVOLVED Process Numb %d and %d th element of my list is %d\n",rank,i+1,convolvedEach[i] );
        }*/
  }

MPI_Barrier(MPI_COMM_WORLD);

//  if(rank!=1 && rank!=size-1 && rank!=0){
    printf("MPIGATHER\n");
    MPI_Gather(&convolvedEach,(N-4)*a,MPI_INT,&all,(N-4)*a,MPI_INT,0,MPI_COMM_WORLD);
    printf("OUT*//////////////////////////////////////////////////////////////////////////////////////////n");
//  }

MPI_Barrier(MPI_COMM_WORLD);

  if(rank == 0){
    printf("writing to \n");
    for(int i=1;i<(N-4)*(N-4)+1;i++){
      //printf("Process Numb %d and %d th element of my list is %d\n",rank,i+1,all[i] );
      fprintf(f, "%d ",all[i] );
      if(i%196==0 && i!=0)
        fprintf(f, "\n");
    }
    fclose(f);
  }
  /*  int masterSignal = 1;
    while(masterSignal){

    	if(rank!= 0){
    		int i = 0;
            for(; i<N*a;i++){
              //  printf("Process Numb %d and %d th element of my list is %d\n",rank,i+1,pref[i] );
            }
    	}

    	if(rank==0){
    		masterSignal=0;
    	}

    	MPI_Bcast(&masterSignal, 1, MPI_INT, 0, MPI_COMM_WORLD);


    }*/


    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
       return 0;
}
