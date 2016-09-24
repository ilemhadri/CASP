#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <mpi.h>
using namespace std;


/* ////////////////////
Calcul en MPI des dRSMD entre tous les couples de prediction
Le resultat est une matrice de distances enregistree dans "MatriceDistances.txt"
//////////////////////////  */

/* creation d'une matrice de maniere contigue dans la memoire
utile pour les operations MPI */
double **alloc_2d_double(int rows, int cols){
   double* data = (double *)malloc(rows*cols*sizeof(double));
   double ** array = (double **)malloc(rows*sizeof(double));
   for (int i =0; i<rows;++i)
      array[i] = &(data[cols*i]);
   return array;
}


/* utile pour le clustering hierarchique par Single Linkage */
double cal(double a,double b){ 
   return min(a,b);
}


/* meme fonction que dans calcul_dist.cpp
Remarque : nous avons essaye de creer un header pour la fonction calcul_dist.cpp
afin d'eviter d'avoir a repeter son code dans chaque programme qui l'appelle,
toutefois il semble y avoir une incompatibilite du header avec le type string */
double dist(string s1, string s2){
   s1 = "T0651_results/"+s1;
   s2 = "T0651_results/"+s2;
   ifstream input1(s1.c_str());
   ifstream input2(s2.c_str());
   int r1, r2;//les residus courants des lignes lues dans les fichiers 1 et 2.
   string line1, line2;
   getline(input1, line1);
   getline(input2,line2);
   istringstream stream1 (line1);
   istringstream stream2 (line2);
   string dummy_string;
   stream1>>dummy_string>>dummy_string>>r1;      
   stream2>>dummy_string>>dummy_string>>r2;
   input1.clear();
   input1.seekg(0, ios::beg);
   input2.clear();
   input2.seekg(0, ios::beg);
   int max_r = max(r1,r2);
   /*calculer ce max sert a calculer le N qui apparait dans la formule de la dRMSD,
   car nous avons remarque que deux structures n'ont pas forcement le meme nombre de residus (certains
   fichiers commencent avec un residu plus grand que 1).
   */
   int N = 95 - max_r +1;//on prend donc les residus communs aux 2 fichiers.
   double ** m1 = new double*[N];
   for (int k = 0; k<N; ++k){
      m1[k] = new double[3];
   }//m1 est une matrice donc chaque ligne comporte les coordonnees du CA d'un residu, du fichier 1.
   double ** m2 = new double*[N];
   for (int k = 0; k<N; ++k){
      m2[k] = new double[3];
   }//m2 est une matrice donc chaque ligne comporte les coordonnees du CA d'un residu, du fichier 2.
   
   int compt = 0;//sert pour savoir quel residu on lit.
   while ( getline(input1,line1) ){
      while (r1<max_r){
         getline(input1,line1);
         r1 +=1;
      } /*si le fichier 2 a moins de residus que le 1(ie r1<max_r), on fait une boucle pour lire des lignes
      du fichier 1 afin de se mettre au meme niveau dans les deux fichiers.
      */
      /*
      au plus une des deux boucles while (r1<max_r) et while (r2<max_r) (plus bas) sera lue.
      */
      istringstream stream1 (line1);
      stream1>>dummy_string>>dummy_string>>r1>>m1[compt][0]>>m1[compt][1]>>m1[compt][2];
      compt += 1;  
   }

   compt = 0;
   while ( getline(input2,line2) ){
      while (r2<max_r){
         getline(input2,line2);
         r2 +=1;
      } /*si le fichier 1 a moins de residus que le 2(ie r2<max_r), on fait une boucle pour lire des lignes
      du fichier 2 afin de se mettre au meme niveau dans les deux fichiers.
      */
      istringstream stream2 (line2);
      stream2>>dummy_string>>dummy_string>>r2>>m2[compt][0]>>m2[compt][1]>>m2[compt][2];
      compt += 1;  
   }
     
    double sum = 0.;
    for (int i =0; i<N;++i){
       for (int j =i+1; j<N;++j){
          double d1 = sqrt( (m1[i][0]-m1[j][0])*(m1[i][0]-m1[j][0]) +
                (m1[i][1]-m1[j][1])*(m1[i][1]-m1[j][1])+(m1[i][2]-m1[j][2])*(m1[i][2]-m1[j][2]) );
          double d2 = sqrt( (m2[i][0]-m2[j][0])*(m2[i][0]-m2[j][0]) +
                (m2[i][1]-m2[j][1])*(m2[i][1]-m2[j][1])+(m2[i][2]-m2[j][2])*(m2[i][2]-m2[j][2]) );
          sum += (d1 - d2)*(d1 - d2);
       }
    }
    double r = ((double) N)*(((double) N) - 1.);
    return sum/r;
}

//////////////////////////////

int main(int argc, char* argv[]){   
   const int N = 456;
         
   string * predictionId = new string[N];
   ifstream filenames("filenames_results.txt");
   int i = 0;
   string filename;
   while ( getline(filenames,filename)  ){
      predictionId[i] = filename;
      i +=1;
   }
   
   double **D = alloc_2d_double(N,N); //Matrice des distances dRMSD   
   int taskid, numtasks;
   MPI_Init(&argc,&argv);
   MPI_Status status;
   MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
   MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
   
   int myLength[numtasks];
   for (int i=0;i<numtasks;++i){
      myLength[i] = N/numtasks;
   }
   myLength[numtasks-1] += N%numtasks;

   int myStart = taskid*(N/numtasks);
   
   double**myMatrix = alloc_2d_double(myLength[taskid],N);
   
   if (taskid ==0){   

	for (int i=0;i<myLength[taskid];++i){
         for (int j=0;j<N;++j){
            D[i][j] = dist(predictionId[i], predictionId[j]);
         }
        }
	
	for (int i =1;i<numtasks;++i){
         MPI_Recv( &(D[i*N/numtasks][0]), myLength[i]*N, MPI_DOUBLE, i, 0,
            MPI_COMM_WORLD,&status);   
        }

      ofstream output("MatriceDistances.txt");
      for (int i=0; i<N;++i){
         for (int j=0; j<N;++j){
            output << D[i][j] << " ";
         }
         output << "\n";
      }  
  }
   
   if (taskid > 0){
     for (int i = 0; i<myLength[taskid];++i){
         for (int j=0; j<N;++j){
            myMatrix[i][j] = dist(predictionId[myStart+i],predictionId[j]);
         }
      } 
      MPI_Send( &(myMatrix[0][0]),myLength[taskid]*N,MPI_DOUBLE,0,0,MPI_COMM_WORLD);     
   }
   
   //////////////
   
   MPI_Finalize();
   return 0;
}
