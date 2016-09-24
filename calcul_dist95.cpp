#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
using namespace std;

double dist(std::string s1, std::string s2){
   s1 = "T0651_results95/"+s1;
   s2 = "T0651_results95/"+s2;
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
   int max_r = max(r1,r2); /*calculer ce max sert a calculer le N qui apparait dans la formule de la dRMSD,
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

int main(){
   double x = dist("T0651TS428_1_result95.txt","target.pdb_result95.txt");
   //double x = dist("T0651TS130_2_result95.txt","T0651TS043_3_result95.txt");
   cout << x<<endl;
   return -1;
}

