#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
using namespace std;
//Clustering hierarchique par Single Linkage

double cal(double a,double b){
   return min(a,b);
}

const int N = 456;
bool used[N];  //clusters qui ont ete supprimes par fusion
double D[N][N];   // matrice des distances

int main(){
   // Recuperer la matrice des distances
   ifstream input("MatriceDistances.txt");
   string line;
   for (int i=0;i<N;++i){
      getline(input,line);
      istringstream stream (line);
      for (int j=0;j<N;++j){
      stream >> D[i][j];
      }
   }
   
   int num_clusters = N; 
   int cluster[N];   
   int num_clusters_end = 90; // nombre de clusters a la fin du regroupement
   for (int i=0;i<N;++i)
      cluster[i] = i;

   while(num_clusters > num_clusters_end){
      int bi = -1,bj;
      for(int i=0;i<N;++i) {
			if(used[i]) continue;
			for(int j=i+1;j<N;++j) {
				if(used[j]) continue;
				if(bi == -1 || D[i][j] < D[bi][bj]) {
					bi = i;
					bj = j;
				}
			}
		}      
// la paire la plus rapprochee est (bi,bj), on la fusionne en nc = min(bi,bj)
   int nc = min(bi,bj); //nouveau cluster
   int m = max(bi,bj);
   //reattribuer les points du cluster m a nc
   for (int i=0;i<N;++i){
      if (cluster[i]==m)
            cluster[i]=nc;
   }
   used[m] = 1;
   --num_clusters;
   if (num_clusters ==1) break;
   //mise a jour de la distance de nc aux autres clusters
   for (int i =0; i<N;++i){
      if (used[i]) continue;
      D[nc][i]=D[i][nc] = cal(D[bi][i], D[bj][i]);
      }
   }
   
   // affichage des clusters
   for (int i=0;i<N;++i){
      bool non_empty_cluster = 0;
      for (int j =0; j<N;++j){
         if (cluster[j] == i){
            if (non_empty_cluster == 0) {
               non_empty_cluster =1;
               cout<<"\n Les points du cluster "<<i<<" sont :"<<j;
            }
            else
               cout<<" , "<<j;
         }
      }
   }
   cout <<endl;

   /////////
   /*affichage des clusters des dix meilleures predictions*/
   /////////
   /*recuperer les noms des fichiers*/
   int DataSize = 425;
   string * predictionId = new string[DataSize];
   ifstream filenames("filenames_results95.txt");
   int i = 0;
   string filename;
   while ( getline(filenames,filename) ){
      predictionId[i] = filename;
      i +=1;
   }
   /*recuperer les identifiants des dix meilleures predictions*/
   int* best_10_ID = new int[10];
   ifstream predictions_ID("10_Best95.txt");
   /*afficher les dix meilleures predictions et leurs clusters*/
   i =0;
   string ID;
   cout << "\n  Les dix meilleures predictions et leurs clusters:\n";
   while( getline(predictions_ID,ID)){
      istringstream stream (ID);
      stream >> best_10_ID[i];
      cout << predictionId[best_10_ID[i]] << " "<<cluster[best_10_ID[i]]<<endl;
      i+=1;
   }
   
}
