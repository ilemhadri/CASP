#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;



void results(string s){
   string s1;
   if (s.compare("target.pdb")==0){
      s1 = s;
   }else{
      s1 = "T0651/"+s;
   }
   // T0651 : dossier contenant les fichiers de prediction
   ifstream input(s1.c_str());
   ofstream myfile;
   string s2 = "T0651_results/"+s+"_result.txt";
   /* T0651_results : dossier contenant les resultats (pour la cible et les
   predictions) */
   myfile.open(s2.c_str());
   int dummy_int;

   string line;
   while( getline(input,line) ){
      istringstream stream (line);
      double x,y,z;
      string dummy_string; // contient les informations non necessaires
      string atom_string;
      string name;
      int residue;
      string residue_name;
      stream>>atom_string>>dummy_string>>name>>residue_name>>residue>>x>>y>>z;
      // on recupere les carbones alpha du domaine 1
      if (atom_string == "ATOM" && name=="CA" && (residue<=95) ){
         myfile <<name+" "<<residue_name+" "<<residue<<" "<<x<<" "<<y<<" "<<z<<endl;
      }
   }
}




int main(int argc, char * argv[]){
   ifstream filenames("filenames.txt");
   // filenames contient tous les noms de fichiers a transformer
   string filename;
   while( getline(filenames,filename) ){
      results(filename);
   }
   results("target.pdb");
   return 0;
}
         
      
