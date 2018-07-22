/*
	Elvin Rivera
	CSci 191T
	Spring 17

  - This program replicates the Smith-Waterman local alignment algorithm
    = uses linear gap scoring scheme

  - Linear gap scoring
    = match : +2
    = mismatch: -1
    = deletion: -1

  + input query and database DNA files

  - to compile in cygwin: g++ Prog4CS191ER.cpp -o "executable name"
  -  to run executable in cygwin: ./"executable name" Prog4-database.fa Prog4-query.fa
    > chr1.67101626.67101698.NM_03_2291_exon_3_0_chr1_67101627_f.+
*/



#include <iostream>
#include <string>
#include <fstream>  // for open, close files
#include <algorithm>  // for reverse string
using namespace std;

//  ------ FUNCTION PROTOTYPE
void SmithWaterman(string querySeq, string dbSeq);
void showInverted(string Qstr, string DBstr);
void traversePath(string querySeq, string dbSeq, int best_i, int best_j, int **V_table, int **T_table);

void SmithWaterman(string querySeq, string dbSeq) {

  int max = 0;
  int direction = 0;
  int **V_table;
  int **T_table;
  int optimumSmithWaterman = 0;
  int optimumSmithWaterman_i = 0;
  int optimumSmithWaterman_j = 0;

  V_table = new int*[querySeq.size() + 1];
  T_table = new int*[querySeq.size() + 1];

  for(int i = 0; i < querySeq.size() + 1; i++) {
    V_table[i] = new int[dbSeq.size() + 1];
    T_table[i] = new int[dbSeq.size() + 1];
  }

  // initialize V_table[0][col] = 0
  for(int col = 0; col < dbSeq.size(); col++) {
    V_table[0][col] = 0;
    T_table[0][col] = 0;
  }

  // initialize T_table[row][0] = 0
  for (int row = 0; row < querySeq.size(); row++) {
    T_table[row][0] = 0;
    V_table[row][0] = 0;
  }

  for(int i = 1; i < querySeq.size(); i++) {
    for(int j = 1; j < dbSeq.size(); j++) {

      // check sequence if matched add 2 to score
      // if mismatch/delete subtract 2 from score
      if(querySeq[i] == dbSeq[j]) {
        V_table[i][j] = V_table[i-1][j-1] + 2; // match
        T_table[i][j] = 1;  // score for diagonal
      }
      else {
        // obtain left value;  direction = 3
        max = V_table[i][j-1]; // mismatch
        direction = 3;

        // check diagonal with max; if larger, direction -> 1
        if(max < V_table[i-1][j-1]) {
          max = V_table[i-1][j-1];
          direction = 1;
        }
        // check top with max, if larger, direction -> 2
        if(max < V_table[i-1][j]) {
          max = V_table[i-1][j];
          direction = 2;
        }

        // if max is 0 go with left path direction
        if (max == 0) {
          V_table[i][j] = 0;
          T_table[i][j] = direction;
        }
        // otherwise -1 from max value and becomes new direction
        else {
          V_table[i][j] = max - 1;
          T_table[i][j] = direction;
        }
      }

      // obtain and save the best score avaiable
      if(V_table[i][j] > optimumSmithWaterman) {
        optimumSmithWaterman = V_table[i][j];
        optimumSmithWaterman_i = i;
        optimumSmithWaterman_j = j;
      }
    }
  }
  // display best Smith Waterman Score
  cout << "Optimum Smith-Waterman score = " << optimumSmithWaterman << endl;

  // call traversePath
  traversePath(querySeq, dbSeq, optimumSmithWaterman_i, optimumSmithWaterman_j, V_table, T_table);
}


void showInverted(string Qstr, string DBstr) {

  reverse(Qstr.begin(), Qstr.end());
  reverse(DBstr.begin(), DBstr.end());
  cout << Qstr << endl;

  for(int i = 0; i < Qstr.size(); i++) {
    if(Qstr[i] == '-' || DBstr[i] == '-') {
      cout << " ";
    }
    else {
      cout << "|";
    }

  }
  cout << endl;

  cout << DBstr << endl;
}

// -------- BEGIN TRAVERSE PATH ---------------
void traversePath(string querySeq, string dbSeq, int best_i, int best_j, int **V_table, int **T_table) {
  int counter = 0;
  string Qstr = "";
  string DBstr = "";
  int tempI = best_i;
  int tempJ = best_j;

  while(V_table[best_i][best_j] != 0) {
    switch (T_table[best_i][best_j]) {
      case 1:
        // diagonial
        Qstr += querySeq[best_i];
        DBstr += dbSeq[best_j];
        best_i--;
        best_j--;
        break;
      case 2:
        // up to down
        Qstr += querySeq[best_i];
        DBstr += '-';
        best_i--;
        break;
      case 3:
        // left side
        Qstr += '-';
        DBstr += dbSeq[best_j];
        best_j--;
        break;
    }
    counter++;
  }

  cout << "Query:  alignment start index = " << best_i << " end index = " << best_i + (tempI - best_i) << endl;
  cout << "DB seq:  alignment start index = " << best_j << " end index = " << best_j + (tempJ - best_j) << endl;
  showInverted(Qstr, DBstr);
}
// ------- END TRAVERSE PATH ------------

// ------ BEGIN MAIN --------
int main(int argc,char* argv[]) {
  fstream QueryInput;
  fstream databaseInput;
  string label;
  string querySequence;
  string databaseSequence;


  QueryInput.open(argv[1]); // open query

  while(QueryInput >> label) {
    QueryInput >> querySequence;
  }

  QueryInput.close(); // close query

  databaseInput.open(argv[2]); // open database

  while(databaseInput >> label) {
    databaseInput >> databaseSequence;
    cout << label << endl;

    // call S_W(querySeq, DB_seq)
    SmithWaterman(querySequence, databaseSequence);
    cout << endl;
  }

  databaseInput.close(); // close database

  return 0;
}
