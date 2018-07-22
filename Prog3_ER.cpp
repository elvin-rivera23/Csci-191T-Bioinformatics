/*
Elvin Rivera
CSCI 191T MW
Prog 3

  -This program locates the median string and consensus string from DNA Sequences
    and the length of consensus string.
  - report 5 median strings and their motif consensus strings and Positions
  - use HMP-part.fa

  ~ start by finding median string, find motif consensus string based on median string

  =to compile program in cygwin type: g++ Prog3_ER.cpp -o Prog3_ER
  =to run program in cygwin type: ./Prog3_ER
*/


#include <iostream> // input/output
#include <fstream>  // file open/close
#include <sstream> // for stringstream
#include <string> // work with strings
#include <vector> // vector constructs
#include <algorithm>


using namespace std;

// function prototypes
void generateCombinations(vector<string> & vs);
int obtainTotDist(vector<string> genes, string letters);
vector<int> startPos; // used in motifString
vector<string> motifString(vector<string> genes, string medianStrings);
string ConsensusMotif(vector<string> motifs);

//-------------FUNCTIONS--------------

// ------- generate combos function begin------
//combinations of ATCG inside Vector
void generateCombinations(vector<string> & vs) {
  string letters[4] = {"A", "C", "G", "T"};

  for(int i0 = 0; i0 < 4; i0++)
    for(int i1 = 0; i1 < 4; i1++)
      for(int i2 = 0; i2 < 4; i2++)
        for(int i3 = 0; i3 < 4; i3++)
          for(int i4 = 0; i4 < 4; i4++)
            for(int i5 = 0; i5 < 4; i5++)
              vs.push_back(letters[i0] + letters[i1] + letters[i2] + letters[i3] + letters[i4] + letters[i5]);
}
// ------- generate combos function end------

// ----- obtain distances function begin ---
// Obtain total distance Function
int obtainTotDist(vector<string> genes, string letters) {

  string gene; // for getting gene i from vector genes
  string gene_Letter;  // for substr(j, 6) gene within gene

  int count = 0;
  int total = 0;
  int overallDist = 0;

  // by end of this for loop we should have the total distance of
  // gene letters passed in
  // for loop 1
  for(int i = 0; i < genes.size(); i++) {
    overallDist += total; // add to total distance
    gene = genes[i]; // get i gene from vector genes
    total = 6; // check Max number of overall distance; reset total to 6
    // locate hamming distance
    //second for loop
    for(int j = 0; j < gene.size() - 5; j++) {
      count = 0;  // reset count to 0
      // get 6 letter from gene
      gene_Letter = gene.substr(j, 6);
      // compare each letter to obtain dist b/w two genes
      // iterate if char doesn't match count
      for(int k = 0; k < 6; k++) { // for3
        if(gene_Letter[k] != letters[k]) {
          count++;
        }
      }// end for3
      // update string if optimal distance found
      if(count < total) {
        total = count;
      }
    }// end for2
  }// end for1
  return overallDist;
}
// ----- obtain distances function end-----


// -----------motif string function begin------------
// vector of motifs string that resemble median string most
vector<string> motifString(vector<string> genes, string medianStrings) {
  vector<string> consensusStrings;
  string gene;
  string gene_Letter;
  string temp_Gene;
  int count = 0;
  int minimum_Count = 6;
  int temp_Pos;
  vector<string> motif;
  // for 1
  for(int i = 0; i < genes.size() - 1; i++) {
    gene = genes[i];
    minimum_Count = 6;
    // for 2
    for(int j = 0; j < gene.size() - 5; j++) {
      gene_Letter = gene.substr(j, 6); // get (j, 6) string from a string of list
      count = 0;
      // for 3
      // Compare chars and count the difference
      for(int k = 0; k < gene_Letter.size(); k++) {
        if(gene_Letter[k] != medianStrings[k]) {
          count++;
        }
      } // end for3
      // better motif found, replace count and obtain position of gene
      if(minimum_Count > count) {
        minimum_Count = count;
        temp_Pos = j;
        temp_Gene = gene_Letter;
      }
    } // end for2
    startPos.push_back(temp_Pos);  // push back the starting position and motif
    motif.push_back(temp_Gene);
  } // end for1
  return motif;
}
// -----------motif string function end----

// ---------consensus motif function begin-----
//Function for finding Consensus Motif
string ConsensusMotif(vector<string> motifs) {
  string obtainMotif;
  string returnConsensusMotif;
  int count_A[6] = {0};
  int count_T[6] = {0};
  int count_C[6] = {0};
  int count_G[6] = {0};
  int countConsensus = 0;
  stringstream ss;  // for converting int to string

  for(int i = 0; i < motifs.size(); i++) {  //for1
    obtainMotif = motifs[i];
    // locate how many times ATCG appears in motif from gene
    for(int j = 0; j < obtainMotif.size(); j++) {  //for2
      if(obtainMotif[j] == 'A')
        count_A[j] += 1;
      else if(obtainMotif[j] == 'T')
        count_T[j] += 1;
      else if(obtainMotif[j] == 'C')
        count_C[j] += 1;
      else
        count_G[j] += 1;
    } //end for2
  }// end for1
  // Finds which char appears more and track countConsensus score
  for(int x = 0; x < obtainMotif.size(); x++) {
    if(count_A[x] > count_T[x] && count_A[x] > count_C[x] && count_A[x] > count_G[x]) {
      returnConsensusMotif += "A";
      countConsensus += count_A[x];
    }
    else if(count_T[x] > count_A[x] && count_T[x] > count_C[x] && count_T[x] > count_G[x]) {
      returnConsensusMotif += "T";
      countConsensus += count_T[x];
    }
    else if(count_C[x] > count_A[x] && count_C[x] > count_T[x] && count_C[x] > count_G[x]) {
      returnConsensusMotif += "C";
      countConsensus += count_C[x];
    }
    else {
      returnConsensusMotif += "G";
      countConsensus += count_G[x];
    }
  }
  ss << " " << countConsensus << "/";   // get consensus score
  ss << (obtainMotif.size() * motifs.size());   // Get max total disance and convert from int to string
  returnConsensusMotif += ss.str(); // pass in ss string to returnConsensusMotif
  return returnConsensusMotif;
}
// ---------consensus motif function end-----


//----Main Begin--------------------

int main() {

  vector<int> gene_total_Dist;
  vector<int> medianStringDistance;
  vector<string> genes;
  vector<string> vector_Letters;
  vector<string> track_Gene;
  vector<string> medianStringName;
  vector<string> motifStrings;
  vector<string> tempVector;

  fstream inputFile;
  inputFile.open("HMP-part.fa");
  string label1, label2, gene;
  string ConsensusMotif;
  int totalDistance;
  // pass file input into genes vector

  cout << "Before file" << endl;
  while(!inputFile.eof()) {
    inputFile >> label1 >> label2 >> gene;
    genes.push_back(gene);
  }
  inputFile.close();  // close file
  cout << "after file" << endl;

  // reference vector_Letters
  // will create all possible combinations of ATCG
  generateCombinations(vector_Letters);

  // get the total distances and genes
  for(int i = 0; i < vector_Letters.size(); i++) {
    totalDistance = otainTotDist(genes, vector_Letters[i]);
    gene_total_Dist.push_back(totalDistance);
    track_Gene.push_back(vector_Letters[i]);
  }
  // find the best 5 distances out of the list of distances
  int find = 0;
  int count = 0;
  while(count != 5) {
    for(int j = 0; j < gene_total_Dist.size(); j++) {
      if(find == gene_total_Dist[j] && count != 5) {
        medianStringDistance.push_back(find);
        medianStringName.push_back(track_Gene[j]);
        count++;
      }
    }
    find++;
  }

  // find motifs
  for(int i = 0; i < medianStringName.size(); i++) {
    motifStrings = motifString(genes, medianStringName[i]);
    ConsensusMotif = ConsensusMotif(motifStrings);

    //show median string and total distance
    cout << "Median String: "<< medianStringName[i] << " (Total Distance:" <<
    medianStringDistance[i] << ")" << endl;

    // Display this iteration's consensus motif
    cout << "motif consensus string: "<< ConsensusMotif << endl;

    // Display starting Positions and motif Strings
    // Clear starting pos for next iteration
    for(int j = 0; j < motifStrings.size(); j++) {
      cout << "(" << startingPosition[j] << ") "<< motifStrings[j] << ", ";
      startingPosition.clear();
    }
    cout << "\n\n";
  }
  return 0;
}
