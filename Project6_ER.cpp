/*
  Elvin Rivera

  This program exhibits the motif finding probelm through the
  branch and bound algorithm

  -take in fasta file and build unordered_map
  -find median string using a priority_q, keep best 5 median strings
    along with their total distances
  -use the best 5 median strings, find corresponding motif consensus

  *** currently not working, can't compile or run

*/

#include <iostream>
#include <unordered_map>  // for unordered_map
#include <string>
#include <fstream>  // for open, close files

using namespace std;

void succeeding_V(string v_string, &level, Lmer, K);
void bypass(string v_string, &level, Lmer, K);
int totalDistance(string *prefix_string, unordered_map<string,string> umap, int level);


// succeeding_V begin
void succeeding_V(string *v_string, &level, Lmer, K)
{
  if(level != 0)
  {
    level -= 1;
  }

  string vString = v_string[level];

  //change vString
  for(int i = 0; i < vString.size(); i++)
  {
    if(vString[i] == 'A')
    {
      vString[i] = 'G'
      break;
    }
    else if(vString[i] == 'G')
    {
      vString[i] = 'C';
      break;
    }
    else if(vString[i] == 'C')
    {
      vString[i] = 'T';
      break;
    }
  }

  //pass back into array
  v_string[level] = vString;
} //end succeeding_V


//bypass begin
void bypass(string *v_string, &level, Lmer, K) //K is alphabet size 4
{
  // max level not achieved
  if(level != (Lmer - 1))
  {
    level += 1; //go up one level
  }

  string vString = v_string[level];

  //find out how to change v_string
  for(int i = 0; i < vString.size(); i++)
  {
    if(vString[i] == 'A')
    {
      vString[i] = 'G'
      break;
    }
    else if(vString[i] == 'G')
    {
      vString[i] = 'C';
      break;
    }
    else if(vString[i] == 'C')
    {
      vString[i] = 'T';
      break;
    }
  }

  //pass vString back to array
  v_string[level] = vString;
}


//should return the total distance of the prefixstring from the unordered
int totalDistance(string *prefix_string, unordered_map<string,string> umap, int level)
{
  int total = 0;
  int distanceTotal = 0;
  int count = 0;
  string gene;
  string geneLetters;

  for (auto local_it = umap.begin(i); local_it != mymap.end(i); local_it++)
  {
    distanceTotal += total;
    total = level;
    gene = local_it->second;

    for(int i = 0; i < gene.size(); i++)
    {
      count = 0;
      geneLetters = gene.substr(i, level);

      for(int j = 0; j < level; j++)
      {
        if(geneLetters[j] != prefix_string[j])
          count++;
      }

      if(count < total)
        total = count;
    }
  }

  return distanceTotal;
}

// begin main
main(int argc, char* argv[])
{
  fstream input;  // for reading files
  string label;
  string sequence;
  string *v_string;
  string *priority_q;
  int *distnace;
  string prefix_string;
  int level;
  int worst_entry_of_priority_q = 0;

  //declare unordered_map<label, seq>; //or, any appropriate datastructure
  unordered_map<string,string> umap; // initalize unordered_map

  //declare priority_q of size 5;
  priority_q = new string[5];
  distance = new int[5];

  string lmerLength = argv[1];

  //open input fasta file;
  input.open(argv[2]);

  // while not eof read a label & sequence, and place that into the unordered_map;

  while(!input.eof())
  {
    input >> label;
    input >> sequence;
  }
  input.close(); // close input fasta file;

  //initialize v_string[Lmer]
  v_string = new string[lmerLength];

  for(int i = 0; i < lmerLength; i++)
  {
    for(int j = -1; j < i; j++)
    {
      v_string[i] += "A";
    }
  }

  level = stoi(lmerLength);
  while(level > 0)
  {
    if (level < lmerLength) //if non_leaf node
    {
      // generate prefix_string from v_string;
      int optimistic_dist = totalDistance(v_string, umap, level);

      if (optimistic_dist >= worst_entry_of_priority_q)
        bypass(v_string, level, lmerLength, 4); //4 is alphabet size (A,C,G,T)
      else
        succeeding_V(v_string, level, lmerLength, 4);
     }
     else //leaf node
     {
       //display this v_string;
       cout << v_string[level] << endl;

       int totalDistance = totalDistance(v_string, umap, lmerLength);
       //update priority_q if total_distance is better than the worst entry;
       if (totalDistance >= worst_entry_of_priority_q)
       {
         for(int x = 0; x < 5; x++)
         {
           if(priority_q[x].empty())
           {
             priority_q[x] = v_string[level];
             distnace[x] = totalDistance;
           }
           else
           {
             if(distance[x] > totalDistance)
             {
               priority_q[x] = v_string[level];
               distance[x] = totalDistance;
             }
           }
         }
       }

       succeeding_V(v_string, level, lmerLength, 4);
     }
  }//while

  //best 5 median string and their total_distances;
  for(int d = 0; d < 5; d++)
  {
    cout << priority_q[d] << endl;
    cout << distance[d] << endl;
  }

}//main
