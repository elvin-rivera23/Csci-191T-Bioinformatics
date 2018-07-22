/*

Elvin Rivera
CSCI 191T

Implement practical algorithm for the "Partial Digest Problem" To locate DNA Restriction sites
	- use recursion method to yield all possible sets of the restriction sites from a sets
	of pairwise distances

-Takes in text file containing a set of distances
	- Expected input: L = { 2, 3, 3, 3, 4, 5, 6, 6, 7, 7, 8, 9, 10, 10, 11, 12, 14, 14, 17, 17, 20 }

-Output: Show all restriction sites, can have more than one list, sort the lists in ascending order

-to compile the program in cygwin: g++ "file"

*/
#include <iostream>	// input output
#include <algorithm>
#include <math.h> //mathematical functions
#include <stdlib.h>	// for variable types
#include <string>
#include <fstream> //input output file stream class
#include <vector>

using namespace std;
fstream fin;
int a;
//function prototype
void Modify(vector<int>, vector<int>);

//----------MODIFY FUCNTION BEGIN----------------
void Modify(vector<int> primaryList, vector<int> secondaryList)
{
	//Vector Lists
	vector<int>distance;
	vector<int>temp;

	//position in list
	int position = 0;

	bool found = true;
	int b;

	//run through to see if list is empty
	if (primaryList.empty())
	{
		sort(secondaryList.begin(), secondaryList.end());
		cout << "X{";

		for (int i = 0; i < secondaryList.size(); i++)
		{
			cout << secondaryList[i] << ", ";
		}

		cout << "} " << endl;
		return;
	}

	//sort the list
	sort(primaryList.begin(), primaryList.end());

	//make c the last element of list
	int c = primaryList.back();

	//creating the distance
	for (int i = 0; i < secondaryList.size(); i++)
	{
		int temp = (int)abs(c - secondaryList[i]);
		distance.push_back(temp);
	}

	//resizing temp list to match that of list 1
	temp.resize(primaryList.size());

	for (int i = 0; i< primaryList.size(); i++)
		temp[i] = primaryList[i];

	found = 1;

	//determine the distance by using variable j
	for (int j = 0; j < distance.size(); j++)
	{

		//creating iterator
		vector <int>::iterator iter = temp.begin();

		iter = find(temp.begin(), temp.end(), distance[j]);
		position = (int)distance(temp.begin(), iter);

		//Ensure lists match sizes
		if (position == temp.size())
			position -= 1;

		//checking if position is found or not
		if (temp[position] != distance[j])
		{
			found = false;
			break;
		}

		else
		{
			temp.erase(temp.begin() + position);
		}
	}

	//If found is true continue on
	if (found) {

		secondaryList.push_back(c);

		//resize primaryList to new size and copy
		primaryList.resize(temp.size());
		for (int i = 0; i < temp.size(); i++)
		{
			primaryList[i] = temp[i];
		}
		Modify(primaryList, secondaryList);


		for (int i = 0; i < secondaryList.size(); i++)
		{

			if (c == secondaryList[i])
			{
				position = i;
				break;
			}

		}

		//erase contents of secondaryList
		secondaryList.erase(secondaryList.begin() + position);

		//Restoring primaryList
		for (int i = 0; i < distance.size(); i++)
			primaryList.push_back(distance[i]);

	}


	//reinitilize both lists to 0
	temp.resize(0);
	distance.resize(0);

	b = a - c;

	//using list to to make distance vector
	for (int i = 0; i < secondaryList.size(); i++)
	{
		int distance2 = (int)abs(b - secondaryList[i]);
		distance.push_back(distance2);
	}

	//setting temp to size of primaryList
	temp.resize(primaryList.size());

	for (int i = 0; i < primaryList.size(); i++)
	{
		temp[i] = primaryList[i];
	}

	found = 1;

	for (int i = 0; i < distance.size(); i++)
	{
		//implemnt iteration method to go through
		vector <int>::iterator iter = temp.begin();
		iter = find(temp.begin(), temp.end(), distance[i]);
		position = (int)distance(temp.begin(), iter);

		//ensure no size mismatch
		if (position == temp.size())
			position -= 1;

		//false case if position does not exist
		if (temp[position] != distance[i])
		{
			found = false;
			break;
		}

		//otherwise
		else
			temp.erase(temp.begin() + position);
	}

	//if found continues to be true
	if (found) {

		secondaryList.push_back(b);

		////setting temp to size of primaryList
		primaryList.resize(temp.size());

		for (int i = 0; i < primaryList.size(); i++)
		{
			primaryList[i] = temp[i];
		}

		Modify(primaryList, secondaryList);


		//Checks for position
		for (int i = 0; i < secondaryList.size(); i++)

		{

			if (b == secondaryList[i])

			{
				position = i;
				break;
			}

		}
		//erase if needed
		secondaryList.erase(secondaryList.begin() + position);

		//Reconstruct
		for (int j = 0; j < distance.size(); j++)
			primaryList.push_back(distance[j]);
	}

	return;
}
//-----------MODIFY FUNCTION END -------------



//--------MAIN BEGIN--------------------
//follow sample_vector program
int main(int argc, const char * argv[])
{

	//L on output
	vector<int>primaryList;
	//X on output
	vector<int>secondaryList;

	fin.open("distances.txt");

	//dealing with any spaces if present
	while (!(fin.eof()))
	{
		string str;
		fin >> str;
		primaryList.push_back((atoi(str.c_str())));
	}

	fin.close();

	//Begin to display in appropriate format
	cout << "L = { ";
	for (int i = 0; i<primaryList.size(); i++)
	{
		if (i == primaryList.size() - 1)
		{
			cout << primaryList[i] << " }" << endl;
		}

		else
			cout << primaryList[i] << ", ";
	}

	//sort through list
	sort(primaryList.begin(), primaryList.end());
	a = primaryList.back();
	primaryList.pop_back();
	secondaryList.push_back(0);
	secondaryList.push_back(a);

	//recursive implementation
	Modify(primaryList, secondaryList);

	system("pause");
	return 0;
}
//--------------MAIN END----------
