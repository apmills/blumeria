// C++

#include <iostream>
#include <fstream>
using namespace std;

int main() {
  char character;
  // string tring;
  ifstream instream;
  // ofstream outstream;
  instream.open("printo.cc");
  while (!instream.eof()) {
    instream.get(character);
    // instream >> tring;
    cout << character;
  }
  instream.close();
  return 0;
}
