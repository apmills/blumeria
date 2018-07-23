// C++

#include <iostream>
#include <fstream>
using namespace std;

int main() {
  /
  cout << "Testing: " << 16/2 << " = " << 4*2 << ".\n\n";
  /
  char ch;
  bool in_com, potential;
  in_com = false;
  potential = false;
  ifstream instream;
  ofstream outstream;
  instream.open("comment.cc");
  outstream.open("WithoutComments.cc");

  /
  while (!instream.eof()) {
    if (in_com) {
      instream.get(ch);
      /
      if (static_cast<int>(ch) == 42) {
        potential = true;
      }
      if (static_cast<int>(ch) == 47 && potential) {
        in_com = false;
      }
    } else {
      instream.get(ch);
      if (static_cast<int>(ch) == 47){
        potential = true;
      } else if (potential && static_cast<int>(ch) == 42){
        in_com = true;
        potential = false;
        continue;
      } else {
        potential = false;
      }
      cout << ch;
      outstream.put(ch);
    }
  }
  return 0;
}

