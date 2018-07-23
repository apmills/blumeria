// C++

#include <iostream>
using namespace std;

int main() {
  int year_now, age_now, another_year, another_age;
  cout << "Enter the current year then press Enter. \n";
  cin >> year_now;

  cout << "Enter your current age in years.\n";
  cin >> age_now;

  cout << "Enter the year for which you wish to know your age.\n";
  cin >> another_year;

  another_age = another_year - (year_now - age_now);
  int difference;
  float lifes;
  difference = another_year - year_now;

  lifes = static_cast<double>(difference) / 85;

  if (another_age >= 0) {
    cout << "Your age in " << another_year <<": " << another_age << "\n";
    cout << "That's " << lifes << " lifetimes away\n";
  } else {
    cout << "You weren't even born in " << another_year << "!\n";
  }
  return 0;
}
