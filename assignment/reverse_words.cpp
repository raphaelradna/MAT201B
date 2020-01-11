//Raphael Radna
//MAT-201B W20
//Assignment 0

#include <iostream>
#include <string>

using namespace std;

int main() {
  while (true) {
    printf("Enter a sentence to reverse word-wise: ");
    string input;
     getline(cin, input);
    if (!cin.good()) {
      printf("Done\n");
      return 0;
    }
    string output, word;
    for (int i = 0; i < input.length(); i++) {
      if (input.at(i) != ' ') {
        word.insert(0, input, i, 1);
      } else {
        output.append(word + " ");
        word.clear();
      }
    }
    cout << output << word << endl;
  }
}
