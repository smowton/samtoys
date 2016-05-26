
#include <functional>
#include <iostream>

#include <stdlib.h>

const char* ws = " \t\n";

bool first_2_fields(std::string& l) {

  size_t boundary = l.find_first_of(ws);
  if(boundary == std::string::npos)
    return false;

  size_t field2start = l.find_first_not_of(ws, boundary);
  if(field2start == std::string::npos)
    return false;

  size_t field2end = l.find_first_of(ws, field2start);
  if(field2end == std::string::npos)
    return false;

  l.erase(field2end);
  return true;

}

int main(int argc, char** argv) {

  if(argc != 2) {
    std::cerr << "Usage: sample_pileup proportion_to_keep < in.pileup > out.pileup\n";
    exit(1);
  }

  std::string propstr = argv[1];
  double prop = std::stof(propstr);
  size_t hash_threshold = (size_t)(((double)SIZE_MAX) * prop);
  std::cerr << "Keeping records with hash <= " << hash_threshold << "\n";

  std::string thisline;
  std::string thisline_chroffset;
  std::hash<std::string> hashfn;

  while(true) {

    std::getline(std::cin, thisline);
    size_t firstc = thisline.find_first_not_of(ws);
    if(firstc == std::string::npos)
      break;

    thisline_chroffset = thisline;

    if(!first_2_fields(thisline_chroffset)) {
      std::cerr << "Malformed line: " << thisline << "\n";
      exit(1);
    }

    size_t hashval = hashfn(thisline_chroffset);
    if(hashval <= hash_threshold)
      std::cout << thisline << "\n";

  }

  return 0;

}
