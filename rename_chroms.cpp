
#include <htslib/hts.h>
#include <htslib/sam.h>

#include <string.h>

#include <map>
#include <string>
#include <iostream>

int main(int argc, char** argv) {

  std::map<std::string, std::string> new_chroms;
  for(int i = 1; i < 23; ++i) {

    char chrbuf[32];
    sprintf(chrbuf, "chr%d", i);
    new_chroms[std::string(chrbuf)] = std::string(chrbuf+3);

  }
    
  new_chroms["chrX"] = "X";
  new_chroms["chrY"] = "Y";
  new_chroms["chrM_rCRS"] = "MT";

  htsFile* hf = hts_open(argv[1], "r");
  if(!hf) {
    fprintf(stderr, "Failed to open %s\n", argv[1]);
    exit(1);
  }

  bam_hdr_t* header = sam_hdr_read(hf);

  htsFile* hfo = hts_open("-", "w");
  if(!hfo) {
    fprintf(stderr, "Failed to open stdout\n");
    exit(1);
  }

  for(int32_t i = 0; i < header->n_targets; ++i) {

    std::map<std::string, std::string>::iterator it = new_chroms.find(header->target_name[i]);

    if(it == new_chroms.end())
      continue;

    free(header->target_name[i]);
    header->target_name[i] = strdup(it->second.c_str());
	 
  }

  // Void cache if any (leaks a khash table; never mind)
  if(header->sdict)
    header->sdict = 0;

  // Remove any @SQ lines from the header text so that our rewrite will take effect:
  std::string htext = header->text;
  size_t sqstart;
  while((sqstart = htext.find("@SQ")) != std::string::npos) {
    size_t nextnl = htext.find("\n", sqstart);
    if(nextnl == std::string::npos) {
      std::cerr << "Malformed header SQ line\n";
      exit(1);
    }
    htext.erase(sqstart, (nextnl - sqstart) + 1);
  }

  free(header->text);
  header->text = strdup(htext.c_str());
  header->l_text = htext.length();
    
  sam_hdr_write(hfo, header);

  hts_close(hf);
  hts_close(hfo);

}
