
#include <htslib/hts.h>
#include <htslib/sam.h>

#include <string.h>

#include <map>
#include <string>
#include <iostream>
#include <algorithm>

int main(int argc, char** argv) {

  std::vector<std::string> newOrder;
  for(int i = 1; i < 23; ++i) {

    char chrbuf[32];
    sprintf(chrbuf, "%d", i);
    newOrder.push_back(std::string(chrbuf));

  }
    
  newOrder.push_back("X");
  newOrder.push_back("Y");
  newOrder.push_back("MT");

  htsFile* hf = hts_open(argv[1], "r");
  if(!hf) {
    fprintf(stderr, "Failed to open %s\n", argv[1]);
    exit(1);
  }

  bam_hdr_t* header = sam_hdr_read(hf);

  std::vector<std::string> oldOrder;

  for(int i = 0; i < header->n_targets; ++i)
    oldOrder.push_back(header->target_name[i]);

  std::map<int32_t, int32_t> oldToNew;

  for(int32_t i = 0; i < oldOrder.size(); ++i) {

    std::vector<std::string>::iterator findit = std::find(newOrder.begin(), newOrder.end(), oldOrder[i]);
    if(findit == newOrder.end()) {
      if(i < newOrder.size()) {
	std::cerr << oldOrder[i] << " clashes with new order\n";
	exit(1);
      }
      oldToNew[i] = i;
    }
    else {
      int32_t newIdx = std::distance(newOrder.begin(), findit);
      oldToNew[i] = std::distance(newOrder.begin(), findit);
    }

  }

  htsFile* hfo = hts_open("-", "wb0");
  if(!hfo) {
    fprintf(stderr, "Failed to open stdout\n");
    exit(1);
  }
  
  // Rearrange the header:
  bam_hdr_t* newheader = bam_hdr_dup(header);

  for(int32_t i = 0; i < header->n_targets; ++i) {

    int32_t newIdx = oldToNew[i];
    free(newheader->target_name[newIdx]);
    newheader->target_name[newIdx] = strdup(header->target_name[i]);
    newheader->target_len[newIdx] = header->target_len[i];

  }

  // Void cache if any (leaks a khash table; never mind)
  if(newheader->sdict)
    newheader->sdict = 0;

  // Remove any @SQ lines from the header text so that our rewrite will take effect:
  std::string htext = newheader->text;
  size_t sqstart;
  while((sqstart = htext.find("@SQ")) != std::string::npos) {
    size_t nextnl = htext.find("\n", sqstart);
    if(nextnl == std::string::npos) {
      std::cerr << "Malformed header SQ line\n";
      exit(1);
    }
    htext.erase(sqstart, (nextnl - sqstart) + 1);
  }

  // Replace @SQ lines, for legacy readers:
  if(htext.size() > 0 && htext[htext.size()-1] != '\n')
    htext += "\n";

  for(int32_t i = 0; i < newheader->n_targets; ++i) {
    char buf[128];
    if(snprintf(buf, 128, "@SQ\tSN:%s\tLN:%u\n", newheader->target_name[i], newheader->target_len[i]) > 128) {
      fprintf(stderr, "Sequence name too long!\n");
      exit(1);
    }
    htext += std::string(buf);
  }

  free(newheader->text);
  newheader->text = strdup(htext.c_str());
  newheader->l_text = htext.length();

  sam_hdr_write(hfo, newheader);

  bam1_t *rec = bam_init1();

  // Leave unmapped reads alone
  oldToNew[-1] = -1;

  while(sam_read1(hf, header, rec) >= 0) {

    if(!oldToNew.count(rec->core.tid)) {
      std::cerr << "Unknown contig ID " << rec->core.tid << "!\n";
      exit(1);
    }
    if(!oldToNew.count(rec->core.mtid)) {
      std::cerr << "Unknown contig ID " << rec->core.mtid << "!\n";
      exit(1);
    }

    rec->core.tid = oldToNew[rec->core.tid];
    rec->core.mtid = oldToNew[rec->core.mtid];
    sam_write1(hfo, newheader, rec);

  }

  hts_close(hf);
  hts_close(hfo);

}
