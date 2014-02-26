#include <typeinfo>
#include <iostream>
#include <string>

#include <sys/times.h>
#include <limits.h>
#include <unistd.h>

#include "stdinc/types.H"
#include "stdinc/defines.H"

#include "io/fastq-mira.H"
#include "io/fastq-lh.H"

#include "util/misc.H"
#include "util/machineinfo.H"

#include "errorhandling/errorhandling.H"


#include <boost/filesystem.hpp>

KSEQ_INIT(gzFile, gzread)

using namespace std;


#include "mira/readpool.H"
#include "mira/skim.H"
#include "mira/bloomfilter.H"
#include "mira/seqtohash.H"
#include "util/dptools.H"
#include "mira/hashstats.H"
void ttt()
{
  FUNCSTART("ttt");
  auto rgid = ReadGroupLib::getReadGroupID(0);

  vector<MIRAParameters> Pv;
  MIRAParameters::setupStdMIRAParameters(Pv);
  ReadPool rp(&Pv);

  cout << "Have " << MachineInfo::getMemAvail() << " mem avail" << endl;
  rp.loadData_rgid(
    "fastq",
    "bla.fastq",
    "",
    rgid,
    false,
    nullptr);
  cout << "Have " << MachineInfo::getMemAvail() << " mem avail" << endl;

  NHashStatistics nhs;
  nhs.setupNewAnalysis(32,4,31,2);
  //nhs.setupNewAnalysis(26,4,31);
  nhs.analyseReadPool(rp);

/*
  nhs.setupNewAnalysis(24,4,31);

 Trimmed 736566 hashes, 3543078 remaining
Localtime: Thu Dec 13 15:04:27 2012
BloomFilter:
memory: 4194304
keys per kmer: 4
l1 occupancy: 12738052
l2 occupancy: 9946103
probable num kmers   : 5898406
 thereof num kmers>=2: 4279644
 thereof num kmers>=3: 3957220

FP=100/(3543078+736566)*736566 = 17.2%


  nhs.setupNewAnalysis(25,4,31);

 done. Trimmed 244140 hashes, 3791601 remaining
Localtime: Thu Dec 13 15:08:27 2012
BloomFilter:
memory: 8388608
keys per kmer: 4
l1 occupancy: 18357631
l2 occupancy: 12031047
probable num kmers   : 7046122
 thereof num kmers>=2: 4035741
 thereof num kmers>=3: 3053407

FP=100/(3791601+244140)*244140 = 6.0%


  nhs.setupNewAnalysis(26,4,31);

 done. Trimmed 39949 hashes, 3820645 remaining
Localtime: Thu Dec 13 15:12:37 2012
BloomFilter:
memory: 16777216
keys per kmer: 4
l1 occupancy: 22986230
l2 occupancy: 13388583
probable num kmers   : 7386566
 thereof num kmers>=2: 3860594
 thereof num kmers>=3: 2888110

FP=100/(3820645+39949)*39949 = 1.0%


*/

  nhs.dumpHealth(cout);
  cout << "hash distrib:\n";
  nhs.dumpHashDistrib(cout);
  //nhs.sortLexicographically();
  //cout << "hash count:\n";
  //nhs.dumpHashCount(cout);
  cout << "trim by 2 2 -1\n";
  nhs.trimHashStatsByFrequency(2,2,-1);
  //cout << "hash count:\n";
  //nhs.dumpHashCount(cout);
  cout << "hash distrib:\n";
  nhs.dumpHashDistrib(cout);

  cout << "trim by 3 3 -1\n";
  nhs.trimHashStatsByFrequency(3,3,-1);
  cout << "hash distrib:\n";
  nhs.dumpHashDistrib(cout);

  cout << "trim by 4 4 -1\n";
  nhs.trimHashStatsByFrequency(4,4,-1);
  cout << "hash distrib:\n";
  nhs.dumpHashDistrib(cout);


  exit(0);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
int main(int argc, char ** argv)
{
  FUNCSTART("int main(int argc, char ** argv)");

  try{
    ttt();
  }
  catch(Notify n){
    n.handleError("main");
  }

  exit(0);

  cout << "Have " << MachineInfo::getCoresTotal() << " cores" << endl;
  cout << "Have " << MachineInfo::getMemTotal() << " mem total" << endl;
  cout << "Have " << MachineInfo::getMemAvail() << " mem avail" << endl;

  srand(1234567);

  string bla("bla");
  boost::system::error_code ec;
  try{
    boost::filesystem::remove_all(bla,ec);
  }
  catch(boost::filesystem::filesystem_error fse){
    cout << fse.what() << endl;
  }

  exit(0);

  try {
    timeval tv;
    suseconds_t sus=0;

//    gzFile fp;
//    kseq_t *seq;
//    fp = gzopen("bla.fastq", "r");
//    seq = kseq_init(fp);
//    gettimeofday(&tv,nullptr);
//    int l;
//    while ((l = kseq_read(seq)) >= 0);
//    sus=diffsuseconds(tv);
//    cout << "timing fastq-lh: " << sus << endl;


    FastQ fq;

    gettimeofday(&tv,nullptr);
    fq.openFile("bla.fastq");
    while(fq.loadNext()>=0){
      //cout << fq.getLineCount() << endl;
      //if(fq.getLineCount()%10000==0) cout << fq.getLineCount() << endl;
    }
    cout << fq.getLineCount() << endl;
    sus=diffsuseconds(tv);
    cout << "timing fastq-mira: " << sus << endl;



  }
  catch(Notify n){
    n.handleError("main");
  }

  return 0;
}
