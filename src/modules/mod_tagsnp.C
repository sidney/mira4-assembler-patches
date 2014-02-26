/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 1997-2000 by the German Cancer Research Center (Deutsches
 *   Krebsforschungszentrum, DKFZ Heidelberg) and Bastien Chevreux
 * Copyright (C) 2000 and later by Bastien Chevreux
 *
 * All rights reserved.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the
 * Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
 *
 */

#include <iostream>

#include "caf/caf.H"
#include "mira/assembly.H"

#include "modules/misc.H"
#include "modules/mod_tagsnp.H"

#include "version.H"


using namespace std;



list<Contig> tagsnp::TS_clist;
ofstream tagsnp::TS_fout;

void tagsnp::usage()
{
  cerr << "tagsnp\t(MIRALIB version " << MIRALIBVERSION << ")\n\n";
  cerr << "Usage:\n";
  cerr << "  tagsnp [-xxx] cafin cafout [optional MIRA settings]\n\n";
  cerr << "Options:\n";
  cerr << "\t-n\t\tnuke all existing SNP and RMB tags in file\n";
//  cerr << "\t-s <filename>\tload strain data from file\n";
//  cerr << "\t-a\t\tassume SNPs instead of repeats\n";
//  cerr << "\t-r <int>\tminimum reads per group (default: 1)\n";
//  cerr << "\t-q <int>\tminimum qual for tagging (default: 30)\n";
//  cerr << "\t-n <int>\tminimum neighbour qual for tagging (default: 20)\n";
//  cerr << "\t-e <int>\tend read exclusion area (default: 25)\n";
//  cerr << "\t-g\t\talso mark gap bases\n";
//  cerr << "\t-m\t\talso mark multicolumn gap bases\n";

}


//void tagsnp::save(string & cafout)
//{
//
//  if(!cafout.empty()){
//    assout::saveAsCAF(contigs,cafout);
//  } else {
//    assout::dumpAsCAF(contigs,cout);
//  }
//
//  //
//  //filename="out.tcs";
//  //assout::saveAsTCS(contigs,filename);
//
//  //filename="tagsnp_out.gap4da";
//  //assout::saveAsGAP4DA(contigs,filename);
//
//  //filename="featureanalysis.txt";
//  //assout::saveFeatureAnalysis(400,100,contigs,readpool,
//  //				filename,
//  //				"featuresummary.txt",
//  //				"featureprot.txt");
//
//  //{
//  //  string filename="out.html";
//  //  cout << "Saving contigs to htmlfile: " << filename << endl;
//  //  ofstream out(filename.c_str(), ios::out | ios::trunc);
//  //  assout::dumpContigListAsHTML(contigs, "Super project", out);
//  //  out.close();
//  //}
//}

//void tagsnp::load (MIRAParameters * mp, string & cafin, string & strainin)
//{
//  cerr << "Loading project from CAF file: " << cafin << endl;
//
//  CAF tcaf(readpool, contigs, mp);
//  tcaf.load(cafin);
//
//  if(!strainin.empty()){
//    cerr << "Loading strain data";
//    readpool.loadStrainData(strainin);
//  }
//
//  Assembly::refreshContigAndReadpoolValuesAfterLoading(readpool,contigs);
//}

//void tagsnp::doit(list<Contig> & contigs)
//{
//  cout << "Tagging reads ..." << endl;
//  list<Contig>::iterator I=contigs.begin();
//  for(;I!=contigs.end(); I++){
//    //I->setParams(&P);
//    //
//    //uint32 numSRMB=0;
//    //uint32 numWRMB=0;
//    //uint32 numSNP=0;
//    //I->transposeReadSRMTagsToContig();
//    ////I->markPossibleRepeats(numSRMB, numWRMB, numSNP);
//    //vector<bool> readsmarkedsrm;
//    //I->newMarkPossibleRepeats(numSRMB,readsmarkedsrm);
//  }
//}


void tagsnp::saveCList(list<Contig> & clist, ReadPool & rp)
{
  Contig::setCoutType(Contig::AS_CAF);
  list<Contig>::iterator I=clist.begin();
  for(;I!=clist.end(); I++){
    TS_fout << *I;
  }
}

void tagsnp::cafload_callback(list<Contig> & clist, ReadPool & rp)
{
  bool dooutput=true;

  Assembly::refreshContigAndReadpoolValuesAfterLoading(rp,clist);
  clist.back().trashConsensusCache(false);

  doit2(clist);
  saveCList(clist, rp);

  clist.clear();
  rp.discard();
}

void tagsnp::doit2(list<Contig> & contigs)
{
  cout << "Tagging reads ..." << endl;
  list<Contig>::iterator I=contigs.begin();
  for(;I!=contigs.end(); I++){
    I->trashConsensusCache(false);

    Contig::repeatmarker_stats_t repstats;
    vector<bool> readsmarkedsrm;
    I->newMarkPossibleRepeats(repstats, readsmarkedsrm);

    I->markFeaturesByConsensus(true,true,true);
  }
}



int tagsnp::mainTagSNP(int argc, char ** argv)
{
  FUNCSTART("int mainTagSNP(int argc, char ** argv)");

  vector<MIRAParameters> Pv;
  MIRAParameters::setupStdMIRAParameters(Pv);

  const_cast<contig_parameters &>(Pv[0].getContigParams()).con_disregard_spurious_rmb_mismatches=false;

  ReadPool thepool(&Pv);

  int c;
  extern char *optarg;
  extern int optind;


  string cafin="";
  string strainin="";

  while (1){
    c = getopt(argc, argv, "+h");
    if(c == -1) break;

    switch (c) {
    case 'h':
    case '?': {
      usage();
      exit(0);
    }
    default : {}
    }
  }

  if(argc-optind < 2) {
    cerr << argv[0] << ": " << "Missing at least infile or outfile as argument!\n";
    usage();
    exit(1);
  }

  string infile=argv[optind++];
  string outfile=argv[optind++];

  if(argc-optind > 0) {
    stringstream tss;
    for(int32 i=optind; i<argc; i++) tss << argv[i] << "  *=BEGIN0=*";
    MIRAParameters::parse(tss,Pv,false);
  }

  MIRAParameters::dumpAllParams(Pv, cout);

  TS_fout.open(outfile.c_str(), ios::out);
  if(!TS_fout){
    MIRANOTIFY(Notify::FATAL, "Could not open file for saving: " << outfile);
  }

  try{
    vector<uint32> dummy;
    CAF tcaf(&thepool, &TS_clist, &Pv);
    tcaf.load(infile,
	      ReadGroupLib::SEQTYPE_SANGER,
	      1,
	      dummy,
	      false,
	      cafload_callback
      );

    //load(&P, infile, strainin);
    //
    //doit();
    //save(outfile);

  }
  catch(Notify n){
    n.handleError("main");
  }
  catch(Flow f){
    cerr << "Unexpected exception: Flow()\n";
  }

  TS_fout.close();

  FUNCEND();
  return 0;
}



