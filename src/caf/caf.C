/*
 * Written by Thomas Pfisterer
 *
 * Copyright (C) 1997-2000 by the German Cancer Research Center (Deutsches
 *   Krebsforschungszentrum, DKFZ Heidelberg) and Thomas Pfisterer
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
 *
 * Reading CAF-Files  Version 1.6   19.11.1999
 *
 * Written by Thomas Pfisterer 13.08.1998
 *
 * Load a CAF-File
 * Create a CAF-Object. This Object can load CAF-Files via the load
 * method. Reads and Contigs from the CAF-file are stored in the readpool
 * and the contig-list.
 *
 *
 * 10.05.99  Some missing keywords added (stolen, asped, dye, processStatus,
 *           primer)
 * 22.06.99  Error with "Primer" keyword fixed
 * 27.07.99  Flex is now using environments; handling for comments improved
 *           (problem with '\n' fixed) 0.91
 * 13.08.99  Using environments in flex for most keword constructions (0.92)
 * 15.09.99  Scanner uses \0 to detect non-caf files
 * 13.10.99  Komma and underscores allowed in identifiers.
 *           "Primer" keyword simplified (Universal/Custom/Unknown dropped)
 * 02.05.03  Small changes by BaCh
 * later     See comments in CVS
 */


#include "caf/caf_tokens.h"
#include "caf/caf.H"

#include "io/annotationmappings.H"
#include "util/progressindic.H"



using namespace std;


void showIntervall(caf_intervall i, ostream &output)
{
  output << "[" << i.leftBorder << "; " << i.rightBorder << "] ";
}


void showIntervallList(list<caf_intervall> l, ostream &output)
{
  list<caf_intervall>::iterator I = l.begin();
  cout << "Intervallliste: " << endl;
  while (I != l.end()) {
    showIntervall(*I, output);
    output << endl;
    I++;
  }
  output << endl;
}


caf_intervall schnittmenge(caf_intervall a, caf_intervall b)
{
  caf_intervall schnittmenge;

  // works also if one of the borders is open (i.e. -1)
  schnittmenge.leftBorder  = max(a.leftBorder, b.leftBorder);

  // if one of the borders is open, the other is the new border and the max
  if (a.rightBorder >= 0 &&  b.rightBorder >= 0) {
    schnittmenge.rightBorder = min(a.rightBorder, b.rightBorder);
  } else {
    if (a.rightBorder < 0) {
      schnittmenge.rightBorder = b.rightBorder;
    } else {
      schnittmenge.rightBorder = a.rightBorder;
    }
  }


  // showIntervall(a, cout); showIntervall(b, cout);
  // cout << " => ";
  // showIntervall(schnittmenge, cout);
  // cout << endl;

  return schnittmenge;
}



bool intervall_leer(caf_intervall a)
{
  // if one of the borders is open, the intervall is not empty
  if (a.leftBorder < 0 || a.rightBorder < 0) {
    return false;
  }
  return (a.leftBorder > a.rightBorder);
}



// ============================================================



CAF::CAF(ReadPool * aPool, list<Contig> * aContiglist, vector<MIRAParameters> * mp) {
  thePool    = aPool;
  theContigs = aContiglist;
  CAF_miraparams = mp;

  CAF_isPadded = -1;
  CAF_insert_size_min = -1;
  CAF_insert_size_max = -1;
  CAF_strand_given = 0;

  strict_seq_vector = true;
  cleanupContig();
}


CAF::~CAF() {
  discard();
}



void CAF::setStrictSeqVector(bool strict)
{
  strict_seq_vector = strict;
}



void CAF::discard() {
  CAF_state.clear();
  CAF_readname.clear();
  CAF_asped_date.clear();
  CAF_scf_file.clear();
  CAF_primer.clear();
  CAF_template.clear();
  CAF_ligation.clear();
  CAF_clone_vector_text.clear();
  CAF_base_caller.clear();
  CAF_stolen.clear();
  CAF_clone.clear();
  CAF_sequencing_vector.clear();

  cleanupContig();
}


void CAF::cleanupContig() {
  cleanup();
  CAF_assembledFrom.clear();
}


void CAF::cleanup() {
  caf_intervall i;

  CAF_type   = caftype_unknown;
  CAF_dye    = DYE_UNDEFINED;
  CAF_length = 0;
  CAF_staden_id = 0;

  CAF_qual_clip_max   = -1;
  CAF_qual_clip_min   = -1;
  CAF_insert_size_min = -1;
  CAF_insert_size_max = -1;
  CAF_strand_given    = 0;
  CAF_clone_max = -1;
  CAF_clone_min = -1;
  CAF_seq_max = -1;
  CAF_seq_min = -1;
  CAF_dye     = -1;

  CAF_align_scf.clear();
  CAF_DNA.clear();
  CAF_quality.clear();
  CAF_unprocessed.clear();
  CAF_taglist.clear();
  CAF_assembledFrom.clear();

  CAF_no_seq_intervall.clear();
  i.leftBorder = -1;
  i.rightBorder = -1;
  CAF_no_seq_intervall.push_back(i);

  CAF_state.clear();
  CAF_readname.clear();
  CAF_asped_date.clear();
  CAF_scf_file.clear();
  CAF_primer.clear();
  CAF_template.clear();
  CAF_ligation.clear();
  CAF_clone_vector_text.clear();
  CAF_base_caller.clear();
  CAF_stolen.clear();
  CAF_clone.clear();
  CAF_sequencing_vector.clear();
}



// ------------------------------------------------
//      Main-loop for consuming the CAF-file
// ------------------------------------------------


// Read the next token; mind the reReadToken flag

caf_token CAF::readToken() {

  if (reReadToken == 0) {
    token = static_cast<caf_token>(lexer->yylex());
  }

  if (verbose && reReadToken == 0) {
    cout << "Read Token: " << static_cast<int32>(token) << " "
	 << lexer->YYText() << endl;
  }

  reReadToken = 0;

  return token;
}


// Load the CAF-File fileName into the ReadPool and the Contig-List of the
// CAF_object.

size_t CAF::countReadsBeforeLoad(const string & fileName)
{
  FUNCSTART("void CAF::countReadsBeforeLoad()");

  ifstream cafin;

  size_t numseqsloaded=0;

  cafin.open(fileName.c_str(), ios::in|ios::ate);
  if(!cafin) {
    MIRANOTIFY(Notify::FATAL, "CAF file not found for loading: " << fileName);
  }
  if(!cafin.tellg() ) {
    MIRANOTIFY(Notify::FATAL, "CAF file is empty: " << fileName);
  }

  ProgressIndicator<streamsize> P(0, cafin.tellg(),5000);

  cafin.seekg(0, ios::beg);

  string actline;
  string acttoken;

  while(!cafin.eof()){
    if(cafin.eof()) break;
    getline(cafin,actline);

    string::size_type tokenstart=actline.find_first_not_of(" \t");
    if(tokenstart!=string::npos){
      string::size_type tokenend=actline.find_first_of(" \t:", tokenstart);
      if(tokenend==string::npos) {
	acttoken=actline.substr(tokenstart,tokenend);
      }else{
	acttoken=actline.substr(tokenstart,tokenend-tokenstart);
      }
      if(acttoken=="Is_read") numseqsloaded++;
    }
    if(P.delaytrigger()) P.progress(cafin.tellg());
  }
  P.finishAtOnce();

  cafin.close();

  FUNCEND();

  return numseqsloaded;
}


/*
  seqtype = default seqtype of sequences if not encoded in the CAF
  loadaction:
    //  0 = count only
    //  1 = load
  lrperseqtype = longest read per seqtype

  returns:
    1) number of sequences loaded
    2) if loading: size of longest read per seqtype in lrperseqtype
 */
size_t CAF::load(const string & fileName, const uint8 seqtype, const uint8 loadaction, vector<uint32> & lrperseqtype, bool recalcconsensus, void (*ccallback)(list<Contig> &, ReadPool &), void (*rcallback)(ReadPool &), bool isVerbose)
{
  FUNCSTART("void CAF::load()");

  BUGIFTHROW(loadaction>1,"loadaction>1??");

  //cout << "gnuuu: " << lrperseqtype.size() << endl;
  //cout << (uint16) loadaction << endl;

  theCCallback=ccallback;
  theRCallback=rcallback;
  CAF_recalcconsensus=recalcconsensus;

  size_t numseqsloaded=0;
  if(loadaction==0){
    cout << "Counting reads:\n";
    numseqsloaded=countReadsBeforeLoad(fileName);
    return numseqsloaded;
  }

  cout << "\nLoading CAF " << fileName << " :\n";

  CAF_defaultseqtype=seqtype;
  CAF_lrperseqtype.clear();
  CAF_lrperseqtype.resize(ReadGroupLib::SEQTYPE_END,0);

  verbose = isVerbose;
  //verbose=true;

  CAF_listofreadnames.clear();
  CAF_hashedreadnames.clear();

  ifstream cafFile;

  cafFile.open(fileName.c_str(), ios::in|ios::ate);
  if(!cafFile) {
    MIRANOTIFY(Notify::FATAL, "CAF file not found for loading: " << fileName);
  }
  if(!cafFile.tellg() ) {
    MIRANOTIFY(Notify::FATAL, "CAF file is empty: " << fileName);
  }

  ProgressIndicator<streamsize> P(0, cafFile.tellg(),5000);

  cafFile.seekg(0, ios::beg);

  lexer = new CAFFlexLexer(&cafFile);
  if (lexer == nullptr) {
    MIRANOTIFY(Notify::WARNING, "Unable to create flexer!");
  }

  reReadToken = 0;


  while(0 != (token = readToken())) {
    if(P.delaytrigger()) P.progress(cafFile.tellg());

    switch(token) {

    case token_sequencename:  {
      getCafReadname();
      CAF_type = caftype_read;   break;
    }
    case token_type_read:{
      numseqsloaded++;
      CAF_type = caftype_read;   break;
    }
    case token_type_contig:   CAF_type = caftype_contig; break;
    case token_type_assembly: CAF_type = caftype_assembly; break;
    case token_type_group:    CAF_type = caftype_group;    break;
    case token_padded:        CAF_isPadded = true; break;
    case token_unpadded:      CAF_isPadded = false; break;
    case token_pstatus:       getCafProcessStatus(); break;
    case token_asped:         getCafAspedDate();  break;
    case token_dye:           getCafDye();        break;
    case token_scf_file:      getCafSCFFile();    break;
    case token_primer:        getCafPrimer();     break;
    case token_template:      getCafTemplate();   break;
    case token_insert_size:   getCafInsertSize(); break;
    case token_ligation:      getCafLigation();   break;
    case token_stolen:        getCafStolen();     break;
    case token_forward:       CAF_strand_given = 1; break;
    case token_reverse:       CAF_strand_given = -1; break;
    case token_base_caller:   getCafBaseCaller();   break;
    case token_seq_vector:    getCafSeqVector();    break;
    case token_sequencing_vector: getCafSequencingVector(); break;
    case token_clone_vector:  getCafCloneVector();  break;
    case token_clipping:      getCafClipping();   break;
    case token_sequence:      getCafSequence();   break;
    case token_dna:           getCafDNA();        break;
    case token_quality:       getCafQuality();    break;
    case token_align_SCF:     getCafAlignSCF();   break;
    case token_tag:           getCafTag();        break;
    case token_assembled:     getCafAssembledFrom(); break;
    case token_staden_id:     getCafStadenId();    break;
    case token_clone:         getCafClone();       break;
    case token_error:         {
      cerr << "ARGH!: " << lexer->YYText() << endl;
      throw  Notify(Notify::FATAL,
		    "CAF::load(const char *filename, bool isVerbose)",
		    "Illegal character in CAF-file!");
    }
    default: addUnprocessed();
    }
  }

  createCafObject();
  delete lexer;

  P.finishAtOnce();
  cout << endl;

  lrperseqtype=CAF_lrperseqtype;

  FUNCEND();

  return numseqsloaded;
}


// --------------------------------
//          Hilfsroutinen
// --------------------------------

/* Read the next token which is not optional */
/* Does not check for reReadToken flag!        */

caf_token CAF::getNextToken() {
  caf_token new_token;

  FUNCSTART("CAF::getNextToken()");

  new_token = static_cast<caf_token>(lexer->yylex());
  if (new_token == 0) {
    MIRANOTIFY(Notify::FATAL, "Unexpected end of CAF-File");
  }

  FUNCEND();
  return new_token;
}


// Read an identifier or a quoted text; allocate memory for it and return
// a pointer to the text

string CAF::getText()
{
  string name;

  FUNCSTART("CAF::getText()");

  token = getNextToken();
  if (token == token_identifier  ||
      token == token_quoted_text) {
    name = lexer->YYText();
    deescapeString(name);

    FUNCEND();
    return name;
  }

  FUNCEND();
  return "";
}


// Read an identifier and allocate memory for it; return pointer to
// text.

char * CAF::getIdentifier()
{
  caf_token new_token;
  char *name;

  FUNCSTART("CAF::getIdentifier()");

  new_token = getNextToken();

  if (new_token == token_identifier || new_token == token_number) {
    name = new char[lexer->YYLeng() + 1];
    // why is valgrind moaning on strlen(lexer->YYText()) here?
    //cout << "Token string (" << lexer->YYLeng() << "): " << lexer->YYText() << endl;
    //cout << (char) lexer->YYText()[0];
    //cout << (char) lexer->YYText()[1];
    //cout << (char) lexer->YYText()[2];
    //cout << (char) lexer->YYText()[3];
    //cout << (char) lexer->YYText()[4];
    //cout << (char) lexer->YYText()[5];
    //cout << (char) lexer->YYText()[6];
    //cout << (char) lexer->YYText()[7];
    //cout << (char) lexer->YYText()[8];
    //cout << (char) lexer->YYText()[9];
    //cout << (char) lexer->YYText()[10];
    //cout << (char) lexer->YYText()[11];
    //cout << (char) lexer->YYText()[12];
    //cout << (char) lexer->YYText()[13];
    //cout << (char) lexer->YYText()[14];
    //cout << (char) lexer->YYText()[15];
    //cout << (char) lexer->YYText()[16];
    //cout << strlen(lexer->YYText()) << endl;
    strcpy(name, lexer->YYText());

    return name;
  }

  cout << "Token found: " << new_token << endl;
  cout << "Token string: " << lexer->YYText() << endl;
  //delete [] name;
  MIRANOTIFY(Notify::FATAL, ": Expected Name");

  FUNCEND();

}


// Reads an identifier and compares it with the text in "old_name"
// Call "dump" if the value of old_name has changed. Returns the new name.
// old_name ist deleted!

string CAF::getReadname(const string old_name) {
  string new_name;

  FUNCSTART("CAF::getReadname(char *old_name)");

  new_name = getIdentifier();

  if (new_name.empty()) {
    return "";
  }

  if (old_name.empty())
    return new_name;

  if (new_name!=old_name) {
    if (verbose) {
      cout << "name has changed : " << old_name << "->" << new_name << endl;
    }

    createCafObject();

    return new_name;
  }

  return old_name;
}


// Reads two numbers an throws an exception if there are any problems.

int32 CAF::get2Numbers(int32 &x1, int32 &x2)
{
  FUNCSTART("CAF::get2Numbers()");

  if (token_number == token || token_quoted_text == token) {
    x1 = atoi(lexer->YYText());
  } else {
    MIRANOTIFY(Notify::FATAL, "Expected a number");
  }

  token = getNextToken();
  if (token_number == token || token_quoted_text == token) {
    x2 = atoi(lexer->YYText());
  } else {
    MIRANOTIFY(Notify::FATAL, "Expected a second number");
  }

  if (verbose) {
    cout << "Get Numbers " << x1 << " " << x2 << endl;
  }
  FUNCEND();
  return 0;
}


void CAF::addUnprocessed() {
  static int line = 0;

  FUNCSTART("CAF::addUnprocessed()");
  line = lexer->lineno();

  while (line == lexer->lineno()) {
    for (int loop=0; loop<lexer->YYLeng(); loop++) {
      CAF_unprocessed.push_back(lexer->YYText()[loop]);
    }
    CAF_unprocessed.push_back(' ');
    readToken();
  }
  CAF_unprocessed.push_back('\n');

  reReadToken = 1;
  FUNCEND();
}



// ---------------------------------------
//   Behandlung der Schlüsselworte
// ---------------------------------------

int32 CAF::getCafReadname()
{
  CAF_readname = getReadname(CAF_readname);


  return (CAF_readname.empty());
}


int32 CAF::getCafSCFFile()
{
  CAF_scf_file = getText();
  return (CAF_scf_file.empty());
}


int32 CAF::getCafAspedDate()
{
  CAF_asped_date = getText();
  return (CAF_asped_date.empty());
}


int32 CAF::getCafBaseCaller()
{
  CAF_base_caller = getText();
  return (CAF_base_caller.empty());
}


int32 CAF::getCafClone()
{
 CAF_clone = getText();
 return (CAF_clone.empty());
}


int32 CAF::getCafSequencingVector()
{
  CAF_sequencing_vector = getText();

  return (CAF_sequencing_vector.empty());
}


int32 CAF::getCafStolen()
{
  CAF_stolen = getText();
  return (CAF_stolen.empty());
}


int CAF::getCafTemplate()
{
  CAF_template = getText();
  return (CAF_template.empty());
}


int32 CAF::getCafStadenId()
{
  caf_token new_token;

  FUNCSTART("CAF::getCafStadenId()");

  new_token = static_cast<caf_token>(lexer->yylex());

  if (token_number == new_token || token_identifier == new_token || token_quoted_text == new_token) {
    CAF_staden_id = atoi(lexer->YYText());
  } else {
    cout <<lexer->YYText() << endl;
    MIRANOTIFY(Notify::FATAL, "Expected a number");
  }

  FUNCEND();
  return CAF_staden_id;
}


int32 CAF::getCafProcessStatus() {
  FUNCSTART("CAF::getCafProcessStatus()");

  CAF_state = getText();
  FUNCEND();
  return (CAF_state.empty());
}


int32 CAF::getCafDye() {
  caf_token new_token;
  new_token = getNextToken();

  if (new_token == token_dye_primer) {
     CAF_dye = DYE_PRIMER;
     return 0;
  }
  if (new_token == token_dye_terminator) {
     CAF_dye = DYE_TERMINATOR;
     return 0;
  }
  return 1;
}


int32 CAF::getCafPrimer()
{
   CAF_primer = getText();

   return 0;
}


int32 CAF::getCafInsertSize()
{
   token = getNextToken();
   return get2Numbers(CAF_insert_size_min, CAF_insert_size_max);
}


int32 CAF::getCafLigation()
{
  FUNCSTART("int32 CAF::getCafLigation()");

  CAF_ligation = getText();

  FUNCEND();
  return (CAF_ligation.empty());
}


int32 CAF::getCafDNA()
{
  CAF_readname = getReadname(CAF_readname);
  if (CAF_readname.empty()) {
    return 1;
  }

  return 0;
}


int32 CAF::getCafSeqVector()
{
  int dummy_min, dummy_max;
  list<caf_intervall>::iterator I = CAF_no_seq_intervall.begin();
  list<caf_intervall> newList;


  FUNCSTART("CAF::getCafSeqVector()");

  token = getNextToken();
  if (token == token_identifier || token == token_state) {
    token = getNextToken(); // ignore type
  }

  get2Numbers(dummy_min, dummy_max);

  if (dummy_min > dummy_max) swap(dummy_min,dummy_max);

  /*

  BaCh: 23.10.2007

  Bug triggered by 454 data

  This was some code which Thomas and I put in to catch broken SVEC
  tags of the type (assuming a sequence length of, say, 1000): SVEC 3 4

  The lines below assume that no one would tag SVEC with only two bases.
  Well, 454 does not do it neither, but they currently also do not
  distinguish between QUAL clippings and SVEC clippings. And the
  NCBI therefore sets SVEC *and* QUAL clippings (same values), but this
  results into SVECs with a length of <=2.

  The reason that we did not check whether the coordinates of SVEC reach
  boundaries of a sequence is that officially the CAF model does not define
  any order of elements in a CAF file. That is, theoretically, "Sequence :"
  (and with it Seq_vec) could appear before "DNA :" or "BaseQuality :" ... and
  then the length of the sequence itself would not be known.

  Anyway, I'm taking this out now and define SVECs in the middle of the
  sequence to be "broken" CAF. That is, SVECs that are not reaching to either
  the left or right border of the sequence are invalid and MIRA not
  responsible for any mishaps that result out of it.

  // Ignore SeqVector definitions in the middle of the sequence
  if (dummy_max - dummy_min < 2 && dummy_min > 1)
    return 0;
  */

  newList.clear();

  if (dummy_min == 1) {
    // We have only one good Intervall [dummy_max...oo]
    caf_intervall neu, schnitt;

    neu.leftBorder  = dummy_max;
    neu.rightBorder = -1;

    while (I != CAF_no_seq_intervall.end()) {
      schnitt = schnittmenge(*I, neu);
      if (!intervall_leer(schnitt)) {
	newList.push_back(schnitt);
      }
      I++;
    }
  } else {
    // We have the SeqVec Intervalls [0...dummy_min-1] and [dummy_max...oo]
    caf_intervall neu_links, neu_rechts, schnitt;

    neu_links.leftBorder   = 0;
    neu_links.rightBorder  = dummy_min-1;
    neu_rechts.leftBorder  = dummy_max;
    neu_rechts.rightBorder = -1;

    while (I != CAF_no_seq_intervall.end()) {
      schnitt = schnittmenge(*I, neu_links);

      if (!intervall_leer(schnitt)) {
	newList.push_back(schnitt);
      }
      schnitt = schnittmenge(*I, neu_rechts);
      if (!intervall_leer(schnitt)) {
	newList.push_back(schnitt);
      }
      I++;
    }
  }

  CAF_no_seq_intervall.swap(newList);
  newList.clear();

  if (verbose) {
    cout << endl;
    showIntervallList(CAF_no_seq_intervall, cout);
  }

  if (dummy_min == 1) {
    CAF_seq_min = dummy_max;
  } else {
    CAF_seq_max = dummy_min - 1;
  }


  token = getNextToken();
  if (token != token_quoted_text && token != token_identifier) {
    reReadToken = 1;
  } else {
    reReadToken = 0;
  }
  FUNCEND();
  return 0;
}



int32 CAF::getCafCloneVector() {
  int dummy_min, dummy_max;

  token = getNextToken();
  //  if (token == token_identifier)
  //   token = getNextToken();  //ignore Type;

  get2Numbers(dummy_min, dummy_max);

  if (dummy_min == 1) {
    CAF_clone_min = dummy_max;
  } else {
    CAF_clone_max = dummy_min - 1;
  }

  token = getNextToken();
  if (token != token_quoted_text && token != token_identifier) {
    reReadToken = 1;
  } else {

    CAF_clone_vector_text = lexer->YYText();
    deescapeString(CAF_clone_vector_text);
  }
  return 0;
}


int32 CAF::getCafClipping() {
  int dummy_min, dummy_max;


  FUNCSTART("CAF::getCafClipping()");

  token = getNextToken();

  get2Numbers(dummy_min, dummy_max);

  if (CAF_qual_clip_min == -1 || CAF_qual_clip_min < dummy_min) {
    CAF_qual_clip_min = dummy_min;
  }
  if (CAF_qual_clip_max == -1 || CAF_qual_clip_max > dummy_max) {
    CAF_qual_clip_max = dummy_max;
  }

  if (CAF_qual_clip_min == 0) {
    cerr << CAF_readname << ": quality clip left is zero??? Should not be, base numbering in CAF starts at 1. Changed to 1, but you DO want to check this." << endl;
    CAF_qual_clip_min = 1;
  }
  if (CAF_qual_clip_max == 0) {
    cerr << CAF_readname << ": quality clip right is zero??? Should not be, base numbering in CAF starts at 1. Changed to 1, but you DO want to check this." << endl;
    CAF_qual_clip_max = 1;
  }

  token = getNextToken(); // quoted_text

  FUNCEND();
  return 0;
}

// Einlesen der Sequenz. Es wird so lange gelesen, bis ein Schlüsselwort
// kommt. Der Sequenzname ersetzt einen evtl. vorher gesetzten Namen.
// Sequenzen die länger als MAX_SEQUENCE_LENGTH sind können nicht eingelesen
// werden -> Align_to_SCF

int32 CAF::getCafSequence() {

  CAF_readname = getReadname(CAF_readname);
  if (CAF_readname.empty()) {
    return 1;
  }

  token = getNextToken();

  CAF_DNA.clear();

  while (token == token_identifier) {
    // BaCh: added IUPAC bases
    for (int32 loop = 0; loop < lexer->YYLeng(); loop++) {
      if(dptools::isValidIUPACBase(lexer->YYText()[loop])) {
	CAF_DNA.push_back((lexer->YYText())[loop]);
      }else if((lexer->YYText())[loop] == '-'){
	CAF_DNA.push_back('*');
      }else{
	CAF_DNA.push_back('N');
      }
    }
    token = readToken();
  }

  CAF_length = CAF_DNA.size();

  reReadToken = (token != token_ende);
  return 0;
}


// Einlesen der Liste der Basenqualitäten. Es wird so lange gelesen bis
// keine Zahl mehr kommt und das zuviel gelesene Token zurückgesetzt und nicht
// bis eine Leerzeile kommt!
// Der Sequenzname ersetzt einen evtl. vorher (beim DNA Einlesen) gelesenen
// Namen.

int32 CAF::getCafQuality() {
  CAF_quality.clear();
  CAF_quality.reserve(CAF_length+5);
  CAF_readname = getReadname(CAF_readname);
  if (CAF_readname.empty())
    return 1;

  token = getNextToken();

  while (token == token_number) {
    CAF_quality.push_back(static_cast<uint8>(atoi(lexer->YYText())));
    token = getNextToken();
  }

  reReadToken = (token != token_ende);
  return 0;
}


// Verarbeiten einer Align_to_SCF Zeile. Die dort stehede Zuordnung wird
// in das CAF_align_scf-Array eingetragen. Dies passiert noch bevor die
// Sequenz überhaupt gelesen wurde und damit die Länge ermittelt werden
// könnte! Daher wird mit einem Array fester Größe gearbeitet.

int32 CAF::getCafAlignSCF()
{
  FUNCSTART("int32 CAF::getCafAlignSCF()");

  int seq_von, seq_bis, scf_von, scf_bis;
  int seq_iterator = 1;
  int scf_iterator = 1;
  int anzahl;

  token = getNextToken();
  get2Numbers(seq_von, seq_bis);
  token = getNextToken();
  get2Numbers(scf_von, scf_bis);

  if (seq_bis < seq_von)
    seq_iterator = -1;
  if (scf_bis < scf_von)
    scf_iterator = -1;

  if (abs(seq_bis - seq_von) != abs(scf_bis - scf_von)) {
    return 0;
  }

  if (max(seq_von, seq_bis) > static_cast<int32>(CAF_align_scf.size())) {
    CAF_align_scf.resize(max(seq_von, seq_bis), -1);
  }

  anzahl = abs(seq_bis - seq_von) + 1;
  int seqi=seq_von;
  int scfi=scf_von;
  for (int loop = 0; loop < anzahl; loop++) {
    if(seqi-1 <0 || seqi-1 >= CAF_align_scf.size()) {
      cout << "While reading read: " << CAF_readname << endl;
      cout << "While reading line: Align_to_SCF " << seq_von << " " << seq_bis << " " <<scf_von << " " << scf_bis << endl;
      cout << "This makes for an illegal alignment of the sequence to the scf, wrong values in this line?\n";
      MIRANOTIFY(Notify::FATAL, "Illegal alignment to SCF: " << CAF_readname);
    }
    CAF_align_scf[seqi-1] = scfi - 1;
    seqi += seq_iterator;
    scfi += scf_iterator;
  }
  return 0;
}


/*************************************************************************
 *
 * function to replace the strings \n, \" etc in the string by their
 *  true char equivalent
 *
 * well, just for the main escaped chars
 *
 *************************************************************************/
void CAF::deescapeString(string & s)
{
  uint32 from=0;
  uint32 to=0;

  for(; from<s.size(); from++, to++) {
    if(s[from]=='\\'){
      switch(s[from+1]) {
      case '"' : {
	s[to]='"';
	from++;
	break;
      }
      case 'n' : {
	s[to]='\n';
	from++;
	break;
      }
      case '\\' : {
	s[to]='\\';
	from++;
	break;
      case 't' : {
	s[to]='\t';
	from++;
	break;
      }
      }
      default : {
	s[to]=s[from];
      }
      }
    }else{
      s[to]=s[from];
    }
  }

  // shorten string to the required size
  s.resize(s.size()-(from-to));
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

int32 CAF::getCafTag() {
  multitag_t tmptag;
  int   tag_von, tag_bis;

  FUNCSTART("CAF::getCafTag()");

  token = getNextToken();

  if (token != token_number) {
    string tmps(lexer->YYText());
    string identifier(AnnotationMappings::translateGAP4feat2SOfeat(tmps));
    if(identifier.empty()){
      identifier=AnnotationMappings::translateXGAP4feat2SOfeat(tmps);
      if(identifier.empty()) identifier=tmps;
    }
    tmptag.setIdentifierStr(identifier);
    token = getNextToken();
  }

  get2Numbers(tag_von, tag_bis);
  // Verarbeite die Positionen
  tmptag.from = tag_von - 1;
  tmptag.to   = tag_bis - 1;

  if(tmptag.from>tmptag.to) {
    tmptag.setStrand('-');
    swap(tmptag.from, tmptag.to);
  }else{
    tmptag.setStrand('+');
  }

  if (verbose) {
    cout << "Tag Positition: " << tag_von << "-" <<  tag_bis << endl;
  }

  token = getNextToken();

  if (token == token_quoted_text) {
    if (lexer->YYLeng() < 1) {
      reReadToken = 1;
    } else {
      // Verarbeite Kommentar
      string tmps(lexer->YYText());
      deescapeString(tmps);
      tmptag.setCommentStr(tmps);
    }
  } else {
    if (verbose) {
      cout << "Read token : " << token << endl;
      cout << "YYText: " << lexer->YYText() << endl;
    }
    reReadToken = 1;
  }

  multitag_t newtag(tmptag);
  Read::upgradeOldTagToMultitagWithGFF3(tmptag,newtag);
  CAF_taglist.push_back(newtag);

  FUNCEND();
  return 0;
}


int32 CAF::getCafAssembledFrom()
{
  FUNCSTART("int32 CAF::getCafAssembledFrom()");
  char *readname;
  int  contig_von, contig_bis, read_von, read_bis;
  int richtung;
  Contig::contig_init_read_t *contigRead;


  readname = getIdentifier();
  if (readname == nullptr) return 1;

  token = getNextToken();
  get2Numbers(contig_von, contig_bis);

  contigRead = new Contig::contig_init_read_t;

  if (contig_von > contig_bis) {
    richtung = -1;
    contigRead->offset_in_contig = contig_bis-1;
 } else {
    richtung = 1;
    contigRead->offset_in_contig = contig_von-1;
  }

  token = getNextToken();
  get2Numbers(read_von, read_bis);

  if (read_von > read_bis) {
    contigRead->read_rclip = read_von;
    contigRead->read_lclip = read_bis-1;
    contigRead->direction = - richtung;
  } else {
    contigRead->read_rclip = read_bis;
    contigRead->read_lclip = read_von-1;
    contigRead->direction = richtung;
  }

  try {
    contigRead->id=-1;
    {
      strmap::iterator rnI;

      rnI=CAF_hashedreadnames.find(readname);
      if(rnI==CAF_hashedreadnames.end()){
	cerr << "Searched for read " << readname << " but did not find it?" << endl;
	MIRANOTIFY(Notify::INTERNAL, "Readname not found in hash as expected?"<< readname);
      }else{
	contigRead->id=rnI->second;
      }
    }

    if(contigRead->id < 0) {
      MIRANOTIFY(Notify::INTERNAL, "Read ID not initialised: " << readname);
    }

    CAF_assembledFrom.push_back(*contigRead);
  }
  catch(Notify n) {
    cout << "Unable to find Read in Pool\n";
    n.handleError("int32 CAF::getCafAssembledFrom()");
    return 0;
  }

  delete [] readname;
  delete contigRead;

  return 0;
}



Read & CAF::createCafRead()
{
  FUNCSTART("Read & CAF::createCafRead()");

  if (verbose) {
    cout << "Searching read in pool." << endl;
    cout.flush();
  }

  {
    uint32 readid=thePool->size();
    strmap::iterator rnI;

    rnI=CAF_hashedreadnames.find(CAF_readname);
    if(rnI!=CAF_hashedreadnames.end()){
      cerr << "Read " << CAF_readname << " already loaded before!" << endl;
      MIRANOTIFY(Notify::FATAL, "Duplicate readname in CAF file: " << CAF_readname);
    }else{
      CAF_listofreadnames.push_back(CAF_readname);
      CAF_hashedreadnames[CAF_readname]=readid;
    }
  }

  if (verbose) {
    cout << "Adding read to pool." << endl;
    cout.flush();
  }
  Read &aRead = thePool->getRead(thePool->provideEmptyRead());
  if (verbose) {
    cout << "Read added." << endl;
    cout.flush();
  }

  if (CAF_length == 0) {
    CAF_length = CAF_DNA.size();
  }

  caf_intervall neu={0,0};
  caf_intervall schnitt=neu;
  list<caf_intervall>::iterator I = CAF_no_seq_intervall.begin();

  neu.leftBorder  = CAF_qual_clip_min-1;
  neu.rightBorder = CAF_qual_clip_max;

  // Start Version B: Ignore SeqVecs not at the ends

  if (false) { // true == strict_seq_vector) {
    caf_intervall alles;
    alles.leftBorder  = neu.rightBorder;
    alles.rightBorder = neu.leftBorder;

    while (I != CAF_no_seq_intervall.end()) {
      schnitt = schnittmenge(*I, neu);
      if (!intervall_leer(schnitt)) {
	if (alles.leftBorder > schnitt.leftBorder)
	  alles.leftBorder = schnitt.leftBorder;
	if (alles.rightBorder < schnitt.rightBorder)
	  alles.rightBorder = schnitt.rightBorder;
      }
      I++;
    }

    //showIntervall(alles, cout);
    //cout << endl;

    schnitt = schnittmenge(alles, neu);

    //showIntervall(schnitt, cout);
    //cout << endl;

    if (!intervall_leer(schnitt)) {
      CAF_seq_min = schnitt.leftBorder;
      CAF_seq_max = schnitt.rightBorder;
    }

  } else {
    // Version A: Select first not empty intervall

    while (I != CAF_no_seq_intervall.end()) {
      schnitt = schnittmenge(*I, neu);
      if (!intervall_leer(schnitt)) break;
      I++;
    }

    if (!intervall_leer(schnitt)) {
      CAF_seq_min = I->leftBorder;
      CAF_seq_max = I->rightBorder;
    }
  }


  if (CAF_clone_max == -1) CAF_clone_max = CAF_length;
  if (CAF_clone_min == -1) CAF_clone_min = 0;


  CAF_align_scf.resize(CAF_length, -1);
  //CAF_quality.resize(CAF_length, 1);
  //cout << "Check: " << CAF_DNA.size() << "\t" << CAF_quality.size() << endl;
  //exit(0);
  //CAF_quality.resize(CAF_length-1);

  if (CAF_seq_max == -1) CAF_seq_max = CAF_length;
  if (CAF_seq_min <  0) CAF_seq_min = 0;
  //if (CAF_seq_min <   1) CAF_seq_min = 1;


  if (verbose) {
    cout << "Create Read " << CAF_readname <<endl;
    cout << "clip left: " << CAF_qual_clip_min -1
	 << "\t clip right: " << CAF_qual_clip_max << endl;
    cout << "Seq left:  " << CAF_seq_min
	 << "\t Seq right : " << CAF_seq_max << endl;

    vector<multitag_t>::iterator J = CAF_taglist.begin();
    while (J != CAF_taglist.end()) {
      cout << "Tag " << J->from << "  " << J->to << endl;
      cout << "Comment: " << J->getCommentStr() << endl;
      J++;
    }
    //    dump(cout);
  }

  // BaCh 04.06.06
  //  Problem: if SVEC is "visible in alignment, this will complete
  //  destroy the assembly reconstructed lateron as the read object
  //  automatically (and understandably) hides away SVEC
  // Two possible resolutions:
  //  - extend qualclipmin back to svec clip
  //  - cut down svec clip to qualclip min
  //
  // Chosen the later (if SVEC was uncovered during editing in GAP4 or
  //  other, then probably because people know what they're doing
  //  (at least I hope so)

  // BaCh 25.06.2006
  // actually, a bad bad bad decision, this breaks a number of different
  //  assertions made by gap4/gap2caf
  //
  // Thomas and I seem to have run upon that kind of problem earlier,
  //  the "initialiseContig()" function has traces of tries for resolving
  //  that problem, but there is no real good solution
  //
  // Only resolution: sorry, none. Just forbid people to have non-covered
  //  SVEC sequence in gap4 projects, no other possibility.
  //
  // if(CAF_seq_min > CAF_qual_clip_min) CAF_seq_min=CAF_qual_clip_min-1;


  // BaCh 13.05.2011
  // Deleting consensus in gap4 and with that completely deleting reads leaves reads of length
  //  0 in the database. gap2caf will convert that to reads having a length of 1, but will get
  //  the quality clips wrong: "2 1"  (left clip > right clip and leftclip 2 > length of read)
  //
  // Therefore, check for these kind of things.

  if(CAF_qual_clip_min>CAF_DNA.size()){
    cout << "WARNING: sequence " << CAF_readname << " in CAF file has quality clip left > length of sequence, that should not be!\n";
    CAF_qual_clip_min=CAF_DNA.size();
  }
  if(CAF_qual_clip_max>CAF_DNA.size()){
    cout << "WARNING: sequence " << CAF_readname << " in CAF file has quality clip right > length of sequence, that should not be!\n";
    CAF_qual_clip_max=CAF_DNA.size();
  }
  if(CAF_seq_min>CAF_DNA.size()){
    cout << "WARNING: sequence " << CAF_readname << " in CAF file has seq clip left > length of sequence, that should not be!\n";
    CAF_seq_min=CAF_DNA.size();
  }
  if(CAF_seq_max>CAF_DNA.size()){
    cout << "WARNING: sequence " << CAF_readname << " in CAF file has seq clip right > length of sequence, that should not be!\n";
    CAF_seq_max=CAF_DNA.size();
  }
  if(CAF_clone_min>CAF_DNA.size()){
    cout << "WARNING: sequence " << CAF_readname << " in CAF file has clone clip left > length of sequence, that should not be!\n";
    CAF_clone_min=CAF_DNA.size();
  }
  if(CAF_clone_max>CAF_DNA.size()){
    cout << "WARNING: sequence " << CAF_readname << " in CAF file has clone clip right > length of sequence, that should not be!\n";
    CAF_clone_max=CAF_DNA.size();
  }

  try {
//    aRead.setReadNamingScheme((*CAF_miraparams)[aRead.getSequencingType()].getAssemblyParams().as_readnaming_scheme);

    // TODO: further adapt to read groups?
//    if (!CAF_clone_vector_text.empty()) {
//      aRead.setClonevecName(CAF_clone_vector_text);
//    }
//
//    if (CAF_dye == DYE_PRIMER) {
//      aRead.setDye("Dye_primer");
//    }
//    if (CAF_dye == DYE_TERMINATOR) {
//      aRead.setDye("Dye_terminator");
//    }
//

    string minft_strainname;
    string minft_seqtypename;
    string minft_machinetype;
    int8 minft_splacementcode;
    bool minft_isbb;
    bool minft_israil;
    bool minft_isCER;

    // treat CAF clones as strains
    if (!CAF_clone.empty() && CAF_clone != "unknown") {
      minft_strainname=CAF_clone;
    }


    (void) Read::extractMINFTagInfo(CAF_taglist,
				    CAF_readname,
				    minft_strainname,
				    minft_seqtypename,
				    minft_machinetype,
				    minft_splacementcode,
				    minft_isbb,
				    minft_israil,
				    minft_isCER);


    string dummy_empty;

    uint8 st=ReadGroupLib::stringToSeqType(minft_seqtypename);
    if(st==ReadGroupLib::SEQTYPE_END) st=ReadGroupLib::SEQTYPE_TEXT;

    ReadGroupLib::ReadGroupID rgid=ReadGroupLib::searchExactRGMatch(
      dummy_empty,
      st,
      CAF_insert_size_min,
      CAF_insert_size_max,
      minft_splacementcode,
      minft_strainname,
      minft_isbb,
      minft_israil,
      minft_isCER,
      CAF_sequencing_vector,
      dummy_empty,
      CAF_base_caller);

    if(rgid.isDefaultNonValidReadGroupID()){
      rgid=ReadGroupLib::newReadGroup();
      rgid.setGroupName(dummy_empty);
      rgid.setSequencingType(st);
      rgid.setInsizeFrom(CAF_insert_size_min);
      rgid.setInsizeTo(CAF_insert_size_max);
      rgid.setSegmentPlacementCode(minft_splacementcode);
      rgid.setStrainName(minft_strainname);
      rgid.setBackbone(minft_isbb);
      rgid.setRail(minft_israil);
      rgid.setCoverageEquivalentRead(minft_isCER);
      rgid.setSeqVecName(CAF_sequencing_vector);
      rgid.setMachineType(minft_machinetype);
      rgid.setBaseCaller(CAF_base_caller);

      rgid.setCloneVecName(CAF_clone_vector_text);
      rgid.setPrimer(CAF_primer);
    }

    bool noqual=false;
    if(CAF_quality.empty()){
      noqual=true;
      cout << "WARNING: sequence " << CAF_readname << " in CAF file has no qualities given!\n";
      CAF_quality.resize(CAF_DNA.size(), rgid.getDefaultQual());
    }

    aRead.initialiseRead(false,
			 false,
			 CAF_isPadded,  // bool
			 rgid,
			 CAF_DNA, CAF_quality, CAF_align_scf,
			 CAF_taglist,
			 CAF_readname, CAF_scf_file,
			 CAF_qual_clip_min -1, CAF_qual_clip_max,
			 CAF_seq_min, CAF_seq_max,
			 CAF_clone_min, CAF_clone_max);

    //cout << aRead;

    if(CAF_strand_given>0){
      aRead.setTemplateSegment(1);
    }else if(CAF_strand_given<0){
      aRead.setTemplateSegment(255);
    }

    if (!CAF_state.empty()) {
      aRead.setProcessStatus(CAF_state);
    }

    if (!CAF_asped_date.empty()) {
      aRead.setAsped(CAF_asped_date);
    }

    if (!CAF_template.empty()) {
      aRead.setTemplate(CAF_template);
    }
  }
  catch (Notify n) {
    n.handleError("Read & CAF::createCafRead()");
  }

  if(!aRead.isBackbone() && !aRead.isCoverageEquivalentRead()){
    if(aRead.getSequencingType() >= ReadGroupLib::SEQTYPE_END){
      MIRANOTIFY(Notify::FATAL, aRead.getName() << " has unknown sequencing type in MINF tag?\n");
    }
    CAF_lrperseqtype[aRead.getSequencingType()]=max(CAF_lrperseqtype[aRead.getSequencingType()],
						    aRead.getLenClippedSeq());
  }

  cleanup();

  if(theRCallback!=nullptr) {
    (*theRCallback)(*thePool);
  }

  //Read::setCoutType(Read::AS_TEXTSHORT);
  //cout << aRead;

  FUNCEND();
  return aRead;
}


void CAF::createCafContig()
{
  if(theContigs!=nullptr){
    {
      Contig dummy(CAF_miraparams, *thePool);
      theContigs->push_back(dummy);
    }

    if (verbose) {
      cout << "Create Contig " << CAF_readname << "\n";
    }

    // WTF??? BaCh 23.08.10§§
    // This destroys all strand information of consensus tags?!
    // Why was this ever done? And by whom (Thomas, me?)
    //
    //// change all tags in consensus to have strand '='
    //vector<tag_t>::iterator tI=CAF_taglist.begin();
    //for(; tI!=CAF_taglist.end(); tI++){
    //  tI->strand='=';
    //}

    if(CAF_recalcconsensus){
      string dummy1;
      vector<base_quality_t> dummy2;
      theContigs->back().initialiseContig(CAF_assembledFrom,
					  CAF_taglist,
					  CAF_readname,
					  dummy1,dummy2);
    }else{
      string dummy1;
      dummy1.reserve(CAF_DNA.size()+2);
      {
	vector<char>::const_iterator dnaI=CAF_DNA.begin();
	for(; dnaI != CAF_DNA.end(); dnaI++){
	  dummy1+=*dnaI;
	}
      }
      theContigs->back().initialiseContig(CAF_assembledFrom,
					  CAF_taglist,
					  CAF_readname,
					  dummy1,
					  CAF_quality);
    }

    if(theCCallback!=nullptr) {
      (*theCCallback)(*theContigs, *thePool);
    }
  }
  cleanupContig();
}


void CAF::createCafObject()
{

  if (CAF_type == caftype_unknown) {
    cerr << "CAF-Object " << CAF_readname << " is not ready! Type missing.\n";
    return ;
  }

  try {
    switch (CAF_type) {
    case caftype_read: {
      createCafRead();
      break;
    }
    case caftype_contig: {
      createCafContig();
      break;
    }
    default: {
      cout << "No Object created! Object-Type not supported\n";
      cleanup();
    }
    }
  }
  catch (Notify n) {
    n.handleError("Error while creating CAF-Object.\n");
  }

}




ostream &operator<<(ostream &output, CAF const &i)
{
  output << "\nReadname   : " << i.CAF_readname << endl;
  output << "Length     : " << i.CAF_length << endl;
  output << "Objekttyp  : ";
  switch (i.CAF_type) {
  case caftype_read:     output << "read\n";     break;
  case caftype_contig:   output << "contig\n";   break;
  case caftype_group:    output << "group\n";    break;
  case caftype_assembly: output << "assembly\n"; break;
  default: output << "unknown\n";
  }
  output << "Padded     : " << i.CAF_isPadded << endl;
  output << "SCF_File   : " << i.CAF_scf_file << endl;
  output << "Primer     : " << i.CAF_primer << endl;
  output << "Q-Clipping : " << i.CAF_qual_clip_min;
  output << " - " << i.CAF_qual_clip_max << endl;
  output << "Seq-Vector : " << i.CAF_seq_min
	 << " - " << i.CAF_seq_max << endl;
  output << "Cloning V. : " << i.CAF_clone_min << " - "
	 << i.CAF_clone_max << endl;
  output << endl;

  /*
  output << "Sequence : \n";

  output << "Quality: \n";
  vector<uint8>::iterator I = CAF_quality.begin();
  while (I != CAF_quality.end()) {
    cout << (int)*I << " ";
    I++;
  }
  */

  return output;
}
