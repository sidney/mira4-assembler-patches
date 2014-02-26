/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2009 and later by Bastien Chevreux
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
 */


#ifndef _maf_parse_h
#define _maf_parse_h


#include <string>

#include "stdinc/stlincludes.H"

#include "mira/read.H"
#include "mira/readpool.H"
#include "mira/contig.H"



class MAFParse
{
  // Variables
private:
  std::vector<MIRAParameters> *   MAF_miraparams;   // parameters for the contigs

  // Result Variables
  ReadPool*         MAF_readpool;         // where to store the reads
  std::list<Contig>*     MAF_contiglist;      //
  void (*MAF_ccallbackfunc)(std::list<Contig> &, ReadPool &);
  void (*MAF_rcallbackfunc)(ReadPool &);


  std::vector<uint32> MAF_lrperseqtype;  // longest read per seqtype

  bool  MAF_recalcconsensus;

  int32 MAF_vmajor;
  int32 MAF_vminor;

  ReadGroupLib::ReadGroupID MAF_readgroup_rgid;  // rgid of current readgroup
  std::vector<ReadGroupLib::ReadGroupID> MAF_readgroup_externalidmapper;

  std::string MAF_read_name;
  ReadGroupLib::ReadGroupID MAF_read_rgid;  // rgid of current read
  std::string MAF_read_scf_file;
  std::string MAF_read_template;

  std::vector<char>  MAF_read_sequence;
  std::vector<base_quality_t> MAF_read_qualities;
  std::vector<int32> MAF_read_align_origin;
  std::vector<multitag_t> MAF_read_taglist;

  int32 MAF_read_len;

  int32 MAF_read_insert_size_min;
  int32 MAF_read_insert_size_max;
  int32 MAF_read_ql;
  int32 MAF_read_qr;
  int32 MAF_read_cl;
  int32 MAF_read_cr;
  int32 MAF_read_sl;
  int32 MAF_read_sr;

  uint8 MAF_read_tsegment_given;

  ///////////////////// V1 ///////////////////////////
  std::string MAF_read_base_caller;
  std::string MAF_read_sequencing_vector;
  std::string MAF_read_strain;
  std::string MAF_read_machinetype;

  char  MAF_read_strand_given;
  int8  MAF_read_seqtype;

  bool  MAF_read_isbackbone;
  bool  MAF_read_israil;
  bool  MAF_read_isCER;

  bool  MAF_read_seenATline;
  ////////////////////////////////////////////////////



  std::list<Contig::contig_init_read_t> MAF_contig_assembledfrom;
  std::vector<char>  MAF_contig_sequence;
  std::vector<base_quality_t> MAF_contig_qualities;
  std::vector<multitag_t> MAF_contig_taglist;

  std::string MAF_contig_name;

  int32  MAF_contig_numreads;
  int32  MAF_contig_len;      // read in from LC



  bool   MAF_isinread;
  bool   MAF_isincontig;
  bool   MAF_isinreadgroup;

  int32  MAF_readpoolid;


  // temp std::string to avoid creating throwaway local std::strings
  std::string MAF_tmp_str;
  // temporary vector of std::string to avoid creating throwaway std::vector<std::strings>
  std::vector<std::string> MAF_tmp_strv;

private:
  void deescapeString(std::string & s);

  size_t countReadsBeforeLoad(const std::string & fileName, size_t & maxlinelen);


  void cleanupHeaderData();
  void cleanupReadData();
  void cleanupContigData();

  void checkParseIsInReadGroup(std::string & acttoken);
  void checkParseIsNotInReadGroup(std::string & acttoken);
  void checkParseIsInRead(std::string & acttoken);
  void checkParseIsNotInRead(std::string & acttoken);
  void checkParseIsInContig(std::string & acttoken);

  void parseLineHeaderVersion(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineHeaderReadGroup(std::ifstream & mafin, std::string & acttoken, size_t & linenumber);

  void parseLineAQ(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineRG(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineRS(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineRQ(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineCS(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineCQ(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineRD(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineLR(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineSV(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineTN(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineDI(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineTF(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineTT(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineTS(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineSF(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineBC(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineSL(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineSR(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineQL(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineQR(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineCL(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineCR(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineAO(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineRT(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineST(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineSN(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineMT(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineIB(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineIC(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineIR(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineAT(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineER(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineCO(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineNR(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineLC(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineCT(std::ifstream & mafin, std::string & acttoken, std::string & actline);
  void parseLineEC(std::ifstream & mafin, std::string & acttoken, std::string & actline);

  void parseTagData(std::ifstream & mafin, std::string & acttoken, multitag_t & tag);
  void parseTagDataV2(std::ifstream & mafin, std::string & acttoken, multitag_t & tag);

  void checkReadData();
  void addReadToReadPool();

public:
  MAFParse(ReadPool * rpool,
	   std::list<Contig>  * clist,
	   std::vector<MIRAParameters> * mp);
  ~MAFParse();

  size_t load(const std::string & fileName,
	      const uint8 seqtype,
	      const uint8 loadaction,
	      std::vector<uint32> & lrperseqtype,
	      bool recalcconsensus = true,
	      void (*ccallback)(std::list<Contig> &, ReadPool &) = nullptr, // callback after each contig loaded
	      void (*rcallback)(ReadPool &) = nullptr, // callback after each read loaded
	      bool isVerbose = false);

  friend std::ostream &operator<<(std::ostream &ostr, MAFParse const &i);

  static void parseReadGroup(std::ifstream & mafin,
			     ReadGroupLib::ReadGroupID & rgid,
			     std::vector<ReadGroupLib::ReadGroupID> & externalidmapper,
			     size_t & linenumber);
};





#endif
