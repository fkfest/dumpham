/*
Copyright (c) 2015, Peter J Knowles.
Copyright (c) 2018, Daniel Kats.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the Peter J Knowles nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL Peter J Knowles BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifdef MOLPRO
#include "common/molpro_config.h"
#else
#define _GNU_SOURCE
#endif
#include "FCIdump.h"
#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <stdexcept>
#define xout std::cout

FCIdump::FCIdump()
{
  _fileName="";
  namelistData=",";
}

FCIdump::FCIdump(const std::string filename)
{
  _fileName = filename;
  std::ifstream s;
  s.open(_fileName.c_str());
  if ( (s.rdstate() & std::ifstream::failbit ) != 0 ) {
    std::cerr << "Error opening " << _fileName <<std::endl;
    throw std::runtime_error("FCIDUMP::parameter file missing");
  }
  // cache the namelist data
  std::string ss;
  s >> ss ; namelistData=","; // forget the first word of namelist
  while (s >> ss && ss != "&END" && ss != "/")
    namelistData.append(ss);
  namelistData.append(",DUMMY_KEY=,"); // dummy entry at end to simplify parsing
//  xout <<"namelistData=" <<namelistData <<std::endl;
}

FCIdump::FCIdump(const std::vector<char> bytestream)
{

}

FCIdump::~FCIdump()
{
  if (outputStream.is_open()) outputStream.close();
}


std::string FCIdump::fileName() const
{
  return _fileName;
}



std::vector<int> FCIdump::parameter(const std::string& key, const std::vector<int> &def) const { // dirty sucking in from FCIDUMP namelist
  std::vector<int> answer;
  std::vector<std::string> strings = parameter(key,std::vector<std::string>(1," "));
  if (strings == std::vector<std::string>(1," ")) return def;
  for (std::vector<std::string>::const_iterator s1=strings.begin(); s1!= strings.end(); s1++) {
    int i;
    std::istringstream b(*s1);
    b >>i;
    answer.push_back(i);
  }
  return answer;
}


std::vector<double> FCIdump::parameter(const std::string& key, const std::vector<double>& def) const { // dirty sucking in from FCIDUMP namelist
  std::vector<double> answer;
  std::vector<std::string> strings = parameter(key,std::vector<std::string>(1," "));
  if (strings == std::vector<std::string>(1," ")) return def;
  for (std::vector<std::string>::const_iterator s1=strings.begin(); s1!= strings.end(); s1++) {
    double i;
    std::istringstream b(*s1);
    b >>i;
    answer.push_back(i);
  }
  return answer;
}

std::vector<std::string> FCIdump::parameter(const std::string& key, const std::vector<std::string>& def) const { // dirty sucking in from FCIDUMP namelist
  std::vector<std::string> answer;
  size_t pos = namelistData.find(","+key+"=");
  if (pos == std::string::npos) return def;
  pos=namelistData.find("=",pos)+1;
  size_t pose=namelistData.find("=",pos);
  pose=namelistData.find_last_of(",",pose)+1;
  for (size_t posNext=pos; posNext < pose; pos=posNext) {
    posNext=namelistData.find(",",pos)+1;
    answer.push_back(namelistData.substr(pos,posNext-pos-1));
  }
  return answer;
}

void FCIdump::addParameter(const std::string& key, const std::vector<std::string>& values)
{
//  xout << "FCIdump::addParameter namelistData originally "<<namelistData<<std::endl;
  namelistData.erase(0,1);
  for (std::vector<std::string>::const_reverse_iterator s=values.rbegin(); s != values.rend(); s++)
    namelistData.insert(0,(*s)+",");
  namelistData.insert(0,","+key+"=");
//  xout << "FCIdump::addParameter namelistData set to "<<namelistData<<std::endl;
}
void FCIdump::addParameter(const std::string& key, const std::vector<int>& values)
{
  std::vector<std::string> valuess;
  for (std::vector<int>::const_iterator v=values.begin(); v!=values.end(); v++) {
    std::ostringstream ss;
    ss << *v;
    valuess.push_back(ss.str());
  }
  addParameter(key,valuess);
}
void FCIdump::addParameter(const std::string& key, const std::vector<double>& values)
{
  std::vector<std::string> valuess;
  for (std::vector<double>::const_iterator v=values.begin(); v!=values.end(); v++) {
    std::ostringstream ss;
    ss << *v;
    valuess.push_back(ss.str());
  }
  addParameter(key,valuess);
}
void FCIdump::addParameter(const std::string &key, const std::string &value)
{
  addParameter(key,std::vector<std::string>(1,value));
}
void FCIdump::addParameter(const std::string &key, const int &value)
{
  addParameter(key,std::vector<int>(1,value));
}
void FCIdump::addParameter(const std::string &key, const double &value)
{
  addParameter(key,std::vector<double>(1,value));
}

void FCIdump::modifyParameter(const std::string& key, const std::vector< std::string >& values)
{
//   xout << "FCIdump::addParameter namelistData originally "<< namelistData << std::endl;
  size_t pos = namelistData.find(","+key+"=");
  if (pos == std::string::npos) return addParameter(key,values);
  size_t posb = pos;
  pos = namelistData.find("=",pos)+1;
  size_t pose = namelistData.find("=",pos);
  pose = namelistData.find_last_of(",",pose)+1;
  namelistData.erase(posb,pose-posb);
  for (std::vector<std::string>::const_reverse_iterator s=values.rbegin(); s != values.rend(); s++)
    namelistData.insert(posb,(*s)+",");
  namelistData.insert(posb,","+key+"=");
//   xout << "FCIdump::addParameter namelistData set to "<< namelistData << std::endl;
}
void FCIdump::modifyParameter(const std::string& key, const std::vector<int>& values)
{
  std::vector<std::string> valuess;
  for (std::vector<int>::const_iterator v=values.begin(); v!=values.end(); v++) {
    std::ostringstream ss;
    ss << *v;
    valuess.push_back(ss.str());
  }
  modifyParameter(key,valuess);
}
void FCIdump::modifyParameter(const std::string& key, const std::vector<double>& values)
{
  std::vector<std::string> valuess;
  for (std::vector<double>::const_iterator v=values.begin(); v!=values.end(); v++) {
    std::ostringstream ss;
    ss << *v;
    valuess.push_back(ss.str());
  }
  modifyParameter(key,valuess);
}
void FCIdump::modifyParameter(const std::string &key, const std::string &value)
{
  modifyParameter(key,std::vector<std::string>(1,value));
}
void FCIdump::modifyParameter(const std::string &key, const int &value)
{
  modifyParameter(key,std::vector<int>(1,value));
}
void FCIdump::modifyParameter(const std::string &key, const double &value)
{
  modifyParameter(key,std::vector<double>(1,value));
}

void FCIdump::rewind() const
{
  stream.open(_fileName.c_str());
  std::string ss;
  while (stream >> ss && ss != "&END" && ss != "/") ;
  uhf = parameter("IUHF").at(0) != 0;
  states.clear();
  states.push_back(I2aa);
  if (uhf) states.push_back(I2bb);
  if (uhf) states.push_back(I2ab);
  states.push_back(I1a);
  if (uhf) states.push_back(I1b);
  states.push_back(I0);
  currentState = states.begin();
}

#include <iomanip>
bool FCIdump::write(std::string filename, fileType type, bool integrals)
{
  outputStream.open(filename.c_str());
  if ( (outputStream.rdstate() & std::ifstream::failbit ) != 0 ) {
    xout << "FCIdump::write failed to open "<<filename<<std::endl;
    outputStream.close();
    return false;
  }
  outputStream << "&FCI" << std::endl;
  size_t pos=0;
  while ( (pos = namelistData.find('=',pos)) != std::string::npos){
    std::string key=namelistData.substr(0,pos);
    size_t pose=namelistData.find("=",pos+1);
    pose=namelistData.find_last_of(",",pose)+1;
    while (namelistData.substr(pose-2,1) == ",") pose--;
    pos=namelistData.find_last_of(",",pos)+1;
    if (namelistData.substr(pos,pose-pos).find("DUMMY_KEY") == std::string::npos)
      outputStream << " " << namelistData.substr(pos,pose-pos) << std::endl;
    pos=pose;
  }
  outputStream<< "/"<<std::endl;
  outputStream<<std::scientific<<std::setprecision(15);
  rewind();
  if (integrals) {
  int i,j,k,l;
  double value;
  while (nextIntegral(i,j,k,l,value) != endOfFile)
    writeIntegral(i,j,k,l,value);
  }
  return true;
}

#include <stdint.h>
std::vector<char> FCIdump::bytestream(bool integrals)
{
  std::vector<char> result;
  for (std::string::const_iterator s=namelistData.begin(); s!=namelistData.end(); s++)
    result.push_back(*s);
  if (integrals) {
    rewind();
    int i,j,k,l;
    union {
      struct {
        int16_t labels[4];
        double value;
      } s;
      char buf[16];
    } u;
    while (nextIntegral(i,j,k,l,u.s.value)!=endOfFile) {
      u.s.labels[0]=i;
      u.s.labels[1]=j;
      u.s.labels[2]=k;
      u.s.labels[3]=l;
    }
    for(i=0; i<16; i++) result.push_back(u.buf[i]);
  }
  return result;
}

void FCIdump::writeIntegral(int i, int j, int k, int l, double value) const
{
//     outputStream << value << " " << i << " " << j << " " << k << " " << l << std::endl;
    outputStream << std::setw(23) << value << " ";
    outputStream << std::setw(3) << i<<" ";
    outputStream << std::setw(3) << j<<" ";
    outputStream << std::setw(3) << k<<" ";
    outputStream << std::setw(3) << l << std::endl;
}

FCIdump::integralType FCIdump::nextIntegral(int &i, int &j, int &k, int &l, double &value) const
{
  integralType result = *currentState;
  if (stream >> value) {
    stream >> i; stream >> j; stream >> k; stream >> l;
  }
  else {
    stream.close();
    return endOfFile;
  }
  // following is tricky stuff reflecting historical structure of UHF and RHF FCIdump files
  if (i == 0) {
//    xout << "zero read uhf="<<uhf<<", *currentState="<<*currentState<<std::endl;
    if (uhf && *currentState != I0) {
      result = endOfRecord;
      currentState++;
//      xout << "end of Record signalled; *currentState="<<*currentState<<std::endl;
    }
    else {
//      xout << "real scalar signalled"<<std::endl;
      result = I0;
    }
  }
  else if (k == 0 && (*currentState != FCIdump::I1a && *currentState != FCIdump::I1b)) {
//    xout << "special state switch to "<<*(currentState+1)<<std::endl;
    result=(*currentState==FCIdump::I2aa)?FCIdump::I1a:FCIdump::endOfRecord; ++currentState;
  }
//  xout << "FCIdump::nextIntegral i,j,k,l,value "<<i<<","<<j<<","<<k<<","<<l<<","<<value<<", result="<<result<<std::endl;
  return result;
}

FCIdump* dump;

void FCIdumpInitialise(char* filename)
{
//  printf("Initialise FCIDUMP from file %s\n",filename);
  dump = new FCIdump(std::string(filename));
}

#include <string.h>
void FCIdumpParameterS(char* key, char *value)
{
  std::vector<std::string> vals(1);
  vals[0]=value;
  vals = dump->parameter(std::string(key), vals);
  if (vals.size()>0) value=strdup(vals[0].c_str());
}

void FCIdumpParameterI(char* key, int* values, int n)
{
  std::vector<int> vals(n);
  for (size_t i=0; i<(size_t)n; i++) vals[i]=values[i];
  vals = dump->parameter(std::string(key), vals);
  for (size_t i=0; i<vals.size() && i < (size_t)n; i++)
    values[i]=vals[i];
}

void FCIdumpParameterF(char* key, double* values, int n)
{
  std::vector<double> vals(n);
  for (size_t i=0; i<(size_t)n; i++) vals[i]=values[i];
  vals = dump->parameter(std::string(key), vals);
  for (size_t i=0; i<vals.size() && i < (size_t)n; i++)
    values[i]=vals[i];
}

void FCIdumpRewind()
{
  dump->rewind();
}

int FCIdumpNextIntegral(int* i, int* j, int* k, int* l, double* value)
{
  int ii, jj, kk, ll;
  double vv;
  FCIdump::integralType type = dump->nextIntegral(ii, jj, kk, ll, vv);
  *i = ii; *j = jj; *k = kk; *l = ll; *value = vv;
  return (int) type;
}

void FCIdumpAddParameterS(char *key, char *value)
{
  std::string keys(key);
  std::vector<std::string> valuess(1,value);
  dump->addParameter(keys,valuess);
}

void FCIdumpAddParameterI(char *key, int values[], int n)
{
  std::string keys(key);
  std::vector<int> valuess;
  for (size_t i=0; i<(size_t)n; i++)
    valuess.push_back(values[i]);
  dump->addParameter(keys,valuess);
}

void FCIdumpAddParameterF(char *key, double values[], int n)
{
  std::string keys(key);
  std::vector<double> valuess;
  for (size_t i=0; i<(size_t)n; i++)
    valuess.push_back(values[i]);
  dump->addParameter(keys,valuess);
}

int FCIdumpWrite(char* filename, int type)
{
  return dump->write(std::string(filename),(FCIdump::fileType) type) ? 1 : 0 ;
}
