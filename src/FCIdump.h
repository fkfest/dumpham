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
#ifndef FCIDUMP_H
#define FCIDUMP_H
#ifdef __cplusplus
#include <string>
#include <vector>
#include <fstream>
/**
 @if MOLPRO
 @page FCIdump FCIdump access through C++, C and Fortran
 @else
 @mainpage FCIdump access through C++, C and Fortran
 @endif

 @section Introduction
 This is a C++ class, together with non-object-oriented C and Fortran 90 bindings, that implement access
 to an FCIDUMP file as produced by the Full CI code (Computer Physics Communications
Volume 54, Issue 1, April 1989, Pages 75–83,
<a href="http://dx.doi.org/10.1016/0010-4655(89)90033-7">
doi:10.1016/0010-4655(89)90033-7</a>;
<a href="http://bitbucket.org:pjknowles/fci">
http://bitbucket.org/pjknowles/fci</a>),
as well as <a href="http://www.molpro.net/">Molpro</a> via its fci;dump facility

@author Peter Knowles

 */

/*!
 * \brief C++ class that provides access to FCIdump files
 */
class FCIdump
{
public:
  /*!
     * \brief Construct an empty FCIdump object
     */
  FCIdump();

  /*!
     * \brief Construct FCIdump object
     * \param filename The file containing the FCIDUMP data
     */
  FCIdump(const std::string filename);

  ~FCIdump();

  FCIdump(FCIdump&&);

  /*!
     * \brief Construct FCIdump object from bytestream
     * \param bytestream The serialised representation
     */
  FCIdump(const std::vector<char> bytestream);

  /*!
     * \brief Obtain an integer namelist parameter from the FCIDUMP data.
     * \param key The name of the parameter
     * \param def Default value if the parameter is not found.
     * \return  The result as a vector of integers.
     */
  std::vector<int> parameter(const std::string &key, const std::vector<int>& def=std::vector<int>(1,0)) const;

  /*!
     * \brief Obtain a real namelist parameter from the FCIDUMP data.
     * \param key The name of the parameter
     * \param def Default value if the parameter is not found.
     * \return  The result as a vector of integers.
     */
  std::vector<double> parameter(const std::string &key, const std::vector<double> &def) const;

  /*!
     * \brief Obtain a string namelist parameter from the FCIDUMP data.
     * \param key The name of the parameter
     * \param def Default value if the parameter is not found.
     * \return  The result as a vector of integers.
     */
  std::vector<std::string> parameter(const std::string& key, const std::vector<std::string>& def) const;
  /*!
   * \brief Add a parameter with array values
   * \param key key
   * \param values values
   */
  void addParameter(const std::string& key, const std::vector<std::string>& values);

  /*!
   * \brief Add a parameter with array values
   * \param key key
   * \param values values
   */
  void addParameter(const std::string& key, const std::vector<int>& values);

  /*!
   * \brief Add a parameter with array values
   * \param key key
   * \param values values
   */
  void addParameter(const std::string& key, const std::vector<double>& values);

  /*!
   * \brief Add a parameter with a scalar value
   * \param key key
   * \param value value
   */
  void addParameter(const std::string& key, const std::string& value);

  /*!
   * \brief Add a parameter with a scalar value
   * \param key key
   * \param value value
   */
  void addParameter(const std::string& key, const int& value);

  /*!
   * \brief Add a parameter with a scalar value
   * \param key key
   * \param value value
   */
  void addParameter(const std::string& key, const double& value);

  /*!
   * \brief Modify a parameter with array values
   * \param key key
   * \param values values
   */
  void modifyParameter(const std::string& key, const std::vector<std::string>& values);

  /*!
   * \brief Modify a parameter with array values
   * \param key key
   * \param values values
   */
  void modifyParameter(const std::string& key, const std::vector<int>& values);

  /*!
   * \brief Modify a parameter with array values
   * \param key key
   * \param values values
   */
  void modifyParameter(const std::string& key, const std::vector<double>& values);

  /*!
   * \brief Modify a parameter with a scalar value
   * \param key key
   * \param value value
   */
  void modifyParameter(const std::string& key, const std::string& value);

  /*!
   * \brief Modify a parameter with a scalar value
   * \param key key
   * \param value value
   */
  void modifyParameter(const std::string& key, const int& value);

  /*!
   * \brief Modify a parameter with a scalar value
   * \param key key
   * \param value value
   */
  void modifyParameter(const std::string& key, const double& value);
  
   /*!
     * \brief The file containing the FCIDUMP data
     */
  std::string fileName() const;

  /*!
   * \brief The possible external file formats
    */
  typedef enum {
    FileFormatted ///< formatted ASCII file
  } fileType;

  /*!
   * \brief Serialise the object to a stream of bytes
   * \param integrals If true, write out the integrals as well as the headers
   * \return the serialised representation of the object
   */
  std::vector<char> bytestream(bool integrals=true);

  /*!
   * \brief Write the object to an external file
   * \param filename The relative or absolute path name of the file
   * \param type The desired format of the file
   * \param integrals If true, write out the integrals; otherwise leave the file open and positioned ready to write integrals later
   * \return true if OK, false if not
   */
  bool write(std::string filename, fileType type=FileFormatted, bool integrals=true);

  /*!
   * \brief writeIntegral Write an integral to the output stream. write() must already have been called.
   * \param i Orbital index
   * \param j Orbital index
   * \param k Orbital index
   * \param l Orbital index
   * \param value The integral
   */
  void writeIntegral(int i, int j, int k, int l, double value) const;

  /*!
   * \brief Indicator of the type of integral record (core, 1-electron, 2-electron integrals; end of record; end of file)
    */
  typedef enum {
    I0, ///< scalar part of energy, ie nuclear repulsion plus any included core
    I1a, ///< 1-electron integrals, alpha spin
    I1b, ///< 1-electron integrals, beta spin
    I2aa, ///< 2-electron integrals, alpha-alpha spin
    I2ab, ///< 2-electron integrals, alpha-beta spin
    I2bb, ///< 2-electron integrals, beta-beta spin
    endOfFile, ///< End of file
    endOfRecord ///< Separator between different spin cases
  } integralType;

  /*!
   * \brief Position the file so that the next call to nextIntegral will deliver the first integral
   */
  void rewind() const;

  /*!
   * \brief Read the next integral from the file
   * \param i orbital label (zero indicates not 1-electron or 2-electron)
   * \param j orbital label
   * \param k orbital label(zero indicates not 2-electron)
   * \param l orbital label
   * \param value numerical value of the integral
   * \return indicator of the type of entry (core, 1-electron, 2-electron integrals; end of record; end of file)
   */
  integralType nextIntegral(int& i, int& j, int& k, int& l, double& value) const;
 std::string data() const { return namelistData;}


private:
  std::string namelistData;
  std::string _fileName;
  mutable std::ifstream stream;
  mutable std::ofstream outputStream;
  mutable bool uhf;
  mutable std::vector<integralType> states;
  mutable std::vector<integralType>::const_iterator currentState;
};
#endif

// C binding
#ifdef __cplusplus
extern "C" {
#endif
/*!
 * \brief C binding of FCIdump: initialise access to an FCIDUMP
 * \param filename The file containing the FCIDUMP data
 */
void FCIdumpInitialise(char* filename);
/*!
 * \brief C binding of FCIdump:  Obtain a string namelist parameter from the FCIDUMP data.
 * \param key The name of the parameter
 * \param value  The result as a char* (arrays not supported). If not present in the FCIDUMP, value is not overwritten, i.e. value serves as a default
 */
void FCIdumpParameterS(char* key, char* value);
/*!
 * \brief C binding of FCIdump:  Obtain an integer namelist parameter from the FCIDUMP data.
 * \param key The name of the parameter
 * \param values  The result as a vector of integers. Any elements not present in the FCIDUMP are not overwritten, i.e. values serves as a list of defaults
 * \param n The length of values
 */
void FCIdumpParameterI(char* key, int* values, int n);
/*!
 * \brief C binding of FCIdump:  Obtain a floating-point namelist parameter from the FCIDUMP data.
 * \param key The name of the parameter
 * \param values  The result as a vector of doubles. Any elements not present in the FCIDUMP are not overwritten, i.e. values serves as a list of defaults
 * \param n The length of values
 */
void FCIdumpParameterF(char* key, double* values, int n);
/*!
 * \brief C binding of FCIdump: Position the file so that the next call to FCIdumpNextIntegral will deliver the first integral
 */
void FCIdumpRewind();
  /*!
 * \brief C binding of FCIdump: Read the next integral from the file
   * \param i orbital label (zero indicates not 1-electron or 2-electron)
   * \param j orbital label
   * \param k orbital label(zero indicates not 2-electron)
   * \param l orbital label
   * \param value numerical value of the integral
   * \return indicator of the type of entry (core, 1-electron, 2-electron integrals; end of record; end of file)
   */
int FCIdumpNextIntegral(int* i, int* j, int* k, int* l, double* value);
/*!
 * \brief C binding of FCIdump: add a parameter
 * \param key key
 * \param value value. Note that through this interface only a single string, not an array, can be given
 */
void FCIdumpAddParameterS(char* key, char* value);
/*!
 * \brief C binding of FCIdump: add a parameter
 * \param key key
 * \param values values
 * \param n length of values
 */
void FCIdumpAddParameterI(char* key, int values[], int n);
/*!
 * \brief C binding of FCIdump: add a parameter
 * \param key key
 * \param values values
 * \param n length of values
 */
void FCIdumpAddParameterF(char* key, double values[],int n);
  /*!
   * \brief C binding of FCIdump: write the data to an external file
   * \param filename The relative or absolute path name of the file
   * \param type The desired format of the file
   * \return 1 if OK, 0 if not
   */
int FCIdumpWrite(char* filename, int type);
#ifdef __cplusplus
}
#endif

#endif // FCIDUMP_H
