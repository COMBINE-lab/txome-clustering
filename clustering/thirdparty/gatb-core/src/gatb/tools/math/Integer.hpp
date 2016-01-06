/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

/** \file Integer.hpp
 *  \date 01/03/2013
 *  \author edrezen
 *  \brief Entry point class for large integer usage
 */
#ifndef _GATB_CORE_TOOLS_MATH_INTEGER_HPP_
#define _GATB_CORE_TOOLS_MATH_INTEGER_HPP_

/********************************************************************************/
#include <gatb/tools/math/LargeInt.hpp>
#include <gatb/system/api/Exception.hpp>

#include <boost/variant.hpp>

/********************************************************************************/
namespace gatb  {  namespace core  { namespace tools {  namespace math  {
/********************************************************************************/

/** \brief Class for large integers calculus
 *
 * The IntegerTemplate is implemented as a boost variant, which means that it can act like T1, T2, T3 or T4 type
 * according to the configuration.
 *
 * The IntegerTemplate should be specialized with 4 different LargeInt implementation
 * classes.
 *
 * All the methods are implemented through a boost variant visitor.
 *
 *  According to the INTEGER_KIND compilation flag, we define the Integer class
 *  as an alias of one from several possible implementations.
 *
 *  Note that we have 2 possible native implementations (NativeInt64 and NativeInt128)
 *  that rely on native types uint64_t and __uint128_t.
 *
 *  For larger integer, a multi-precision LargeInt is used.
 *
 *  From the user point of view, [s]he has just to include this file and use the Integer
 *  class.
 *
 */
template<typename T1, typename T2, typename T3, typename T4>
class IntegerTemplate
{
private:

    typedef boost::variant<T1,T2,T3,T4> Type;
    Type v;

          Type& operator *()       { return v; }
    const Type& operator *() const { return v; }

public:

    /** Get a number telling which template class is actually used by the variant.
     * \return an integer n, meaning that the Tn template type is used.
     */
    static char& getType()  {  static char instance = 0; return instance; }

    /** Set the type of which template class has to be used by the variant.
     * \param[in] type : integer value
     */
    static void setType (char type)  {  getType() = type; }

    /** Constructor. Note that the type (see getType and setType) has to be first initialized
     * otherwise no instance can be created (exception thrown).
     * \param[in] n : value for initialization of the integer.
     */
    IntegerTemplate (int64_t n=0)
    {
        switch (getType())
        {
        case PREC_1: v = T1(n); break;
        case PREC_2: v = T2(n); break;
        case PREC_3: v = T3(n); break;
        case PREC_4: v = T4(n); break;
        default:  if (getType()<=PREC_4) { v = T4(n);  break; }
                  else { throw system::Exception ("class Integer not initialized"); }
        }
    }

    /** Copy constructor. Relies on the copy constructor of boost variant
     * \param[in] t : the object to be used for initialization
     */
    template<typename T>  explicit IntegerTemplate (const T& t) : v (t)  {}

    /** Affectation operator. Relies on the affectation operator of boost variant
     * \param[in] t : object to be copied
     * \return the current object.
     */
    template<typename T>
    IntegerTemplate& operator=(const T& t)
    {
        v = t;
        return *this;
    }

    /** Get the name of the class used by the variant (ie. one of the Ti template class parameters)
     * \return the class name.
     */
    const char*  getName ()         { return boost::apply_visitor (Integer_name(),  *(*this)); }

    /** Get the size of an instance of the class used by the variant  (ie. one of the Ti template class parameters)
     * \return the size of an object (in bits).
     */
    const size_t getSize ()         { return boost::apply_visitor (Integer_size(),  *(*this)); }

    /** Get the HDF5 type for the the class used by the variant  (ie. one of the Ti template class parameters)
     * \param[in] isCompound : tells whether the type is composed or not
     * \return a HDF5 identifier for the type.
     */
    hid_t hdf5 (bool& isCompound)   { return boost::apply_visitor (Integer_hdf5(isCompound),  *(*this)); }

    /** Operator +
     * \param[in] a : first operand
     * \param[in] b : second operand
     * \return sum of the two operands.
     */
    inline friend IntegerTemplate   operator+  (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_plus(),  *a, *b);  }

    /** Operator -
     * \param[in] a : first operand
     * \param[in] b : second operand
     * \return substraction of the two operands.
     */
    inline friend IntegerTemplate   operator-  (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_minus(), *a, *b);  }

    /** Operator |
     * \param[in] a : first operand
     * \param[in] b : second operand
     * \return 'or' of the two operands.
     */
    inline friend IntegerTemplate   operator|  (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_or(),    *a, *b);  }

    /** Operator ^
     * \param[in] a : first operand
     * \param[in] b : second operand
     * \return 'xor' of the two operands.
     */
    inline friend IntegerTemplate   operator^  (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_xor(),   *a, *b);  }

    /** Operator &
     * \param[in] a : first operand
     * \param[in] b : second operand
     * \return 'and' of the two operands.
     */
    inline friend IntegerTemplate   operator&  (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_and(),   *a, *b);  }

    /** Operator ~
     * \param[in] a : operand
     * \return negation of the operand
     */
    inline friend IntegerTemplate   operator~  (const IntegerTemplate& a)                    {  Integer_compl v; return  boost::apply_visitor (v, *a);  }

    /** Operator ==
     * \param[in] a : first operand
     * \param[in] b : second operand
     * \return equality of the two operands.
     */
    inline friend bool      operator== (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_equals(), *a, *b);  }

    /** Operator !=
     * \param[in] a : first operand
     * \param[in] b : second operand
     * \return inequality of the two operands.
     */
    inline friend bool      operator!= (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  ! (*a==*b);  }

    /** Operator <
     * \param[in] a : first operand
     * \param[in] b : second operand
     * \return '<' of the two operands.
     */
    inline friend bool      operator<  (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_less(),   *a, *b);  }

    /** Operator <=
     * \param[in] a : first operand
     * \param[in] b : second operand
     * \return '<=' of the two operands.
     */
    inline friend bool      operator<= (const IntegerTemplate& a, const IntegerTemplate& b)  {  return  boost::apply_visitor (Integer_lesseq(), *a, *b);  }

    /** Operator *
     * \param[in] a : first operand
     * \param[in] c : second operand
     * \return multiplication of the two operands.
     */
    inline friend IntegerTemplate   operator*  (const IntegerTemplate& a, const int&       c)  {  return  boost::apply_visitor (Integer_mult(c), *a);  }

    /** Operator /
     * \param[in] a : first operand
     * \param[in] c : second operand
     * \return division of the two operands.
     */
    inline friend IntegerTemplate   operator/  (const IntegerTemplate& a, const u_int32_t& c)  {  return  boost::apply_visitor (Integer_div(c),  *a);  }

    /** Operator %
     * \param[in] a : first operand
     * \param[in] c : second operand
     * \return modulo of the two operands.
     */
    inline friend u_int32_t operator%  (const IntegerTemplate& a, const u_int32_t& c)  {  return  boost::apply_visitor (Integer_mod(c),  *a);  }

    /** Operator >>
     * \param[in] a : first operand
     * \param[in] c : second operand
     * \return right shift of the two operands.
     */
    inline friend IntegerTemplate   operator>> (const IntegerTemplate& a, const int& c)  {  return  boost::apply_visitor (Integer_shiftLeft(c),   *a);  }

    /** Operator <<
     * \param[in] a : first operand
     * \param[in] c : second operand
     * \return left shift of the two operands.
     */
    inline friend IntegerTemplate   operator<< (const IntegerTemplate& a, const int& c)  {  return  boost::apply_visitor (Integer_shiftRight(c),  *a);  }

    /** Operator +=
     * \param[in] a : first operand
     * \return addition and affectation.
     */
    void operator+= (const IntegerTemplate& a)  {  boost::apply_visitor (Integer_plusaffect(),  *(*this), *a);  }

    /** Operator ^=
     * \param[in] a : first operand
     * \return xor and affectation.
     */
    void operator^= (const IntegerTemplate& a)  {  boost::apply_visitor (Integer_xoraffect(),   *(*this), *a);  }

    /** Operator[] access the ith nucleotide in the given integer. For instance a[4] get the 5th nucleotide of
     * a kmer encoded as an Integer object.
     * \param[in] idx : index of the nucleotide to be retrieved
     * \return the nucleotide value as follow: A=0, C=1, T=2 and G=3
     */
    u_int8_t  operator[]  (size_t idx) const   { return  boost::apply_visitor (Integer_value_at(idx), *(*this)); }

    /** Get the reverse complement of a kmer encoded as an IntegerTemplate object. Note that the kmer size must be known.
     * \param[in] a : kmer value to be reversed-complemented
     * \param[in] sizeKmer : size of the kmer
     * \return the reverse complement kmer as a IntegerTemplate value
     */
    friend IntegerTemplate revcomp (const IntegerTemplate& a,  size_t sizeKmer)  {  return  boost::apply_visitor (Integer_revomp(sizeKmer),  *a);  }

    /** Get a hash value on 64 bits for a given IntegerTemplate object.
     * \param[in] a : the integer value
     * \param[in] seed : some seed value used for the hash computation.
     * \return the hash value on 64 bits.
     */
    friend u_int64_t hash1        (const IntegerTemplate& a,  u_int64_t seed)  {  return  boost::apply_visitor (Integer_hash1(seed),  *a);          }

    /** Get a hash value on 64 bits for a given IntegerTemplate object.
     * \param[in] a : the integer value
     * \return the hash value on 64 bits.
     */
    friend u_int64_t oahash       (const IntegerTemplate& a)                   {  return  boost::apply_visitor (Integer_oahash(), *a);              }

    /** Get a hash value on 64 bits for a given IntegerTemplate object.
     * Note: although we return 64 bits, only the first 16 bits are set
     * \param[in] a : the integer value
     * \param[in] shift : some value used for the hash computation.
     * \return the hash value on 64 bits.
     */
    friend u_int64_t simplehash16 (const IntegerTemplate& a,  int shift)       {  return  boost::apply_visitor (Integer_simplehash16(shift),  *a);  }

    /** Get an ASCII string representation of a kmer encoded as a IntegerTemplate object
     * \param[in] sizeKmer : size of the kmer
     * \return the ASCII representation of the kmer.
     */
    std::string toString (size_t sizeKmer) const  {  return boost::apply_visitor (Integer_toString(sizeKmer), *(*this)); }

    /** Output stream operator for the IntegerTemplate class
     * \param[in] s : the output stream to be used.
     * \param[in] a : the object to output
     * \return the modified output stream.
     */
    friend std::ostream & operator<<(std::ostream & s, const IntegerTemplate& a)  {  s << *a;  return s;  }

    /** Get the value of the IntegerTemplate object as a U type, U being one of the T1,T2,T3,T4
     * template class parameters. This method can be seen as a converter from the IntegerTemplate class
     * to a specific U type (given as a template parameter of this method).
     * \return  the converted value as a U type.
     */
    template<typename U>
    const U& get ()  const  {  return * boost::get<U>(&v);  }

private:

    struct Integer_name : public boost::static_visitor<const char*>    {
        template<typename T>  const char* operator() (const T& a) const { return a.getName();  }};

    struct Integer_size : public boost::static_visitor<const size_t>    {
        template<typename T>  const size_t operator() (const T& a) const  { return a.getSize();  }};

    struct Integer_plus : public boost::static_visitor<IntegerTemplate>    {
        template<typename T>              IntegerTemplate operator() (const T& a, const T& b) const  { return IntegerTemplate(a + b);  }
        template<typename T, typename U>  IntegerTemplate operator() (const T& a, const U& b) const  { return IntegerTemplate();       }
    };

    struct Integer_minus : public boost::static_visitor<IntegerTemplate>    {
        template<typename T>              IntegerTemplate operator() (const T& a, const T& b) const  { return IntegerTemplate(a - b);  }
        template<typename T, typename U>  IntegerTemplate operator() (const T& a, const U& b) const  { return IntegerTemplate();  }
    };

    struct Integer_or : public boost::static_visitor<IntegerTemplate>    {
        template<typename T>              IntegerTemplate operator() (const T& a, const T& b) const  { return IntegerTemplate(a | b);  }
        template<typename T, typename U>  IntegerTemplate operator() (const T& a, const U& b) const  { return IntegerTemplate();  }
    };

    struct Integer_xor : public boost::static_visitor<IntegerTemplate>    {
        template<typename T>              IntegerTemplate operator() (const T& a, const T& b) const  { return IntegerTemplate(a ^ b);  }
        template<typename T, typename U>  IntegerTemplate operator() (const T& a, const U& b) const  { return IntegerTemplate();  }
    };

    struct Integer_and : public boost::static_visitor<IntegerTemplate>    {
        template<typename T>              IntegerTemplate operator() (const T& a, const T& b) const  { return IntegerTemplate(a & b);  }
        template<typename T, typename U>  IntegerTemplate operator() (const T& a, const U& b) const  { return IntegerTemplate();  }
    };

    struct Integer_less : public boost::static_visitor<bool>    {
        template<typename T>              bool operator() (const T& a, const T& b) const  { return a < b;  }
        template<typename T, typename U>  bool operator() (const T& a, const U& b) const  { return false;  }
    };

    struct Integer_lesseq : public boost::static_visitor<bool>    {
        template<typename T>              bool operator() (const T& a, const T& b) const  { return a <= b;  }
        template<typename T, typename U>  bool operator() (const T& a, const U& b) const  { return false;   }
    };

    struct Integer_equals : public boost::static_visitor<bool>    {
        template<typename T>              bool operator() (const T& a, const T& b) const  { return a == b;  }
        template<typename T, typename U>  bool operator() (const T& a, const U& b) const  { return false;   }
    };

    struct Integer_plusaffect : public boost::static_visitor<>    {
        template<typename T>              void operator() ( T& a, const T& b) const  { a += b;  }
        template<typename T, typename U>  void operator() ( T& a, const U& b) const  {   }
    };

    struct Integer_xoraffect : public boost::static_visitor<>    {
        template<typename T>              void operator() ( T& a, const T& b) const  { a ^= b;  }
        template<typename T, typename U>  void operator() ( T& a, const U& b) const  {   }
    };

    struct Integer_compl : public boost::static_visitor<IntegerTemplate>    {
        template<typename T>  IntegerTemplate operator() (const T& a)  { return IntegerTemplate(~a);  }};

    template<typename Result, typename Arg>
    struct Visitor : public boost::static_visitor<Result>
    {
        Visitor (Arg a=Arg()) : arg(a) {}
        Arg arg;
    };



	
    struct Integer_hdf5 : public Visitor<hid_t,bool&>   {
        Integer_hdf5 (bool& c) : Visitor<IntegerTemplate,bool&>(c) {}
        template<typename T>  hid_t operator() (const T& a)  { return a.hdf5 (this->arg);  }};

    struct Integer_mult : public Visitor<IntegerTemplate,const int>    {
        Integer_mult (const int& c) : Visitor<IntegerTemplate,const int>(c) {}
        template<typename T>  IntegerTemplate operator() (const T& a) const  { return IntegerTemplate(a*this->arg);  }};

    struct Integer_div : public Visitor<IntegerTemplate,const u_int32_t>    {
        Integer_div (const u_int32_t& c) : Visitor<IntegerTemplate,const u_int32_t>(c) {}
        template<typename T>  IntegerTemplate operator() (const T& a) const  { return IntegerTemplate(a/this->arg);  }};

    struct Integer_mod : public Visitor<u_int32_t,const u_int32_t>    {
        Integer_mod (const u_int32_t& c) : Visitor<u_int32_t,const u_int32_t>(c) {}
        template<typename T>  u_int32_t operator() (const T& a) const  { return (a%this->arg);  }};

    struct Integer_shiftLeft : public Visitor<IntegerTemplate,const int>    {
        Integer_shiftLeft (const int& c) : Visitor<IntegerTemplate,const int>(c) {}
        template<typename T>  IntegerTemplate operator() (const T& a) const  { return IntegerTemplate (a >> this->arg);  }};

    struct Integer_shiftRight : public Visitor<IntegerTemplate,const int>    {
        Integer_shiftRight (const int& c) : Visitor<IntegerTemplate,const int>(c) {}
        template<typename T>  IntegerTemplate operator() (const T& a) const  { return IntegerTemplate (a << this->arg);  }};

    struct Integer_revomp : public Visitor<IntegerTemplate,size_t>    {
        Integer_revomp (const size_t& c) : Visitor<IntegerTemplate,size_t>(c) {}
        template<typename T>  IntegerTemplate operator() (const T& a) const  { return IntegerTemplate (revcomp(a,this->arg));  }};

    struct Integer_hash1 : public Visitor<u_int64_t,u_int64_t>    {
        Integer_hash1 (const u_int64_t& c) : Visitor<u_int64_t,u_int64_t>(c) {}
        template<typename T>  u_int64_t operator() (const T& a) const  { return (hash1(a,this->arg));  }};

    struct Integer_oahash : public boost::static_visitor<u_int64_t>    {
        template<typename T>  u_int64_t operator() (const T& a) const  { return (oahash(a));  }};

	
    struct Integer_simplehash16 : public Visitor<u_int64_t,int>    {
        Integer_simplehash16 (const int& c) : Visitor<u_int64_t,int>(c) {}
        template<typename T>  u_int64_t operator() (const T& a) const  { return (simplehash16(a,this->arg));  }};

    struct Integer_value_at : public Visitor<u_int8_t,size_t>   {
        Integer_value_at (size_t idx) : Visitor<u_int8_t,size_t>(idx) {}
        template<typename T>  u_int8_t operator() (const T& a) const { return a[this->arg];  }};

    struct Integer_toString : public Visitor<std::string,size_t>   {
        Integer_toString (size_t c) : Visitor<std::string,size_t>(c) {}
        template<typename T>  std::string operator() (const T& a) const  { return a.toString(this->arg);  }};
};

/********************************************************************************/

#define INTEGER_TYPES   LargeInt<PREC_1>,LargeInt<PREC_2>,LargeInt<PREC_3>,LargeInt<PREC_4>

typedef IntegerTemplate <INTEGER_TYPES> Integer;

/********************************************************************************/
}}}};
/********************************************************************************/

#endif /* _GATB_CORE_TOOLS_MATH_INTEGER_HPP_ */
