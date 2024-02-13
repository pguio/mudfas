/**************************************************************************
 *
 * $Id: Parser.cpp,v 1.35 2011/03/26 12:56:40 patrick Exp $
 *
 * Copyright (c) 2000-2011 Patrick Guio <patrick.guio@gmail.com>
 * All Rights Reserved.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2.  of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include <blitz/array.h>
using namespace blitz;
#include <parser.h>
using namespace parser;

#define ID "$Id: Parser.cpp,v 1.35 2011/03/26 12:56:40 patrick Exp $"

#define COPYRIGHT \
"Copyright (c) 2000-2011 Patrick Guio <patrick.guio@gmail.com>\n\n"\
"This is free software; see the source for copying conditions.  There is NO\n"\
"warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n"

class Foo : public Parser {
  enum { boolean=1, character, chararray, integer, real32, real64,
         intvect, strvect
       };

private:
  bool Boolean;
  char Character;
  string String;
  int Integer;
  float Float;
  double Double;
  TinyVector<int,3> IntVect;
  vector<string> StrVect;
public:
  Foo(int nargs, char *args[]);
  ~Foo()
  {}
  void parseParameters();
  friend ostream& operator<<(ostream &os, const Foo &foo);
};

template<class T>
ostream& operator<<(ostream &os, const  vector<T> &v)
{
  if (!v.empty()) {
    typename vector<T>::const_iterator i = v.begin();
    typename vector<T>::const_iterator end = v.end();
    os << *i;
    for ( ++i; i != end; ++i ) os << ',' << *i;
  }
  return os;
}

Foo::Foo(int nargs, char *args[]) : Parser(nargs, args),
  Boolean(true), Character('\0'),
  String(""), Integer(0), Float(0.0), Double(0.0),
  IntVect(0), StrVect(0)
{
  registerClass("Foo");
  registerPackage(PACKAGE, VERSION "\n" ID "\n", COPYRIGHT);

  using parser::types::boolean;
  using parser::types::real;
  using parser::types::charStr;
  using parser::types::intVect;
  using parser::types::stringVect;

  insertOption(boolean       , "boolean"   , boolean   ,"Parse a bool"               , Any(Boolean));
  insertOptionAlias(boolean  , "-boolean");
  insertOption(character     , "character" , character ,"Parse a char"               , Any(Character));
  insertOptionAlias(character, "-character");

  insertOption(chararray     , "string"    , charStr   ,"Parse a char *"             , Any(String));
  insertOptionAlias(chararray, "-string");

  insertOption(integer       , "integer"   , integer   ,"Parse an int"               , Any(Integer));
  insertOptionAlias(integer  , "--integer");

  insertOption(real32        , "real32"    , real      ,"Parse a float"              , Any(Float));
  insertOption(real64        , "real64"    , real      , "Parse a double"            , Any(Double));
  insertOption(intvect       , "TinyVector", intVect   , "Parse an TinyVector<int,3>", Any(IntVect));
  insertOption(strvect       , "vector"    , stringVect, "Parse a vector<string>"    , Any(StrVect));
}

void Foo::parseParameters()
{
  if (parseOption(boolean,Boolean))
    cout << "Boolean   = " << Boolean << " parsed." << endl;
  if (parseOption(character,Character))
    cout << "Character = " << Character << " parsed." << endl;
  if (parseOption(chararray,String))
    cout << "String    = " << String << " parsed." << endl;
  if (parseOption(integer,Integer))
    cout << "Integer   = " << Integer << " parsed." << endl;
  if (parseOption(real32,Float))
    cout << "Float     = " << Float << " parsed." << endl;
  if (parseOption(real64,Double))
    cout << "Double    = " << Double << " parsed." << endl;
  if (parseOption(intvect,IntVect))
    cout << "IntVect   = " << IntVect << " parsed." << endl;
  if (parseOption(strvect,StrVect)) {
    cout << "StrVect   = " << StrVect << " parsed." << endl;
  }
}

inline
ostream& operator<<(ostream &os, const Foo &foo)
{
  os << "**********" << '\n';
  os << "Foo set up" << '\n';
  os << "**********" << '\n';
  os << "Boolean   = " << foo.Boolean << '\n';
  os << "Character = " << foo.Character << '\n';
  os << "String    = " << foo.String << '\n';
  os << "Integer   = " << foo.Integer << '\n';
  os << "Float     = " << foo.Float << '\n';
  os << "Double    = " << foo.Double << '\n';
  os << "IntVect   = " << foo.IntVect << '\n';
  os << "StrVect   = " << foo.StrVect << '\n';

  return os << "**********";
}

int main(int nargs, char *args[])
{
  try {

#if defined(HAVE_MPI)
    MPI_Init(&nargs, &args);
#endif

    Parser parser(nargs, args);
    parser.registerProgram(args[0]);
    parser.registerPackage(PACKAGE, VERSION, COPYRIGHT);

    Foo foo(nargs, args);

#define PARSE(Fun)   \
if (parser.Fun()) {  \
	foo.Fun();         \
	return 0;          \
}                    \
 
    PARSE(parseHelp)
    PARSE(parseVersion)
    PARSE(parseTemplate)

#undef PARSE

    parser.viewArgs();
    foo.viewArgs();

    foo.parseParameters();

    cout << foo << endl;

#if defined(HAVE_MPI)
    MPI_Finalize();
#endif

    return 0;

  } catch(ClassException& x) {
    cerr << x.what() << endl;
    return EXIT_FAILURE;
  }
}
