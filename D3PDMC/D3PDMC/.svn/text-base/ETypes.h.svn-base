#ifndef __ETYPES_HH__
#define __ETYPES_HH__

////////////////////////////////////////////////////////////////////
//
// Simone Pagan Griso <pagan@fnal.gov>
//
// Collection of type with error.
// Every class implements:
//  - correct formatting of output to streams (with error)
//    given a number of significant number of digits
//  - implement basic operations with error propagation
//  - local and global settings
//
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <string>
#include "TMath.h"
#include <math.h>

namespace ETypes {

  //common types
  //enums
  enum Notations {
    kN_Default, //use global (default) one
    kFixed,
    kScientific,
    kN_Notations //number of different options
  };

  enum PMSeparators {
    kS_Default, //use global (default) one
    kS_Custom,
    kTextPM,
    kTextSpacedPM,
    kLatexPM,
    kLatexSpacedPM,
    kLatexNoMathPM,
    kN_PMSeparators
  };

  //Global settings methods
  int SetGlobalSignificantDigits(int significantDigits);
  int SetGlobalNotation(Notations notation);
  void SetGlobalPMString(std::string pmSeparator); //set custom +/- separator
  int SetGlobalPMString(PMSeparators pmSeparator); //choose among pre-defined ones
  bool SetGlobalRoundIntegers(bool roundIntegers); //round to significan digits if number > 10^sig_digits


  class EDouble {
    protected:

      //internal settings
      int _significantDigits; //number of significant digits on error (value accordingly)
      Notations _notation; //from enum Notations
      std::string _pmSeparator; // '+/-' string
      PMSeparators _pmLSeparator; // logic value
      int _roundIntegers;

      //internal state
      bool _valuesChanged;

      //actual values
      double _value;
      double _error;

      //output
      double _roundedValue;
      double _roundedError;
      int _precision;

    protected:
      //protected methods
      void InitSettings();
      void RoundValues();

    public:
      //constructors, destructor
      EDouble();
      EDouble(double value, double error, int significantDigits=-1);
      ~EDouble();

      //Settings methods
      int SetSignificantDigits(int significantDigits);
      int SetNotation(Notations notation);
      void SetPMString(std::string pmSeparator); //set custom +/- separator
      int SetPMString(PMSeparators pmSeparator); //choose among pre-defined ones
      int SetRoundIntegers(int roundIntegers); //round to significan digits if number > 10^sig_digits. 0=false, 1=true,-1=global
      void ForceRecalculation(); //if you changed global settings and want to re-force calculation


      //accessors
      void SetValue(double value);
      void SetValue(double value, double error);
      void SetError(double error);
      double GetValue() {
        return _value;
      }
      double GetError() {
        return _error;
      }
      double& Value() {
        return _value; //direct access
      }
      double& Error() {
        return _error; //direct access
      }
      double GetRoundedValue();
      double GetRoundedError();
      int GetPrecision();
      int GetSignificantDigits();
      Notations GetNotation();
      std::string GetPMString();
      bool GetRoundIntegers();

      //friendly overloaded operators: +,-,/,*,<<,>>,=
      friend std::ostream &operator<<(std::ostream &stream, EDouble d);
  };


} //namespace ETypes

#endif
