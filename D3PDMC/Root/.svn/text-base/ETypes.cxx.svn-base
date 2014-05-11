#include "D3PDMC/ETypes.h"


using namespace ETypes;

//
// Global Settings
//
//static global settings shared by all types
namespace ETypes {
  int SsignificantDigits=2;
  Notations Snotation=kFixed;
  std::string SpmSeparator="+/-";
  bool SroundIntegers=true;
}

int ETypes::SetGlobalSignificantDigits(int significantDigits) {
  return (ETypes::SsignificantDigits = significantDigits);
}

int ETypes::SetGlobalNotation(Notations notation) {
  return (ETypes::Snotation = notation);
}

void ETypes::SetGlobalPMString(std::string pmSeparator) {
  ETypes::SpmSeparator = pmSeparator;
}

int ETypes::SetGlobalPMString(PMSeparators pmSeparator) {
  switch (pmSeparator) {
    case kTextPM:
      SetGlobalPMString(std::string("+/-"));
      break;
    case kTextSpacedPM:
      SetGlobalPMString(std::string(" +/- "));
      break;
    case kLatexPM:
      SetGlobalPMString(std::string("$\\pm$"));
      break;
    case kLatexSpacedPM:
      SetGlobalPMString(std::string("$~\\pm~$"));
      break;
    case kLatexNoMathPM:
      SetGlobalPMString(std::string("\\pm")); //assume already in math mode
      break;
    default:
      SetGlobalPMString(std::string(""));
      break;
  }
  return 1;
}

bool ETypes::SetGlobalRoundIntegers(bool roundIntegers) {
  return (ETypes::SroundIntegers = roundIntegers);
}

//
// Constructors and Co.
//

ETypes::EDouble::EDouble() {
  //init settings
  InitSettings();

  //init values
  _value = 0.0;
  _error = 0.0;
  _precision = 0;
  _roundedError = 0.0;
  _roundedValue = 0.0;
}

ETypes::EDouble::EDouble(double value, double error, int significantDigits) {
  //init settings
  InitSettings();

  //init values
  _value = value;
  _error = error;
  SetSignificantDigits(significantDigits);
  RoundValues();
}

ETypes::EDouble::~EDouble() {
}


//
// Settings
//

void ETypes::EDouble::InitSettings() {
  SetNotation(kN_Default);
  SetSignificantDigits(-1); //use global default one
  SetPMString(kS_Default);
  _valuesChanged = true; //force re-calculating rounded values
  _roundIntegers=-1;//take global one
}

int ETypes::EDouble::SetSignificantDigits(int significantDigits) {
  _valuesChanged = true; //force re-calculation
  return (_significantDigits = significantDigits);
}

int ETypes::EDouble::SetNotation(Notations notation) {
  return (_notation = notation);
}

void EDouble::SetPMString(std::string pmSeparator) {
  _pmSeparator = pmSeparator;
  _pmLSeparator = kS_Custom;
}

int ETypes::EDouble::SetPMString(PMSeparators pmSeparator) {
  _pmLSeparator = pmSeparator;
  switch (pmSeparator) {
    case kTextPM:
      SetPMString(std::string("+/-"));
      break;
    case kTextSpacedPM:
      SetPMString(std::string(" +/- "));
      break;
    case kLatexPM:
      SetPMString(std::string("$\\pm$"));
      break;
    case kLatexSpacedPM:
      SetPMString(std::string("$~\\pm~$"));
      break;
    case kLatexNoMathPM:
      SetPMString(std::string("\\pm")); //assume already in math mode
      break;
    default:
      SetPMString(std::string(""));
      _pmLSeparator = kS_Default;
      break;
  }
  return 1;
}

int ETypes::EDouble::SetRoundIntegers(int roundIntegers) {
  _valuesChanged = true;
  return (_roundIntegers = roundIntegers);
}

void ETypes::EDouble::ForceRecalculation() {
  _valuesChanged = true;
}

void ETypes::EDouble::SetValue(double value) {
  _value = value;
  _valuesChanged = true;
}

void ETypes::EDouble::SetValue(double value, double error) {
  _value = value;
  _error = error;
  _valuesChanged = true;
}

void ETypes::EDouble::SetError(double error) {
  _error = error;
  _valuesChanged = true;
}

double ETypes::EDouble::GetRoundedValue() {
  RoundValues();
  return _roundedValue;
}

double ETypes::EDouble::GetRoundedError() {
  RoundValues();
  return _roundedError;
}

int ETypes::EDouble::GetPrecision() {
  RoundValues();
  return _precision;
}

int ETypes::EDouble::GetSignificantDigits() {
  //  std::cout << "Getting significan digits (global = " << SsignificantDigits << ", local= " << _significantDigits << ")" << std::endl;
  if (_significantDigits < 0) {
    return SsignificantDigits;
  } else {
    return _significantDigits;
  }
}

ETypes::Notations ETypes::EDouble::GetNotation() {
  if (_notation == ETypes::kN_Default) {
    return Snotation;
  } else {
    return _notation;
  }
}

std::string ETypes::EDouble::GetPMString() {
  if (_pmLSeparator == kS_Default) {
    return SpmSeparator;
  } else {
    return _pmSeparator;
  }
}

bool ETypes::EDouble::GetRoundIntegers() {
  if (_roundIntegers >=0 ) { //0=false, 1=true, -1=global(default:true)
    return _roundIntegers;
  } else {
    return SroundIntegers;
  }
}



//
// Internal methods
//
void ETypes::EDouble::RoundValues() {
  if (!_valuesChanged) {
    return;  //avoid double-calculations
  }
  //round values and set precision
  if (_error == 0.0) {
    _roundedError = _error;
    _roundedValue = _value;
    _precision = SsignificantDigits;
    return;
  }
  int precision=0;
  _roundedValue = _value;
  _roundedError = _error;
  int current_significantDigits = GetSignificantDigits();
  //take error in the 10.**(current_significantDigits-1), 10.**(current_significantDigits) interval (extremes included)
  //  std::cout << "Init: " << _roundedValue << " " << _roundedError << std::endl;
  while (_roundedError >= TMath::Power(10., current_significantDigits) ||
         _roundedError <  TMath::Power(10., current_significantDigits-1)) {
    if (_roundedError < TMath::Power(10., current_significantDigits-1)) {
      _roundedError *= 10.;
      precision++;
    } else { // >= TMath::Power(10., current_significantDigits)-1
      _roundedError /= 10.;
      precision--;
    }
    //      std::cout << "Now: Error = "<< _roundedError << ", precision = " << precision << std::endl;
  }

  //convert now also the value to an integer - if needed
  if (_roundIntegers || (precision >=0)) {
    _roundedValue *= TMath::Power(10., precision);
  }

  //round
  _roundedError = (int)round(_roundedError);
  _roundedValue = (int)round(_roundedValue);

  //  std::cout << "Multiplied: val = " << _roundedValue << ", " << _roundedError << std::endl;

  //take the final values back
  if (GetRoundIntegers() || (precision >=0)) {
    _roundedValue /= TMath::Power(10., precision);
  }
  _roundedError /= TMath::Power(10., precision);


  if (precision > 0) {
    _precision = precision;  //digits after '.'
  } else {
    _precision = 0;  //no digits after '.'
  }

  //  std::cout << "Final: val = " << _roundedValue << ", " << _roundedError << ", precision = " << _precision << std::endl;

  _valuesChanged=false;
}

//
// Friendly overloaded operators
//
namespace ETypes {

  std::ostream &operator<<(std::ostream &stream, ETypes::EDouble d) {
    switch (d.GetNotation()) {
      case ETypes::kFixed:
        stream << std::fixed;
        break;
      case ETypes::kScientific:
        stream << std::scientific;
        break;
      default:
        break;
    }
    //std::cout << std::endl << "Writing " << d.GetValue() << "+/-" << d.GetValue() << std::endl;
    stream << std::setprecision(d.GetPrecision()) << d.GetRoundedValue() << d.GetPMString() << d.GetRoundedError();
    return stream;
  }

}
