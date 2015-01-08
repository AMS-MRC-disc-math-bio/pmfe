from pyparam_base import *
from fractions import Fraction
from decimal import Decimal

try:
    from sage.rings.rational_field import RationalField
    QQ = RationalField()
    IN_SAGE = True
except ImportError:
    IN_SAGE = False

def stringify_fraction(frac):
    if int(frac) == frac:
        return str(int(frac))
    else:
        return str(frac.numerator) + "/" + str(frac.denominator)

def pair_as_fraction(pair):
    return Fraction(pair[0], pair[1])

def vector_as_fractions(vec):
    return tuple(pair_as_fraction(term) for term in vec._as_pairs())

def vector_as_rationals(vec):
    if IN_SAGE:
        return [QQ((frac.numerator, frac.denominator)) for frac in vec.as_fractions()]
    else:
        raise NotImplementedError("This method is only available from a Sage environment.")

def pynum_as_pair(num):
    try:
        if int(num) == num: # if num can be cast cleanly to int
            return (int(num), 1)
    except ValueError: # if int(num) is invalid, that's okay--we'll try the other conversions
        True

    try: # if num is a Sage Rational
        return (int(num.numerator()), int(num.denominator()))
    except TypeError: # otherwise, treat num as a Fraction
        return (num.numerator, num.denominator)

def from_py_scores(cls, vec): # Set a ScoreVector from a Python list
    assert len(vec) == 5
    return cls(
        pynum_as_pair(vec[0]), # Multiloops
        pynum_as_pair(vec[1]), # Unpaired
        pynum_as_pair(vec[2]), # Branches
        pynum_as_pair(vec[3]), # w
        pynum_as_pair(vec[4])  # Energy
    )

def from_py_params(cls, vec): # Set a ParameterVector from a Python list
    assert len(vec) == 4
    return cls(
        pynum_as_pair(vec[0]), # Multiloop penalty
        pynum_as_pair(vec[1]), # Unpaired penalty
        pynum_as_pair(vec[2]), # Branch penalty
        pynum_as_pair(vec[3])  # Dummy scaling
    )

from_py_params.__doc__ = "Assign the parameters from a Python list of integers, Fractions, and/or Sage Rationals. Terms must be given in standard order: [multiloop_penalty, unpaired_penalty, branch_penalty, dummy_scaling]."

def stringify_vector(vec):
    return "[" + ", ".join(stringify_fraction(val) for val in vec.as_fractions()) + "]"

ScoreVector.as_fractions = vector_as_fractions
ScoreVector.as_rationals = vector_as_rationals
ScoreVector.__repr__ = stringify_vector
setattr(ScoreVector, 'from_py_scores', classmethod(from_py_scores))

ParameterVector.as_fractions = vector_as_fractions
ParameterVector.as_rationals = vector_as_rationals
ParameterVector.__repr__ = stringify_vector
setattr(ParameterVector, 'from_py_params', classmethod(from_py_params))
