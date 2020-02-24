"""
Abstract supertype of properties that can be checked.

Every concrete subtype should provide the following functions:
  - `dim(𝑃::Property)::Int`
  - `check(𝑃::Property, X::LazySet; witness::Bool=false)`
"""
abstract type Property end
