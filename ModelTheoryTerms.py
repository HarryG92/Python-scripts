# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 11:48:38 2020

Term-building for model theory.

@author: gulli
"""


import copy


class ConstantSymbol:
    def __init__(self, name):
        if not isinstance(name, str):
            raise TypeError("Pass a string")
        else:
            self.string = name


class FunctionSymbol:
    def __init__(self, name, arity, infix=False):
        if not (isinstance(name, str) and isinstance(arity, int)):
            raise TypeError("Pass a string and an int")
        self.arity = arity
        self.string = name
        if infix:
            if arity == 2:
                self.infix = True  # print function in infix notation
            else:
                raise ValueError("can only have infix notation with arity 2")
        else:
            self.infix = False


class RelationSymbol:
    def __init__(self, name, arity, infix=False):
        if not (isinstance(name, str) and isinstance(arity, int)):
            raise TypeError("Pass a string and an int")
        self.arity = arity
        self.string = name
        if infix:
            if arity == 2:
                self.infix = True  # print relation in infix notation
            else:
                raise ValueError("can only have infix notation with arity 2")
        else:
            self.infix = False


# build up terms inductively, given a list of parent terms and a function to
# combine them with
class Term:
    def __init__(self, parent_terms, function_symbol='ATOMIC'):
        self.function_symbol = function_symbol
        if function_symbol == 'ATOMIC':  # if term is a constant or variable
            if isinstance(parent_terms, ConstantSymbol):
                self.parent_terms = [parent_terms]
                self.string = parent_terms.string
                self.free_variables = []  # if constant, no free variables
                self.is_variable = False
            elif (
                  isinstance(parent_terms, list)
                  and len(parent_terms) == 1
                  and isinstance(parent_terms[0], ConstantSymbol)
                  ):  # single parents can be passed in a list or singly
                self.parent_terms = parent_terms
                self.string = parent_terms[0].string
                self.free_variables = []
                self.is_variable = False
            elif isinstance(parent_terms, str):  # if term is a variable
                self.parent_terms = [parent_terms]
                self.string = parent_terms
                self.free_variables = [parent_terms]
                self.is_variable = True
            elif (
                  isinstance(parent_terms, list)
                  and len(parent_terms) == 1
                  and isinstance(parent_terms[0], str)
                  ):  # single parents can be passed in a list or singly
                self.parent_terms = parent_terms
                self.string = parent_terms[0]
                self.free_variables = parent_terms
                self.is_variable = True
            else:
                raise TypeError("""
                                For an atomic term,
                                pass a string or ConstantSymbol"""
                                )
        elif isinstance(function_symbol, FunctionSymbol):
            self.is_variable = False
            if (isinstance(parent_terms, list)
                    and len(parent_terms) == function_symbol.arity):
                self.parent_terms = copy.deepcopy(parent_terms)
                if function_symbol.infix:
                    self.string = (
                                 '('
                                 + self.parent_terms[0].string
                                 + ' '
                                 + function_symbol.string
                                 + ' '
                                 + self.parent_terms[1].string
                                 + ')'
                                 )
                    self.free_variables = (
                            self.parent_terms[0].free_variables
                            + self.parent_terms[1].free_variables
                            )
                else:
                    self.string = "%s(" % function_symbol.string
                    self.free_variables = []
                    for term in range(len(self.parent_terms)):
                        if not isinstance(self.parent_terms[term], Term):
                            raise TypeError("Each parent term must be a Term!")
                        else:
                            self.string += (
                                    "%s, " % self.parent_terms[term].string
                                    )
                            self.free_variables += (
                                    self.parent_terms[term].free_variables)
                    # remove final comma and space; add closing bracket
                    self.string = self.string[0:(len(self.string) - 2)] + ")"
                self.free_variables = list(dict.fromkeys(self.free_variables))
            elif (isinstance(parent_terms, Term)
                    and function_symbol.arity == 1):
                self.parent_terms = [copy.deepcopy(parent_terms)]
                self.string = "%s(%s)" % (
                                        function_symbol.string,
                                        self.parent_terms.string
                                        )
                self.free_variables = self.parent_terms.free_variables
            else:
                raise TypeError("""
                                Must give a list of parent terms matching the
                                arity of the function symbol"""
                                )
        else:
            raise TypeError("Pass a FunctionSymbol to build the term from")

    def number_free_variables(self):
        return len(self.free_variables)

    def substitute(self, old_variable, new_term):
        if (
                isinstance(old_variable, Term)
                and old_variable.is_variable
                and isinstance(new_term, Term)
                ):
            if old_variable.string in self.free_variables:
                if self.function_symbol == 'ATOMIC':
                    new_term_copy = copy.deepcopy(new_term)
                    self.free_variables = new_term_copy.free_variables
                    self.function_symbol = new_term_copy.function_symbol
                    self.string = new_term_copy.string
                    self.parent_terms = new_term_copy.parent_terms
                else:
                    new_parents = []
                    for parent in self.parent_terms:
                        parent_copy = copy.deepcopy(parent)
                        parent_copy.substitute(old_variable, new_term)
                        new_parents.append(parent_copy)
                    Term.__init__(self, new_parents, self.function_symbol)
        else:
            raise TypeError("""
                            pass a variable to be replaced, and a term to
                            replace it with"""
                            )
