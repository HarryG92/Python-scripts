# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 11:47:20 2020

@author: gulli
"""

import ModelTheoryTerms as Terms

import copy


class VariableClashError(Exception):
    def __init__(self, variable_list):
        self.variable_list = variable_list


class InductionMethod:
    def __init__(self):
        self.method = ''

    def string(self, parent_formulae):
        return ''

    def free_variables(self, parent_formulae):
        return []

    def bound_variables(self, parent_formulae):
        return []


class AtomicFormulaBuilder(InductionMethod):
    def __init__(self, relation='='):
        self.relation = relation
        self.method = 'ATOMIC'

    def string(self, parent_formulae):
        if (
                isinstance(self.relation, Terms.RelationSymbol)
                and not self.relation.infix
                ):
            if (
                    isinstance(parent_formulae, list)
                    and len(parent_formulae) == self.relation.arity
                    ):
                formula_string = self.relation.string + '('
                for formula in parent_formulae:
                    if isinstance(formula, Terms.Term):
                        formula_string = (
                                          formula_string
                                          + formula.string
                                          + ', '
                                          )
                    else:
                        raise TypeError("""
                                        all parent formulae must be
                                        Terms"""
                                        )
                formula_string = (
                                  formula_string[0:(len(formula_string) - 2)]
                                  + ')'
                                  )
                return formula_string
        elif self.relation == '=' or (
                                      isinstance(
                                                 self.relation,
                                                 Terms.RelationSymbol
                                                 )
                                      and self.relation.infix
                                      ):
            if self.relation == '=':
                relation_symbol = ' = '
            else:
                relation_symbol = ' ' + self.relation.string + ' '
            if (
                    isinstance(parent_formulae, list)
                    and len(parent_formulae) == 2
                    ):
                if (
                        isinstance(parent_formulae[0], Terms.Term)
                        and isinstance(parent_formulae[1], Terms.Term)
                        ):
                    return (
                            '('
                            + parent_formulae[0].string
                            + relation_symbol
                            + parent_formulae[1].string
                            + ')'
                            )
                else:
                    raise TypeError("must pass Terms")
            else:
                raise TypeError("must pass a list of two Terms")
        else:
            raise TypeError("pass a relation symbol (or leave blank for =)")

    def free_variables(self, parent_formulae):
        if isinstance(parent_formulae, list):
            free_variables = []
            for formula in parent_formulae:
                if isinstance(formula, Terms.Term):
                    free_variables += formula.free_variables
                else:
                    raise TypeError("""
                                    all parent formulae must be
                                    Terms"""
                                    )
            free_variables = list(dict.fromkeys(free_variables))
            return free_variables
        else:
            raise TypeError("pass a list of parent terms")


class BooleanFormulaBuilder(InductionMethod):
    def __init__(self, name, symbol):
        if isinstance(symbol, str):
            self.symbol = symbol
        else:
            raise TypeError("symbol must be a string")
        if isinstance(name, str):
            self.method = name
        else:
            raise TypeError("name must be a string")

    def string(self, parent_formulae):
        if isinstance(parent_formulae, list):
            if len(parent_formulae) == 2:
                variable_clashes = []
                for variable in parent_formulae[0].bound_variables:
                    if variable in parent_formulae[1].free_variables:
                        variable_clashes.append(variable)
                for variable in parent_formulae[1].bound_variables:
                    if variable in parent_formulae[0].free_variables:
                        variable_clashes.append(variable)
                if variable_clashes == []:
                    return (
                            '('
                            + parent_formulae[0].string
                            + ' '
                            + self.symbol
                            + ' '
                            + parent_formulae[1].string
                            + ')'
                            )
                else:
                    raise VariableClashError(variable_clashes)
            else:
                raise ValueError("pass a list of TWO formulae to be connected")
        else:
            raise TypeError("pass a list of two formulae")

    def free_variables(self, parent_formulae):
        if isinstance(parent_formulae, list):
            if len(parent_formulae) == 2:
                return list(dict.fromkeys(
                                          parent_formulae[0].free_variables
                                          + parent_formulae[1].free_variables
                                          ))
            else:
                raise ValueError("pass a list of TWO formulae to be connected")
        else:
            raise TypeError("pass a list of two formulae")

    def bound_variables(self, parent_formulae):
        if isinstance(parent_formulae, list):
            if len(parent_formulae) == 2:
                return list(dict.fromkeys(
                                          parent_formulae[0].bound_variables
                                          + parent_formulae[1].bound_variables
                                          ))
            else:
                raise ValueError("pass a list of TWO formulae to be connected")
        else:
            raise TypeError("pass a list of two formulae")


AND = BooleanFormulaBuilder('AND', 'AND')
OR = BooleanFormulaBuilder('OR', 'OR')
IMPLIES = BooleanFormulaBuilder('IMPLIES', '->')
IFF = BooleanFormulaBuilder('IFF', '<->')


class NotFormulaBuilder(InductionMethod):
    def __init__(self):
        self.method = 'NOT'

    def string(self, parent_formulae):
        if isinstance(parent_formulae, list):
            if len(parent_formulae) == 1:
                return '(' + '¬' + parent_formulae[0].string + ')'
            else:
                raise ValueError("pass a list of ONE formula to be negated")
        elif isinstance(parent_formulae, Formula):
            return '(' + '¬' + parent_formulae.string + ')'
        else:
            raise TypeError("pass a list of one formula")

    def free_variables(self, parent_formulae):
        if isinstance(parent_formulae, list):
            if len(parent_formulae) == 1:
                return parent_formulae[0].free_variables
            else:
                raise ValueError("pass a list of ONE formula to be negated")
        elif isinstance(parent_formulae, Formula):
            return parent_formulae.free_variables
        else:
            raise TypeError("pass a list of one formula")

    def bound_variables(self, parent_formulae):
        if isinstance(parent_formulae, list):
            if len(parent_formulae) == 1:
                return parent_formulae[0].bound_variables
            else:
                raise ValueError("pass a list of ONE formula to be negated")
        elif isinstance(parent_formulae, Formula):
            return parent_formulae.bound_variables
        else:
            raise TypeError("pass a list of one formula")


NOT = NotFormulaBuilder()


class QuantifierFormulaBuilder(InductionMethod):
    def __init__(self, name, symbol):
        if isinstance(symbol, str):
            self.symbol = symbol
        else:
            raise TypeError("symbol must be a string")
        if isinstance(name, str):
            self.method = name
        else:
            raise TypeError("name must be a string")

    def string(self, parent_formulae):
        if isinstance(parent_formulae, list):
            if len(parent_formulae) == 2:
                if (
                        isinstance(parent_formulae[0], Formula)
                        and isinstance(parent_formulae[1], Terms.Term)
                        and parent_formulae[1].is_variable
                        ):
                    variable = parent_formulae[1]
                    formula = parent_formulae[0]
                elif (
                        isinstance(parent_formulae[0], Terms.Term)
                        and parent_formulae[0].is_variable
                        and isinstance(parent_formulae[1], Formula)
                        ):
                    variable = parent_formulae[0]
                    formula = parent_formulae[1]
                else:
                    raise TypeError("""
                                    pass a formula and a variable to be
                                    quantified"""
                                    )
                if variable.string in formula.bound_variables:
                    raise VariableClashError([variable.sring])
                else:
                    return (
                            '('
                            + self.symbol
                            + ' '
                            + variable.string
                            + ' '
                            + formula.string
                            + ')'
                            )
            else:
                raise ValueError("""
                                 list length is wrong! pass a list of a formula
                                 and a variable"""
                                 )
        else:
            raise ValueError("pass a list of a formula and a variable")

    def free_variables(self, parent_formulae):
        if isinstance(parent_formulae, list):
            if len(parent_formulae) == 2:
                if (
                        isinstance(parent_formulae[0], Formula)
                        and isinstance(parent_formulae[1], Terms.Term)
                        and parent_formulae[1].is_variable
                        ):
                    variable = parent_formulae[1]
                    formula = parent_formulae[0]
                elif (
                        isinstance(parent_formulae[0], Terms.Term)
                        and parent_formulae[0].is_variable
                        and isinstance(parent_formulae[1], Formula)
                        ):
                    variable = parent_formulae[0]
                    formula = parent_formulae[1]
                else:
                    raise TypeError("""
                                    pass a formula and a variable to be
                                    quantified"""
                                    )
                variables_list = copy.deepcopy(formula.free_variables)
                variables_list.remove(variable.string)
                return variables_list
            else:
                raise ValueError("""
                                 list length is wrong! pass a list of a formula
                                 and a variable"""
                                 )
        else:
            raise ValueError("pass a list of a formula and a variable")

    def bound_variables(self, parent_formulae):
        if isinstance(parent_formulae, list):
            if len(parent_formulae) == 2:
                if (
                        isinstance(parent_formulae[0], Formula)
                        and isinstance(parent_formulae[1], Terms.Term)
                        and parent_formulae[1].is_variable
                        ):
                    variable = parent_formulae[1]
                    formula = parent_formulae[0]
                elif (
                        isinstance(parent_formulae[0], Terms.Term)
                        and parent_formulae[0].is_variable
                        and isinstance(parent_formulae[1], Formula)
                        ):
                    variable = parent_formulae[0]
                    formula = parent_formulae[1]
                else:
                    raise TypeError("""
                                    pass a formula and a variable to be
                                    quantified"""
                                    )
                variables_list = copy.deepcopy(formula.bound_variables)
                variables_list.append(variable.string)
                return variables_list
            else:
                raise ValueError("""
                                 list length is wrong! pass a list of a formula
                                 and a variable"""
                                 )
        else:
            raise ValueError("pass a list of a formula and a variable")


EXISTS = QuantifierFormulaBuilder('EXISTS', 'EXISTS')
FORALL = QuantifierFormulaBuilder('FORALL', 'FORALL')


class Formula:
    def __init__(self, parent_formulae, induction_method):
        if isinstance(induction_method, InductionMethod):
            try:
                self.string = induction_method.string(parent_formulae)
            except VariableClashError as variable_clash:
                print(
                      """there was a clash of bound variables. please rename
                      the following variables: %s"""
                      % variable_clash.variable_list
                      )
            else:
                if isinstance(parent_formulae, list):
                    self.parent_formulae = copy.deepcopy(parent_formulae)
                elif isinstance(parent_formulae, Formula):
                    self.parent_formulae = [copy.deepcopy(parent_formulae)]
                else:
                    raise TypeError("""
                                    parent_formulae must be a list of
                                    formulae"""
                                    )
                self.free_variables = induction_method.free_variables(
                        self.parent_formulae)
                self.bound_variables = induction_method.bound_variables(
                        self.parent_formulae)
                self.induction_method = induction_method
        else:
            raise TypeError("""
                            pass a valid induction method for building
                            formulae"""
                            )

    def substitute(self, old_variable, new_term):
        if (
                isinstance(old_variable, Terms.Term)
                and old_variable.is_variable
                and isinstance(new_term, Terms.Term)
                ):
            if isinstance(self.induction_method, AtomicFormulaBuilder):
                new_parent_formulae = []
                for term in self.parent_formulae:
                    substituted_term = copy.deepcopy(term)
                    substituted_term.substitute(old_variable, new_term)
                    new_parent_formulae.append(substituted_term)
                Formula.__init__(
                                 self,
                                 new_parent_formulae,
                                 self.induction_method
                                 )
            else:
                new_parent_formulae = []
                for formula in self.parent_formulae:
                    if (
                            isinstance(formula, Terms.Term)
                            and formula.is_variable
                            and formula.string == old_variable.string
                            and (not new_term.is_variable)
                            ):
                        raise VariableClashError([
                                old_variable.string,
                                new_term.string
                                ])
                    formula_copy = copy.deepcopy(formula)
                    formula_copy.substitute(old_variable, new_term)
                    new_parent_formulae.append(formula_copy)
                Formula.__init__(
                                 self,
                                 new_parent_formulae,
                                 self.induction_method
                                 )
        else:
            raise TypeError("""
                            pass a variable to be replaced, and a term to
                            replace it with"""
                            )


EMPTYMETHOD = InductionMethod()
EMPTYFORMULA = Formula([], EMPTYMETHOD)
