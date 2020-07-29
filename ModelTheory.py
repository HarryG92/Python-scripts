# -*- coding: utf-8 -*-
"""
Created on Wed Mar 25 14:44:36 2020

Part of a Model Theory program. Tools for building up terms and formulae are
provided in prior files; this file is concerned with coding the rules of
inference for proofs. The idea is to have a virtual class, Deductions, then
subclasses for each rule of inference (e.g., modus ponens); then an individual
application of the rule is an object of the class. So if A and B are formulae,
we represent the deduction "A; A->B: therefore B" as ModusPonens(A,A->B), which
has B as a "conclusion" attribute.

@author: Harry Gulliver
"""


import copy

import ModelTheoryTerms as Terms

import ModelTheoryFormulae as Formulae


class Deduction:
    def __init__(self):
        self.input_formulae = []
        self.conclusion = Formulae.EMPTYFORMULA


""" propositional calculus rules of inference """


# given a->b and a, deduce b
class ModusPonens(Deduction):
    def __init__(self, a_implies_b, a):
        if (
                isinstance(a_implies_b, Formulae.Formula)
                and isinstance(a, Formulae.Formula)
                ):
            self.input_formulae = [copy.deepcopy(a_implies_b),
                                   copy.deepcopy(a)
                                   ]
            if self.input_formulae[0].induction_method == Formulae.IMPLIES:
                if (
                        self.input_formulae[0].parent_formulae[0].string
                        == self.input_formulae[1].string
                        ):
                    self.conclusion = self.input_formulae[0].parent_formulae[1]
                else:
                    raise ValueError("""
                                     second formula must be antecedent of
                                     first: the A of A->B"""
                                     )
            else:
                raise ValueError("""
                                 first formula must be of method IMPLIES; so
                                 of the form A-> B"""
                                 )
        else:
            raise TypeError("pass two formulae")


# given a->b and ¬b, deduce a
class ModusTollens(Deduction):
    def __init__(self, a_implies_b, not_b):
        if (
                isinstance(a_implies_b, Formulae.Formula)
                and isinstance(not_b, Formulae.Formula)
                ):
            self.input_formulae = [copy.deepcopy(a_implies_b),
                                   copy.deepcopy(not_b)
                                   ]
            if self.input_formulae[0].induction_method == Formulae.IMPLIES:
                implication_negated_consequent = Formulae.Formula(
                        [self.input_formulae[0].parent_formulae[1]],
                        Formulae.NOT
                        )
                if (
                        implication_negated_consequent.string
                        == self.input_formulae[1].string
                        ):
                    self.conclusion = Formulae.Formula(
                            [self.input_formulae[0].parent_formulae[0]],
                            Formulae.NOT
                            )
                else:
                    raise ValueError("""
                                     second formula must be negation of the
                                     consequent of the first: the ¬B of A->B"""
                                     )
            else:
                raise ValueError("""
                                 first formula must be of method IMPLIES; so
                                 of the form A-> B"""
                                 )
        else:
            raise TypeError("pass two formulae")


# given ¬¬a, deduce a
class DoubleNegative(Deduction):
    def __init__(self, double_negative):
        if isinstance(double_negative, Formulae.Formula):
            self.input_formulae = [copy.deepcopy(double_negative)]
            if self.input_formulae[0].induction_method == Formulae.NOT:
                single_negative = self.input_formulae[0].parent_formulae[0]
                if single_negative.induction_method == Formulae.NOT:
                    self.conclusion = single_negative.parent_formulae[0]
                else:
                    raise ValueError("input formula must be a DOUBLE negative")
            else:
                raise ValueError("input formula must be a double negative")
        else:
            raise TypeError("pass a formula")


# given a and b, deduce (a AND b)
class ConjunctionJoin(Deduction):
    def __init__(self, first_conjunct, second_conjunct):
        if (
                isinstance(first_conjunct, Formulae.Formula)
                and isinstance(second_conjunct, Formulae.Formula)
                ):
            self.input_formulae = [copy.deepcopy(first_conjunct),
                                   copy.deepcopy(second_conjunct)
                                   ]
            self.conclusion = Formulae.Formula(
                                               self.input_formulae,
                                               Formulae.AND
                                               )
        else:
            raise TypeError("Pass two formulae to be conjoined")


# given (a AND b), deduce a
class FirstConjunct(Deduction):
    def __init__(self, conjunction):
        if isinstance(conjunction, Formulae.Formula):
            self.input_formulae = [copy.deepcopy(conjunction)]
            if self.input_formulae[0].induction_method == Formulae.AND:
                self.conclusion = self.input_formulae[0].parent_formulae[0]
            else:
                raise ValueError("must pass a conjunction")
        else:
            raise TypeError("pass a conjunction to extract its first conjunct")


# given (a AND b), deduce b
class SecondConjunct(Deduction):
    def __init__(self, conjunction):
        if isinstance(conjunction, Formulae.Formula):
            self.input_formulae = [copy.deepcopy(conjunction)]
            if self.input_formulae[0].induction_method == Formulae.AND:
                self.conclusion = self.input_formulae[0].parent_formulae[1]
            else:
                raise ValueError("must pass a conjunction")
        else:
            raise TypeError("pass a conjunction to extract second conjunct")


# given a and b, deduce (a OR b)
class DisjunctionJoin(Deduction):
    def __init__(self, first_disjunct, second_disjunct):
        if (
                isinstance(first_disjunct, Formulae.Formula)
                and isinstance(second_disjunct, Formulae.Formula)
                ):
            self.input_formulae = [copy.deepcopy(first_disjunct),
                                   copy.deepcopy(second_disjunct)
                                   ]
            self.conclusion = Formulae.Formula(
                                               self.input_formulae,
                                               Formulae.OR
                                               )
        else:
            raise TypeError("Pass two formulae to be disjoined")


# given a->z, b->z, and (a OR b), deduce z
class CaseAnalysis(Deduction):
    def __init__(self, a_implies_z, b_implies_z, a_or_b):
        if (
                isinstance(a_implies_z, Formulae.Formula)
                and isinstance(b_implies_z, Formulae.Formula)
                and isinstance(a_or_b, Formulae.Formula)
                ):
            if (
                    a_implies_z.induction_method == Formulae.IMPLIES
                    and b_implies_z.induction_method == Formulae.IMPLIES
                    ):
                if (
                        a_implies_z.parent_formulae[1].string
                        == b_implies_z.parent_formulae[1].string
                        ):
                    if isinstance(a_or_b, Formulae.Formula):
                        a_or_b_copy = copy.deepcopy(a_or_b)
                        a_implies_z_copy = copy.deepcopy(a_implies_z)
                        b_implies_z_copy = copy.deepcopy(b_implies_z)
                        first_disjunct = a_or_b_copy.parent_formulae[0].string
                        second_disjunct = a_or_b_copy.parent_formulae[1].string
                        first_antecedent = (
                                a_implies_z_copy.parent_formulae[0].string
                                )
                        second_antecedent = (
                                b_implies_z_copy.parent_formulae[0].string
                                )
                        if (
                                (
                                 first_disjunct == first_antecedent
                                 and second_disjunct == second_antecedent
                                 )
                                or (
                                    first_disjunct == second_antecedent
                                    and second_disjunct == first_antecedent
                                    )
                                ):
                            self.input_formulae = [
                                                   a_implies_z_copy,
                                                   b_implies_z_copy,
                                                   a_or_b_copy
                                                   ]
                            self.conclusion = (
                                    a_implies_z_copy.parent_formulae[1]
                                    )
                        else:
                            raise ValueError("""
                                             the disjuncts of the third formula
                                             must be the antecendents of the
                                             implications"""
                                             )
                    else:
                        raise TypeError("final formula must be disjunction")
                else:
                    raise ValueError("""
                                     both implications must have the same
                                     consequent"""
                                     )
            else:
                raise ValueError("first two formulae must be implications")
        else:
            raise TypeError("pass three formulae")


# given (a OR b) and ¬a, deduce b (or, given (a OR b) and ¬b, deduce a)
class DisjunctiveSyllogism(Deduction):
    def __init__(self, a_or_b, not_a):
        if (
                isinstance(a_or_b, Formulae.Formula)
                and isinstance(not_a, Formulae.Formula)
                ):
            if not_a.induction_method == Formulae.NOT:
                not_a_copy = copy.deepcopy(not_a)
                a_or_b_copy = copy.deepcopy(a_or_b)
                negated_formula = not_a_copy.parent_formulae[0].string
                first_disjunct = a_or_b_copy.parent_formulae[0].string
                second_disjunct = a_or_b_copy.parent_formulae[1].string
                if negated_formula == first_disjunct:
                    self.input_formulae = [a_or_b_copy, not_a_copy]
                    self.conclusion = a_or_b_copy.parent_formulae[1]
                elif negated_formula == second_disjunct:
                    self.input_formulae = [a_or_b_copy, not_a_copy]
                    self.conclusion = a_or_b_copy.parent_formulae[0]
                else:
                    raise ValueError("""
                                     second formula must be the negation of a
                                     disjunct of the first formula"""
                                     )
            else:
                raise ValueError("second formula must be a negative")
        else:
            raise TypeError("Pass two formulae")


# given a->x, b->y, and (a OR b), deduce (x OR y)
class ConstructiveDilemma(Deduction):
    def __init__(self, a_implies_x, b_implies_y, a_or_b):
        if (
                isinstance(a_implies_x, Formulae.Formula)
                and isinstance(b_implies_y, Formulae.Formula)
                and isinstance(a_or_b, Formulae.Formula)
                ):
            if (
                    a_implies_x.induction_method == Formulae.IMPLIES
                    and b_implies_y.induction_method == Formulae.IMPLIES
                    and a_or_b.induction_method == Formulae.OR
                    ):
                a_implies_x_copy = copy.deepcopy(a_implies_x)
                b_implies_y_copy = copy.deepcopy(b_implies_y)
                a_or_b_copy = copy.deepcopy(a_or_b)
                first_antecedent = a_implies_x_copy.parent_formulae[0]
                second_antecedent = b_implies_y_copy.parent_formulae[0]
                first_disjunct = a_or_b_copy.parent_formulae[0]
                second_disjunct = a_or_b_copy.parent_formulae[1]
                first_consequent = a_implies_x_copy.parent_formulae[1]
                second_consequent = b_implies_y_copy.parent_formulae[1]
                if (
                        first_antecedent.string == first_disjunct.string
                        and second_antecedent.string == second_disjunct.string
                        ):
                    self.input_formulae = [a_implies_x_copy,
                                           b_implies_y_copy,
                                           a_or_b_copy
                                           ]
                    self.conclusion = Formulae.Formula(
                                                       [
                                                        first_consequent,
                                                        second_consequent
                                                        ],
                                                       Formulae.OR
                                                       )
                elif (
                        first_antecedent.string == second_disjunct.string
                        and second_antecedent.string == first_disjunct.string
                        ):
                    self.input_formulae = [b_implies_y_copy,
                                           a_implies_x_copy,
                                           a_or_b_copy
                                           ]
                    self.conclusion = Formulae.Formula(
                                                       [
                                                        second_consequent,
                                                        first_consequent
                                                        ],
                                                       Formulae.OR
                                                       )
                else:
                    raise ValueError("""
                                     the disjuncts in the third formula must be
                                     the antecedents of the implications"""
                                     )
            else:
                raise ValueError("""first two formulae must be implications,
                                 third must be disjunction"""
                                 )
        else:
            raise TypeError("pass three formulae")


# given a<->b, deduce (a->b AND b->a)
class IffSplitter(Deduction):
    def __init__(self, a_iff_b):
        if isinstance(a_iff_b, Formulae.Formula):
            if a_iff_b.induction_method == Formulae.IFF:
                a_iff_b_copy = copy.deepcopy(a_iff_b)
                forward_implication = Formulae.Formula(
                        a_iff_b_copy.parent_formulae,
                        Formulae.IMPLIES
                        )
                backward_implication = Formulae.Formula(
                        a_iff_b_copy.parent_formulae[::-1],
                        Formulae.IMPLIES
                        )
                self.input_formulae = [a_iff_b_copy]
                self.conclusion = Formulae.Formula(
                                                   [
                                                    forward_implication,
                                                    backward_implication
                                                    ],
                                                   Formulae.AND)
            else:
                raise ValueError("pass an iff formula")
        else:
            raise TypeError("pass a formula")


# given a->b and b->a, deduce a<->b
class IffJoin(Deduction):
    def __init__(self, a_implies_b, b_implies_a):
        if (
                isinstance(a_implies_b, Formulae.Formula)
                and isinstance(b_implies_a, Formulae.Formula)
                ):
            a_implies_b_copy = copy.deepcopy(a_implies_b)
            b_implies_a_copy = copy.deepcopy(b_implies_a)
            first_antecedent = a_implies_b_copy.parent_formulae[0]
            first_consequent = a_implies_b_copy.parent_formulae[1]
            second_antecedent = b_implies_a_copy.parent_formulae[0]
            second_consequent = b_implies_a_copy.parent_formulae[1]
            if (
                    first_consequent.string == second_antecedent.string
                    and first_antecedent.string == second_consequent.string
                    ):
                self.input_formulae = [a_implies_b_copy, b_implies_a_copy]
                self.conclusion = Formulae.Formula(
                                                   [
                                                    first_antecedent,
                                                    first_consequent
                                                    ],
                                                   Formulae.IFF
                                                   )
            else:
                raise ValueError("pass an implication and its converse")
        else:
            raise TypeError("Pass two formulae to be conjoined")


""" predicate calculus rules of inference """


# given [FORALL x phi(x)] and a term b, deduce phi(b)
class UniversalInstantiation(Deduction):
    def __init__(self, forall_x_phi, term):
        if not (
                isinstance(forall_x_phi, Formulae.Formula)
                and isinstance(term, Terms.Term)
                ):
            raise TypeError("pass a formula and a term")
        if forall_x_phi.induction_method != Formulae.FORALL:
            raise ValueError("first argument must be a universal formula")
        self.input_formulae = [copy.deepcopy(forall_x_phi),
                               copy.deepcopy(term)
                               ]
        quantified_variable = self.input_formulae[0].parent_formulae[0]
        self.conclusion = copy.deepcopy(
                self.input_formulae[0].parent_formulae[1]
                )
        self.conclusion.substitute(quantified_variable, self.input_formulae[1])


# given phi(x) for a free variable x, deduce [FORALL x phi(x)]
class UniversalGeneralisation(Deduction):
    def __init__(self, free_variable, general_formula):
        if not (
                isinstance(general_formula, Formulae.Formula)
                and isinstance(free_variable, Terms.Term)
                and free_variable.is_variable
                ):
            raise TypeError("""
                            pass a free variable and a formula to be
                            generalised"""
                            )
        if free_variable.string not in general_formula.free_variables:
            raise ValueError("the variable must be free in the formula")
        self.input_formulae = [
                               copy.deepcopy(free_variable),
                               copy.deepcopy(general_formula)
                               ]
        self.conclusion = Formulae.Formula(self.input_formulae,
                                           Formulae.FORALL
                                           )


# given phi(b) for some term b, deduce [EXISTS x phi(x)]
class ExistentialGeneralisation(Deduction):
    def __init__(self, variable, term, formula_of_variable, formula_of_term):
        if not (
                isinstance(variable, Terms.Term)
                and variable.is_variable
                and isinstance(term, Terms.Term)
                and isinstance(formula_of_variable, Formulae.Formula)
                and isinstance(formula_of_term, Formulae.Formula)
                ):
            raise TypeError("""pass a variable, a term, a formula in the
                            variable, and that formula with the term instead"""
                            )
        formula_of_variable_copy = copy.deepcopy(formula_of_variable)
        formula_of_variable_copy.substitute(variable, term)
        if formula_of_variable_copy.string != formula_of_term.string:
            raise ValueError("""
                             the second formula must be the first, but with the
                             term substituted in place of the string"""
                             )
        self.input_formulae = [copy.deepcopy(variable),
                               copy.deepcopy(term),
                               copy.deepcopy(formula_of_variable),
                               copy.deepcopy(formula_of_term)]
        self.conclusion = Formulae.Formula(
                                           [
                                            self.input_formulae[0],
                                            self.input_formulae[2]
                                            ],
                                           Formulae.EXISTS
                                           )
