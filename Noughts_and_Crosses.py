# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 16:34:57 2019

Creating a noughts and crosses program that learns.
Initially, it picks a move totally at random; at the end of each game it
adjusts the probabilities with which it picks moves based on the outcome; so
after each loss, it becomes less likely to repeat the moves it played in that
game, and after each win, it becomes more likely to repeat those moves.

@author: gulli
"""

import numpy

import pandas

from copy import deepcopy


class NoughtsAndCrossesGameState():
    def __init__(self, game_state):
        if isinstance(game_state, numpy.ndarray):
            self.state = game_state
        elif isinstance(game_state, list):
            self.state = numpy.array(game_state)
        else:
            raise ValueError("game_state must be a list or numpy array")
        self.move_number = 1
        for row in range(3):
            for column in range(3):
                if (self.state[row][column] == 'X'
                        or self.state[row][column] == 'O'):
                    self.move_number += 1
        if self.move_number % 2 == 0:
            self.next_player = 'Crosses'
        else:
            self.next_player = 'Noughts'

    def checkForEnd(self, display=True):  # test for a win or draw
        players = ['Crosses', 'Noughts']
        triples = [['X', 'X', 'X'], ['O', 'O', 'O']]  # winning combinations
        for player in range(2):
            for index in range(3):
                if ((self.state[index] == triples[player]).all()  # row
                        or (self.state[0:3, index]  # column
                            == triples[player]).all()):
                    if display:
                        print("Game over! %s is the winner!" % players[player])
                    self.next_player = 'none'  # game over, so no next player
                    return '%s' % players[player]
            if ((self.state[[0, 1, 2], [0, 1, 2]]  # leading diagonal
                 == triples[player]).all()
                    or (self.state[[0, 1, 2], [2, 1, 0]]  # other diagonal
                        == triples[player]).all()):
                if display:
                    print("Game over! %s is the winner!" % players[player])
                self.next_player = 'none'  # game over, so no next player
                return '%s' % players[player]
        draw = True  # assume it's a draw unless you find a legal move
        for row in range(3):
            for column in range(3):
                if self.state[row][column] == ' ':
                    # if there's an empty cell, it's not yet a draw
                    draw = False
        if draw:
            if display:
                print("It's a draw! Game over!")
            self.next_player = 'none'  # game over, so no next player
            return 'draw'
        return 'still going'

    def displayBoard(self):
        # print out the current grid
        for i in range(3):
            print("--------------")
            print("| %s | %s | %s |" % (self.state[i][0],
                                        self.state[i][1],
                                        self.state[i][2]))
        print("--------------")

    def cellNumberToIndex(self, cell_number):
        # convert human-friendly cell number to array index
        cells = [
                 [0, 0],
                 [0, 1],
                 [0, 2],
                 [1, 0],
                 [1, 1],
                 [1, 2],
                 [2, 0],
                 [2, 1],
                 [2, 2]]
        return cells[cell_number - 1]

    def indexToCell(self, index_pair):
        # convert array index to human-friendly cell number
        cells = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        return cells[index_pair[0]][index_pair[1]]

    def stringDescription(self):  # return a string describing the game state
        description = ''
        NUMBER_OF_COLUMNS = 3
        NUMBER_OF_ROWS = 3
        for row in range(NUMBER_OF_COLUMNS):
            for column in range(NUMBER_OF_ROWS):
                description = description + self.state[row][column]  # append
                # symbol in game state cell
                if column < NUMBER_OF_COLUMNS - 1:
                    description = description + ','  # separate columns with ,
            if row < NUMBER_OF_ROWS - 1:
                description = description + ';'  # separate rows with ;
        return description


class NoughtsAndCrosses(NoughtsAndCrossesGameState):
    def __init__(self, display=True):
        NoughtsAndCrossesGameState.__init__(self, [
                                                   [' ', ' ', ' '],
                                                   [' ', ' ', ' '],
                                                   [' ', ' ', ' ']
                                                   ])  # set up empty game grid
        if display:  # display the empty board and prompt next turn
            self.displayBoard()
            print("%s' turn\n\n" % self.next_player)

    def playerMove(self, cell_number, display=True):
        # cell number is an int from 1 to 9, for human readability
        CELL_INDICES_LIST = [
                             [0, 0],
                             [0, 1],
                             [0, 2],
                             [1, 0],
                             [1, 1],
                             [1, 2],
                             [2, 0],
                             [2, 1],
                             [2, 2]]  # list of array indices for each cell
        # find the row/col indices for the given cell number
        cell_row_index = CELL_INDICES_LIST[cell_number - 1][0]
        cell_column_index = CELL_INDICES_LIST[cell_number - 1][1]

        # check the cell being played into is empty, and if it is, fill with
        # symbol based on current player, and advance the next player and
        # move number
        if self.state[cell_row_index, cell_column_index] == ' ':
            if self.next_player == 'Noughts':
                self.state[cell_row_index, cell_column_index] = 'O'
                self.next_player = 'Crosses'
                self.move_number += 1
            elif self.next_player == 'Crosses':
                self.state[cell_row_index, cell_column_index] = 'X'
                self.next_player = 'Noughts'
                self.move_number += 1
            else:
                print("Game over!")  # in case someone tries to keep playing
                # after the game is over
        else:  # if playing into an already filled cell
            print("Invalid move! Try again!")
            return False
        if display:
            self.displayBoard()
        self.checkForEnd(display)
        if display and (self.next_player != 'none'):
            print("%s' turn\n\n" % self.next_player)
        return True


class SingleStateProbabilityTable(NoughtsAndCrossesGameState):
    # store the possible moves and probabilities of playing them for a single
    # game state
    def __init__(self, game_state):
        if (isinstance(game_state, list)
                or isinstance(game_state, numpy.ndarray)):
            NoughtsAndCrossesGameState.__init__(self, game_state)
        elif isinstance(game_state, NoughtsAndCrossesGameState):
            self.state = game_state.state
            self.move_number = game_state.move_number
            self.next_player = game_state.next_player
        else:
            raise ValueError("""game_state must be a list, numpy array, or
                             NoughtsAndCrossesGameState""")

        moves = [[1], [2], [3], [4], [5], [6], [7], [8], [9]]
        number_legal_moves = 0

        # count how many empty cells are left
        for row in range(3):
            for column in range(3):
                if self.state[row, column] == ' ':
                    number_legal_moves += 1

        self.move_number = 9 - number_legal_moves + 1

        # initially, all legal moves are to be equally likely
        probability = 1/number_legal_moves
        for row in range(3):
            for column in range(3):
                if self.state[row, column] == ' ':
                    # if cell is legal, assign standard probability to it
                    moves[self.indexToCell([row, column])
                          - 1].append(probability)
                else:  # if cell already full, assign probability 0
                    moves[self.indexToCell([row, column]) - 1].append(0)
        self.probability_table = numpy.array(moves)

    def increaseProbability(self, cell_number, probability_change,
                            tolerance=0.0001):
        # increase the probability of playing a certain cell, and decrease all
        # other probabilities accordingly
        cell_index = cell_number - 1
        if self.probability_table[cell_index][1] == 1:
            return None  # if probability to be increased is already 1, finish

        # need to be sure not to increase probability above 1
        actual_change = min(probability_change,
                            1 - self.probability_table[cell_index][1])
        if actual_change == 0:
            return None

        # count up number of other moves whose probability is not yet 0, and
        # the minimum probability among such
        number_allowed_moves = 0
        min_probability = 1
        for entry in self.probability_table:
            if entry[1] != 0 and entry[0] != cell_number:
                number_allowed_moves += 1
                min_probability = min(min_probability, entry[1])

        if number_allowed_moves == 0:
            return None

        # all other non-zero probabilities must decrease by the same amount
        change_to_other_probabilities = (actual_change
                                         / number_allowed_moves)

        # if decreasing all other non-zero probabilities by this amount won't
        # make any of them negative, go ahead
        if change_to_other_probabilities <= min_probability:
            for entry in self.probability_table:
                if entry[0] == cell_number:
                    entry[1] += actual_change  # increase the given probability
                elif entry[1] != 0:
                    entry[1] -= change_to_other_probabilities
                    # subtract off the change
        else:  # another probability would become negative, so
            # the biggest amount we can safely increase the desired prob by
            actual_change = min_probability * number_allowed_moves
            for entry in self.probability_table:
                if entry[0] == cell_number:
                    entry[1] += actual_change  # increase as much as safely can
                elif entry[1] != 0:
                    entry[1] -= min_probability  # decrease others

            # now call the function again with the remaining increase that's
            # still needed, now that one of the other probs has been set to
            # 0 and won't cause a problem. But don't bother if the remaining
            # change is too small, otherwise we could get excessive calls with
            # ever smaller amounts.
            #if probability_change - actual_change >= tolerance:
                #self.increaseProbability(cell_number, (probability_change
                                                       #- actual_change))

    def decreaseProbability(self, cell_number, probability_change):
        # decrease the probability of playing a certain cell, and increase all
        # other probabilities accordingly
        cell_index = cell_number - 1
        # need to be sure not to decrease probability below 0
        actual_change = min(probability_change,
                            self.probability_table[cell_index][1])

        if actual_change == 0:
            return None  # if probability to be decreased is already 0, finish

        number_allowed_moves = 0
        for entry in self.probability_table:
            if entry[1] != 0 and entry[0] != cell_number:
                number_allowed_moves += 1

        if number_allowed_moves == 0:
            return None

        # all other non-zero probabilities must increase by the same amount
        change_to_other_probabilities = (actual_change
                                         / number_allowed_moves)
        # not possible to push another probability above 1, so don't need to
        # worry about that.
        for entry in self.probability_table:
            if entry[0] == cell_number:
                entry[1] -= actual_change  # decrease the given probability
            elif entry[1] != 0:
                entry[1] += change_to_other_probabilities  # add on the change


class SingleMoveProbabilityFrame():
    def __init__(self, move_number, previous_frame=0):
        # self.move_number = move_number

        if move_number == 1:  # first move
            game_state = NoughtsAndCrossesGameState([
                                                      [' ', ' ', ' '],
                                                      [' ', ' ', ' '],
                                                      [' ', ' ', ' ']
                                                      ])  # opening game state
            state_descriptions = [game_state.stringDescription()]
            self.probability_frame = pandas.DataFrame(
                    [SingleStateProbabilityTable(game_state)],
                    index=state_descriptions)  # all moves equiprobable
        else:
            # figure out who played last
            if move_number % 2 == 0:
                last_player = 'O'
            else:
                last_player = 'X'
            possible_states = []
            state_descriptions = []
            # look up previous possible states
            if previous_frame != 0:
                previous_move = deepcopy(previous_frame)
            else:
                previous_move = SingleMoveProbabilityFrame(move_number - 1)

            considered_states = []
            # loop through probability tables for previous move
            for table in previous_move.probability_frame[0]:
                if table.checkForEnd(False) == 'still going':
                    previous_state = table.state
                    for row in range(3):
                        for column in range(3):
                            # if a cell is still unfilled
                            if previous_state[row][column] == ' ':
                                new_state = deepcopy(previous_state)
                                # form a new_state with move from last player
                                # in this cell
                                new_state[row][column] = last_player
                                # take probability table for new_state
                                new_table = SingleStateProbabilityTable(
                                        new_state)
                                if (new_table.stringDescription() not in
                                        considered_states):
                                    considered_states.append(deepcopy(
                                            new_table.stringDescription()))
                                    possible_states.append(new_table)
                                    # take string description of new state
                                    state_descriptions.append(
                                            new_table.stringDescription())
            # gather together all possible states and their probabilities for
            # this move
            self.probability_frame = pandas.DataFrame(possible_states,
                                                      index=state_descriptions)


class FullProbabilityFrame(SingleMoveProbabilityFrame):
    def __init__(self):
        previous_frame = 0
        self.probability_frame = pandas.DataFrame([])
        for move in range(1, 10):
            # compute probability frame for this move
            previous_frame = SingleMoveProbabilityFrame(move, previous_frame)
            self.probability_frame = self.probability_frame.append(
                    previous_frame.probability_frame)


class ComputerPlayer(NoughtsAndCrosses):
    def __init__(self, game, player='Noughts'):
        self.game = game
        self.state = game.state
        self.move_number = game.move_number
        self.next_player = game.next_player
        self.player = player
        self.probability_frame = FullProbabilityFrame().probability_frame
        self.moves_record = []

    def makeMove(self, cell_number):
        if self.game.next_player == self.player:
            return self.game.playerMove(cell_number, False)
        else:
            print("It's not my turn!")

    def newGame(self, game, player='Noughts'):
        self.game = game
        self.state = game.state
        self.player = player
        self.moves_record = []

    def chooseMove(self):
        description = self.stringDescription()
        probability_table = self.probability_frame.loc[description][
                0].probability_table
        cell_numbers = probability_table[:, 0]
        probabilities = probability_table[:, 1]
        choice = numpy.random.choice(cell_numbers, p=probabilities)
        return int(choice)

    def updateRecord(self, string_description, cell_number):
        self.moves_record.append([string_description, cell_number])

    def updateGameState(self, game_state):
        self.state = game_state

    def resetProbabilityFrame(self):
        self.probability_frame = FullProbabilityFrame().probability_frame


class ComputerMatch(NoughtsAndCrosses):
    def __init__(self, noughts, crosses, reward=0.001, penalty=0.001,
                 display=False):
        NoughtsAndCrosses.__init__(self, False)
        self.noughts_player = noughts
        self.crosses_player = crosses
        self.noughts_player.newGame(self, 'Noughts')
        self.crosses_player.newGame(self, 'Crosses')
        while self.checkForEnd(False) == 'still going':
            if self.next_player == 'Noughts':
                next_player = self.noughts_player
            else:
                next_player = self.crosses_player
            next_player.updateGameState(self.state)
            move = next_player.chooseMove()
            next_player.updateRecord(self.stringDescription(), move)
            next_player.makeMove(move)

        if self.checkForEnd(display) == 'Noughts':
            winner = self.noughts_player
            loser = self.crosses_player
            self.outcome = 'Noughts'
        elif self.checkForEnd(False) == 'Crosses':
            winner = self.crosses_player
            loser = self.noughts_player
            self.outcome = 'Crosses'
        else:
            self.outcome = 'draw'
            return None
        for record in winner.moves_record:
            winner.probability_frame.loc[record[0]][0].increaseProbability(
                    record[1], reward)
        for record in loser.moves_record:
            loser.probability_frame.loc[record[0]][0].decreaseProbability(
                    record[1], penalty)
