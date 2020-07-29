"""
To be run with the Noughts_and_Crosses module; this sets up two players, trains
them, then tests them out and checks their win ratios.
With perfect play from both players, noughts and crosses always ends in a draw;
so after training, the final outcome should always be a draw.
"""

from Noughts_and_Crosses import *

import time

import numpy

start_time = time.time()

game = NoughtsAndCrosses(display=False)
player_one = ComputerPlayer(game)
player_two = ComputerPlayer(game)
player_one_wins = 0
player_two_wins = 0
draws = 0
outcomes = []

training_games = 1000
testing_games = 1000


# training with player_one going first
for i in range(training_games):
    match = ComputerMatch(player_one, player_two, reward=0.1, penalty=0.1)
    if match.outcome == 'Noughts':
        player_one_wins += 1
        outcomes.append(1)
    elif match.outcome == 'Crosses':
        player_two_wins += 1
        outcomes.append(2)
    else:
        draws += 1
        outcomes.append(0)



# training with player_two going first
for i in range(training_games):
    match = ComputerMatch(player_two, player_one, reward=0.1, penalty=0.1)
    if match.outcome == 'Noughts':
        player_two_wins += 1
        outcomes.append(2)
    elif match.outcome == 'Crosses':
        player_one_wins += 1
        outcomes.append(1)
    else:
        draws += 1
        outcomes.append(0)


print("Training complete with %d games played\n" % (2*training_games))
print("Player one won %d games in training" % player_one_wins)
print("Player two won %d games in training" % player_two_wins)
print("There were %d draws in training\n\n" % draws)

player_one_wins = 0
player_two_wins = 0
draws = 0
outcomes = []


# testing with player_one going first
for i in range(testing_games):
    match = ComputerMatch(player_one, player_two)
    if match.outcome == 'Noughts':
        player_one_wins += 1
        outcomes.append(1)
    elif match.outcome == 'Crosses':
        player_two_wins += 1
        outcomes.append(2)
    else:
        draws += 1
        outcomes.append(0)



# testing with player_two going first
for i in range(testing_games):
    match = ComputerMatch(player_two, player_one, reward=0.01, penalty=0.01)
    if match.outcome == 'Noughts':
        player_two_wins += 1
        outcomes.append(2)
    elif match.outcome == 'Crosses':
        player_one_wins += 1
        outcomes.append(1)
    else:
        draws += 1
        outcomes.append(0)


print("Testing complete with %d games played\n" % (2*testing_games))
print("Player one won %d games in testing" % player_one_wins)
print("Player two won %d games in testing" % player_two_wins)
print("There were %d draws in testing" % draws)



print("--- %s seconds ---" % (time.time() - start_time))
