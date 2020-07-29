# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 11:00:04 2020

A handful of basic machine learning models that I coded myself for practice;
nothing here that's not available in sklearn. At some point I might go back and
do more on this, but for now it's just a few starter pieces

@author: gulli
"""

import numpy as np

import pandas as pd


def input_formatter(data):
    if isinstance(data, pd.DataFrame):
        return data
    elif isinstance(data, pd.Series):
        return pd.DataFrame(data)
    elif isinstance(data, np.ndarray):
        return pd.DataFrame(data)
    else:
        raise TypeError("Wrong argument type")


class RegressionModel:
    def fit(self, X_train, y_train):
        return None

    def predict(self, x_test):
        return None


class LeastSquaresLinearRegression(RegressionModel):
    def __init__(self):
        self.params = np.array([])

    def fit(self, X_train, y_train):
        x_data = input_formatter((X_train))
        sample_size = x_data.shape[0]
        self.num_features = x_data.shape[1]
        X = np.ones((sample_size, self.num_features + 1))
        for sample in range(sample_size):
            X[sample][1:] = x_data.iloc[sample]
        XT = np.transpose(X)
        self.params = np.dot(
                             np.dot(
                                    np.linalg.inv(np.dot(XT, X)),
                                    XT),
                             y_train)

    def predict(self, X_test):
        x_data = input_formatter(X_test)
        test_size = x_data.shape[0]
        predictions = []
        for test in range(test_size):
            y = self.params[0]
            for feature in range(self.num_features):
                y += self.params[feature + 1] * x_data.iloc[test][feature]
            predictions.append(y)
        return np.array(predictions)


def logistic(x):
    return 1 / (1 + np.exp(-x))



class LogisticRegression(RegressionModel):
    def __init__(self):
        self.params = np.array([])

    def fit(self, X_train, y_train):
        x_data = input_formatter(X_train)
        sample_size = x_data.shape[0]
        self.num_features = x_data.shape[1]
        def grad_logl(params):
            grad = np.zeros(len(params))
            grad_params = np.array(params[1:])
            for datum in range(sample_size):
                grad[0] += (1
                            / (1
                               + np.exp(params[0]
                                        + (grad_params * x_data[datum])
                                        )
                               )
                            - 1
                            + y_train[datum]
                            )
                for feature in range(self.num_features):
                    grad[feature] += (x_data.iloc[datum][feature]
                                      / (1
                                         + np.exp(params[0]
                                                  + (grad_params
                                                     * x_data[datum])
                                                  )
                                         )
                                      - (x_data.iloc[datum][feature]
                                         * (1 - y_train[datum])
                                         )
                                      )
            return grad

        def jacobian(params):
            jacobian = numpy.zeros((len(params), len(params)))
            for row in range(len(params)):
                if row == 0:
                    for datum in range(sample_size):
                        jacobian[row, row] += (np.exp(params[0]
                                                      + (grad_params
                                                         * x_data[datum])
                                                      )
                                               / (1
                                                  + np.exp(params[0]
                                                           + (grad_params
                                                              * x_data[datum]
                                                              )
                                                           )
                                                  ) ** 2
                else:
                    for datum in range(sample_size):
                        jacobian[row, row] += (x_data.iloc[datum][row] ** 2
                                               * np.exp(params[0]
                                                        + (grad_params
                                                           * x_data[datum])
                                                        )
                                               / (1
                                                  + np.exp(params[0]
                                                           + (grad_params
                                                              * x_data[datum]
                                                              )
                                                           )
                                                  ) ** 2
                for column in range(row + 1, len(params)):
                    jacobian[row, column] += (x_data.iloc[datum][row]
                                              * x_data.iloc[datum][column]
                                              * np.exp(params[0]
                                                       + (grad_params
                                                          * x_data[datum])
                                                       )
                                              / (1
                                                 + np.exp(params[0]
                                                          + (grad_params
                                                             * x_data[datum]
                                                             )
                                                          )
                                                 ) ** 2
            return jacobian

        def newton_raphson_iterator(previous_iterate):
            return previous_iterate + np.dot(
                                             np.linalg.inv(
                                                 jacobian(previous_iterate)),
                                             grad_logl(previous_iterate)
                                             )

        
