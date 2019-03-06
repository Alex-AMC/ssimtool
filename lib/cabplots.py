import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression


def drawcabplots(df_data, title=None, label='APIvExp', xlim=(0, -10), ylim=(0, -10)):
    """
    Method plots the CAB plots from the statistical table passed. Function also generates a balance cohesive/adhesive line
    using the linear regression model. This is plotted as black on the graph.

    :param (obj)    df_data     : Takes Dataframe containing the statistical descriptors which are used for plotting
    :param (str)    title       : Passes the name of the title
    :param (str)    label       : Passes the labels for the excipients
    :param (floats) xlim        : Passes a list of two floats to set the x axis
    :param (floats) ylim        : Passes a list of two floats to set the y axis
    :return: Graph and of CAB Plot
    """
    x1 = np.asarray(df_data['Adhesion']).reshape(-1, 1)
    y1 = np.asarray(df_data['Cohesion']).reshape(-1, 1)

    lm1 = LinearRegression(fit_intercept=False)
    lm1.fit(x1, y1)

    x_test = np.linspace(-100, 0)
    y_pred1 = lm1.predict(x_test[:, None])
    plt.figure(figsize=(8, 6))

    plt.scatter(df_data['Adhesion'], df_data['Cohesion'], s=100, label=label, color='b')
    print("Gradient : " + str(lm1.coef_))
    print("R^2 : " + str(lm1.score(x1, y1)))

    plt.xlabel(r'$Adhesion Energy (mJ/m^2)$', fontsize=12)
    plt.ylabel(r'$Cohesion Energy (mJ/m^2)$', fontsize=12)
    plt.plot(x_test, y_pred1, 'blue', linewidth=2, label=label + '-Fit')
    plt.xlim(xlim[0], xlim[1])
    plt.ylim(ylim[0], ylim[1])
    plt.grid(True)
    plt.tick_params(axis='x', labelsize=12)
    plt.tick_params(axis='y', labelsize=12)
    plt.plot(x_test, 1 * x_test, 'black', linewidth=2)
    plt.legend(loc='upper left', fontsize=12)
    plt.title(title)
    plt.show()
