import scipy.optimize as spo




# This is the function that will be tested by the autograder
# The student must update this code to properly implement the functionality
def optimize_portfolio(sd=dt.datetime(2008,1,1), ed=dt.datetime(2009,1,1), \
    syms=['GOOG','AAPL','GLD','XOM'], gen_plot=False):

    # Read in adjusted closing prices for given symbols, date range
    dates = pd.date_range(sd, ed)
    prices_all = get_data(syms, dates)  # automatically adds SPY
    prices = prices_all[syms]  # only portfolio symbols
    prices_SPY = prices_all['SPY']  # only SPY, for comparison later

    if syms == []:
        return None

    # find the allocations for the optimal portfolio
    # note that the values here ARE NOT meant to be correct for a test case
    allocs = np.asarray([1.0/len(syms)]*len(syms)) # add code here to find the allocations
    normlized_prices = normalize(prices)

    constraints = ({'type': 'eq', 'fun': lambda x: 1.0 - np.sum(x)})
    bounds = tuple([(0.0, 1.0) for i in range(len(syms))])

    allocs = spo.minimize(f, allocs, args=(normlized_prices,), method='SLSQP', bounds=bounds, constraints=constraints).x


def distance(data, clusters, allocs):
    """
    :param data: ndarray contains
    :param clusters:
    :param allocs:
    :return:
    """
    predict_signal = np.zeros(clusters.shape[1])
    return (-1.0)*(average - rfr)/std * math.sqrt(sf)