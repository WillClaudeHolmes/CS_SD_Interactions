import math

def format_plot(a,ts):
    """Format a plot for presenation.  Function takes one parameter, a matplotlib.axes object."""

    goldenratio = 1 / 2 * (1 + math.sqrt(5))  # The next few lines are used for the size of plots
    fsx = 7  # Width (in inches) for the figures.
    fsy = fsx / goldenratio  # Height (in inches) for the figures.

    a.tick_params(labelsize=ts)
    a.xaxis.label.set_size(ts+2)
    a.yaxis.label.set_size(ts+2)
    a.title.set_size(ts+4)
    a.legend(fontsize=ts)
    a.get_figure().set_size_inches(fsx, fsy)
    a.grid(1)


def two_pop_model(dt, timemax, r, a, pp, T_Interval = 1.0):
    """Model for two species population growth with interaction.
    Parameters:
        dt - Timestep
        timemax - time to which the simulation will run
        r - 2 x 2 list of coeficients which are multiplied by the population
            terms to determine growth rate.
        a - 2 x 2 list of coeficients which are multiplied by the population squared
            terms to determine growth rate.
        pp - 2 x 1 list of initial populatino values.
        T_Interval - (Optional Keyword argument).  Time between data being returned.
          if not given, defaults to 1.
    Returns:
        Touple of time and population values: (t[n], [p[0],p[1]][n])
    """

    #set initial values
    p = [pp[0],pp[1]]
    t = 0
    t_array = [0]
    population_array = [[p[0],p[1]]]

    while t < timemax:
        #Store current population values so all calculations are based on values from
        #  previous time step.
        ptn1=[p[0],p[1]]

        #Calculate new population values.
        p[0] += (r[0][0]*ptn1[0] + r[0][1]*ptn1[1] - a[0][0]*ptn1[0]**2 - a[0][1]*ptn1[0]*ptn1[1])*dt
        p[1] += (r[1][1]*ptn1[1] + r[1][0]*ptn1[0] - a[1][1]*ptn1[1]**2 - a[1][0]*ptn1[1]*ptn1[0])*dt
        #Increment time step.
        t += dt
        #Only store values ever Dt = 1.
        if abs(t % T_Interval - T_Interval) <= dt:
            t_array.append(t)
            population_array.append([p[0],p[1]])

    return t_array, population_array

def pop_model(dt, timemax, r, a, pp, T_Interval = 1):
        """Model for an N species population model.
        The parameters are the same as two_pop_model"""

        N=len(pp)
        p = [pp[j] for j in range(N)]
        t = 0
        t_array = [0]
        population_array = [p[:]]

        while t < timemax:
            ptn1 = p[:]
            for j in range(N):
                delta = sum([r[j][k] * ptn1[k]           for k in range(N)])
                delta -= sum([a[j][k] * ptn1[j] * ptn1[k] for k in range(N)])
                p[j] += delta * dt
            t += dt
            if abs(t % T_Interval) < dt:
                t_array.append(t)
                population_array.append(p[:])

        return t_array, population_array