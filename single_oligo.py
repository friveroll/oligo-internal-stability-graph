import re
import matplotlib.pyplot as plt


def get_delta_g_pentamer(pentamer):
    """
       Return the deltaG of the oligo's pentamers.

       param: oligo: a string of nucleotides
       return: a list of floats

    """
    if not isinstance(pentamer, str):
        raise TypeError("pentamer must be a string with 5 nucleotides")
    if not bool(re.compile('^[ACTG]+$').match(pentamer)):
        raise ValueError("pentamer must contain only A, C, T, G")
    if len(pentamer) != 5:
        raise ValueError("Pentamer must be 5 nucleotides long")

    pentamer = pentamer.upper()
    delta_g = 0
    near_pair_delta_g = {'AA': -1.9, 'AC': -1.3, 'AG': -1.6, 'AT': -1.5,
                         'CA': -1.9, 'CC': -3.1, 'CG': -3.6, 'CT': -1.6,
                         'GA': -1.6, 'GC': -3.1, 'GG': -3.1, 'GT': -1.3,
                         'TA': -1.0, 'TC': -1.6, 'TG': -1.9, 'TT': -1.9}

    pentamer_two_mers = get_two_mers(pentamer)
    for two_mer in pentamer_two_mers:
        delta_g += pentamer_two_mers.get(two_mer) * \
            near_pair_delta_g.get(two_mer)
    return delta_g


def get_pentamers(oligo):
    """
           Return the oligo's pentamers.

           param: oligo: a string of nucleotides
           return: a list of strings for each pentamer
           in the oligo
        """
    return [oligo[i:i+5] for i in range(len(oligo)-4)]


def get_two_mers(oligo):
    """
           Return the oligo's two mers.

           param: oligo: a string of nucleotides
           return: a dict of strings with counts
           for each two mer in the oligo
        """
    two_mers = [oligo[i:i+2]
                for i in range(len(oligo)) if len(oligo[i:i+2]) >= 2]
    return {x: two_mers.count(x) for x in two_mers}


def get_delta_g_pentamers(oligo):
    """
       Return the deltaG of the oligo.

       param: oligo: a string of nucleotides
       return: list of floats with deltaG for each pentamer
    """
    if not isinstance(oligo, str):
        raise TypeError("oligo must be a string of nucleotides")
    if not bool(re.compile('^[ACTG]+$').match(oligo)):
        raise ValueError("oligo must contain only A, C, T, G")

    oligo = oligo.upper()
    pentamers = get_pentamers(oligo)
    return [get_delta_g_pentamer(pentamer) for pentamer in pentamers]


def get_stability_graph(oligo):
    """
       Return a graph of the oligo's pentamers and their deltaG.

       param: oligo: a string of nucleotides
       return: a matplotlib graph
    """
    if not isinstance(oligo, str):
        raise TypeError("oligo must be a string of nucleotides")
    if not bool(re.compile('^[ACTG]+$').match(oligo)):
        raise ValueError("oligo must contain only A, C, T, G")

    oligo = oligo.upper()
    pentamers = get_pentamers(oligo)
    delta_g_pentamers = get_delta_g_pentamers(oligo)
    x_labels = list(oligo[0:-5])
    x_labels.append(pentamers[-1])

    plt.plot(delta_g_pentamers, marker='o', linestyle=':', color='b')

    plt.ylabel(r'$\sum \Delta$ G $\left(\frac{kcal}{mol}\right)$ pentamers')
    plt.xticks(range(len(x_labels)), x_labels)
    plt.gca().axes.set_ylim([-12, -5])
    plt.gca().invert_yaxis()
    plt.savefig('stability_graph.png', bbox_inches='tight')
    # plt.show()


if __name__ == '__main__':
    get_stability_graph("ACTTGGGATTGGGCT")
