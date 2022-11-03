import re
import matplotlib.pyplot as plt
from decimal import Decimal
import decimal

decimal.getcontext().prec = 3


def get_oligo_kmers(sequence, k):
    """Get all possible substrings of length k contained within a biological sequence  

    :param str sequence a text string with a nucleotide sequence
    :param int k an integer greater than 1 representing the frequency k-mer 
    :return list a list of all possible k-mers

    Usage examples:
    >>> get_oligo_kmers("ACGAT", 2)
    ['AC', 'CG', 'GA', 'AT']
    >>> get_oligo_kmers("ACTTGGGATTGGGCT", 5)
    ['ACTTG', 'CTTGG', 'TTGGG', 'TGGGA', 'GGGAT', 'GGATT', 'GATTG', 'ATTGG', 'TTGGG', 'TGGGC', 'GGGCT']
    """
    if k < 2:
        raise ValueError("k must be greater than 1")
    if k > len(sequence):
        raise ValueError("k must be less than the length of the sequence")

    return [sequence[i:i+k] for i in range(len(sequence)) if len(sequence[i:i+k]) == k]


def get_delta_g(sequence):
    """Obtiene el valor de Delta G para una secuencia de nucleótidos

    :param str sequence a text string with a nucleotide sequence
    :return float Delta G value calculated for the nucleotide sequence

    Usage examples:
    >>> get_delta_g('ACTTG')
    -6.7
    """
    if not isinstance(sequence, str):
        raise TypeError("sequence  must be a string for an oligonucleotide")
    if not bool(re.compile('^[ACTG]+$').match(sequence)):
        raise ValueError("sequence must contain only A, C, T, G")

    near_pairs_delta_g = {'AA': -1.9, 'AC': -1.3, 'AG': -1.6, 'AT': -1.5,
                          'CA': -1.9, 'CC': -3.1, 'CG': -3.6, 'CT': -1.6,
                          'GA': -1.6, 'GC': -3.1, 'GG': -3.1, 'GT': -1.3,
                          'TA': -1.0, 'TC': -1.6, 'TG': -1.9, 'TT': -1.9}
    delta_g = sum([Decimal(near_pairs_delta_g.get(dimer))
                  for dimer in get_oligo_kmers(sequence, 2)])
    return float(delta_g)


def get_delta_g_pentamers(sequence):
    """Obtiene el valor de Delta G para cada uno de los pentámeros

    :param srt sequence a text string with A C T G letters for a 
    oligo sequence with more than 5 nucleotides
    :retrun list a list with Delta G value for each perntamer in the oligo

    Usage example:
    get_delta_g_pentamers("TCTTGT")
    [-7.0, -6.7]

    get_delta_g_pentamers("ACTTGGGATTGGGCT")
    [-6.700000000000001, -8.5, -10.0, -9.7, -9.3, -8.1, -6.9, -8.4, -10.0, -11.2, -10.9]

    get_delta_g_pentamers("TAATACGACTCACTATAGGG")
    [-5.4, -5.7, -7.4, -7.5, -7.8, -8.1, -6.1, -6.4, -6.3999999999999995, -6.4, -5.800000000000001, -5.4, -5.1, -5.1, -7.199999999999999, -8.8]
    """
    if not isinstance(sequence, str):
        raise TypeError("sequence  must be a string for an oligonucleotide")
    if len(sequence) < 5:
        raise TypeError(
            "sequence must be a string with more than 5 nucleotides")
    sequence = sequence.upper()
    if not bool(re.compile('^[ACTG]+$').match(sequence)):
        raise ValueError("pentamer must contain only A, C, T, G")

    return [get_delta_g(pentamer) for pentamer in get_oligo_kmers(sequence, 5)]


def get_stability_graph(oligo):
    """
       Return a graph of the oligo's pentamers and their deltaG.

       param: oligo: a string of nucleotides
       return: a matplotlib graph
    """

    oligo = oligo.upper()
    pentamers = get_oligo_kmers(oligo, 5)
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


def main():
    """
       Main function
    """
    oligo = "ACTTGGGATTGGGCT"
    get_stability_graph(oligo)


if __name__ == '__main__':
    main()
