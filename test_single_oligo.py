import pytest
import single_oligo


def test_get_pentamers():
    oligo = 'ATGCGT'
    pentamers = single_oligo.get_pentamers(oligo)
    assert pentamers == ['ATGCG', 'TGCGT']


def test_get_two_mers():
    oligo = 'ATGCGT'
    two_mers = single_oligo.get_two_mers(oligo)
    assert two_mers == {'AT': 1, 'TG': 1, 'GC': 1, 'CG': 1, 'GT': 1}


def test_get_delta_g_pentamer():
    pentamer = 'AATAA'
    delta_g = single_oligo.get_delta_g_pentamer(pentamer)
    assert delta_g == -6.3


def test_get_delta_g_pentamers():
    oligo = 'ATGCGT'
    delta_g = single_oligo.get_delta_g_pentamers(oligo)
    assert delta_g == [-10.1, -9.9]

