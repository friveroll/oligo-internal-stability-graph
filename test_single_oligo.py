import pytest
import single_oligo


@pytest.mark.parametrize('sequence, k, expected_result', [
    ('ACTTG', 2, ['AC', 'CT', 'TT', 'TG']),
    ('ACTTG', 3, ['ACT', 'CTT', 'TTG']),
    ('AATAAG', 2, ['AA', 'AT', 'TA', 'AA', 'AG']),
    ('AATAAG', 3, ['AAT', 'ATA', 'TAA', 'AAG']),
    ('AATAAG', 5, ['AATAA', 'ATAAG']),
])
def test_get_oligo_kmers(sequence, k, expected_result):
    assert single_oligo.get_oligo_kmers(sequence, k) == expected_result


def test_get_oligo_kmers_raises_type_error():
    with pytest.raises(TypeError):
        single_oligo.get_oligo_kmers(123, 2)


@pytest.mark.parametrize('sequence, k, message', [
    ('ACTTG', 1, 'k must be an integer greater than 1'),
    ('ACTTG', 6, 'k must be less than the length of the sequence'),
])
@pytest.mark.xfail(raises=ValueError)
def test_get_oligo_kmers_raises_value_error(sequence, k, message):
    assert single_oligo.get_oligo_kmers(sequence, k) == message


@pytest.mark.parametrize('sequence, expected_result, dimer_list', [
    ('ACTTG', -6.7, ['AC', 'CT', 'TT', 'TG']),
    ('TCTTGT', -8.3, ['TC', 'CT', 'TT', 'TG', 'GT']),
    ('ACTTGGGATTGGGCT', -24.5, ['AC', 'CT', 'TT',
                                'TG', 'GG', 'GA',
                                'AT', 'TT', 'TG',
                                'GG', 'GC', 'CT']),
])
def test_get_delta_g(mocker, sequence, dimer_list, expected_result):
    mocker.patch('single_oligo.get_oligo_kmers', return_value=dimer_list)
    assert single_oligo.get_delta_g(sequence) == expected_result


def test_get_delta_g_raises_type_error():
    with pytest.raises(TypeError):
        single_oligo.get_delta_g(123)


def test_get_delta_g_raises_value_error():
    with pytest.raises(ValueError):
        single_oligo.get_delta_g('HIJKLM')


@pytest.mark.parametrize('sequence, expected_result', [
    ('ACTTG', [-6.7]),
    ('TCTTGT', [-7.0, -6.7]),
    ('ACTTGGGATTGGGCT', [-6.7, -8.5, -10.0,
                         -9.7, -9.3, -8.1,
                         -6.9, -8.4, -10.0,
                         -11.2, -10.9]),
])
def test_get_delta_g_pentamers(sequence, expected_result):
    assert single_oligo.get_delta_g_pentamers(sequence) == expected_result
