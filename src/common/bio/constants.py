import tensorflow as tf

ID_TO_SMILES_CHARACTER = {0: ' ',
                          1: '(',
                          2: ')',
                          3: '*',
                          4: '+',
                          5: '-',
                          6: '1',
                          7: '2',
                          8: '3',
                          9: '4',
                          10: '5',
                          11: '=',
                          12: 'C',
                          13: 'H',
                          14: 'N',
                          15: 'O',
                          16: 'P',
                          17: 'S',
                          18: '[',
                          19: ']',
                          20: 'c',
                          21: 'n',
                          22: 'o'}

SMILES_CHARACTER_TO_ID = {' ': 0,
                          '(': 1,
                          ')': 2,
                          '*': 3,
                          '+': 4,
                          '-': 5,
                          '1': 6,
                          '2': 7,
                          '3': 8,
                          '4': 9,
                          '5': 10,
                          '=': 11,
                          'C': 12,
                          'H': 13,
                          'N': 14,
                          'O': 15,
                          'P': 16,
                          'S': 17,
                          '[': 18,
                          ']': 19,
                          'c': 20,
                          'n': 21,
                          'o': 22}

ID_TO_AMINO_ACID = {0: 'A',
                    1: 'T',
                    2: 'G',
                    3: 'C'}

AMINO_ACID_TO_ID = {'A': 0,
                    'T': 1,
                    'G': 2,
                    'C': 3}

NON_STANDARD_AMINO_ACIDS = ['B', 'O', 'U', 'X', 'Z', 'J']


def get_lesk_color_mapping():
    """http://www.bioinformatics.nl/~berndb/aacolour.html
    Args:

    Returns:
         amino acid color mapping
    """
    return tf.cast(tf.constant([
        [0, 0, 0],  # 0 - black - 0
        [255, 255, 255],  # A - White 0,0,0 - 1
        [255, 255, 0],  # C - Yellow 255,255,0 - 2
        [255, 0, 0],  # D - Red - 3

    ]), tf.float32)
