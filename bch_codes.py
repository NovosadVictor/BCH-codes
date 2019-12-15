import random
import numpy as np

from bch_class import GF, BCH

random.seed(1)

q = int(input('Write q: '))
t = int(input('Write t: '))
n = 2 ** q - 1

while t > (n - 1) // 2:
    print('\tt should be <= [(2^q - 1) / 2]\n\tWrite again...\n')
    t = int(input('Write t: '))


field = GF(q, useLUT=0)
coder = BCH(field, t)

print('n = ', coder.n, 'k = ', coder.k, 't = ', 'd = ', 2*t + 1)


def get_random_word(length):
    return [random.choice([0, 1]) for i in range(length)]


all_good = True
for n_errors in range(t + 1):
    print(n_errors)
    for _ in range(100):
        random_word = np.array(get_random_word(coder.k)).astype(int)
        random_code_word = coder.Encode(random_word)

        mask = [1] * n_errors + [0] * (coder.n - n_errors)
        random.shuffle(mask)

        broken_code_word = [random_code_word[i] ^ mask[i] for i in range(len(random_code_word))]

        decoded_word = coder.Decode(np.array(broken_code_word).astype(int))

        if decoded_word.tolist() != random_word.tolist():
            print(random_word, random_code_word, decoded_word, 'WRONG')
            all_good = False

print('All went ', 'GOOD :)' if all_good else 'BAD :C')

if input('Do you want to try? type "yes" - ') == 'yes':
    while True:
        if input('Random word? type "yes" - ') == 'yes':
            word = np.array(get_random_word(coder.k)).astype(int)
        else:
            word = np.array(list(input('Ok, Write word to encode: '))).astype(int)
        if len(word) != coder.k:
            print('World should have length = ', coder.k, '... Type again')
            continue

        print('word: ', ''.join(word.astype(str).tolist()))

        code = coder.Encode(word)
        print('code: ', ''.join(code.astype(str).tolist()))

        if input('Random broken code? type "yes" - ') == 'yes':
            n_errors = random.randint(0, t)
            print('random number of errors: ', n_errors)

            mask = [1] * n_errors + [0] * (coder.n - n_errors)
            random.shuffle(mask)

            broken_code_word = np.array([code[i] ^ mask[i] for i in range(len(code))]).astype(int)
        else:
            broken_code_word = np.array(list(input('Ok, Write broken code: '))).astype(int)

        print('broken code: ', ''.join(broken_code_word.astype(str).tolist()))

        decode = coder.Decode(broken_code_word)
        print('decode: ', ''.join(decode.astype(str).tolist()))
