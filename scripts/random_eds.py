from optparse import OptionParser
import csv
import random
import math

DNA_ALPHABET = 'ACGT'
PROTEIN_ALPHABET = 'ACDEFGHIKLMNPQRSTVWY'


class AverageSettings:
    def __init__(self, length, length_dev, number, number_dev):
        self.length = length
        self.length_dev = length_dev
        self.number = number
        self.number_dev = number_dev

    def create_reduced_copy(self):
        return AverageSettings(math.log2(self.length) if self.length >= 1 else self.length,
                                self.length_dev,
                                math.log2(self.number) if self.number >= 1 else self.number,
                                self.number_dev)
    
    def get_number(self):
        number = int(math.fabs(random.gauss(self.number, self.number_dev)))
        # print('number - {}'.format(number))
        return number
    
    def get_length(self, max_length=None):
        rnd_length = int(math.fabs(random.gauss(self.length, self.length_dev)))
        if max_length is not None and rnd_length > max_length:
            return max_length

        return rnd_length


def get_random_string(alphabet, weights, length):
    return ''.join(random.choices(alphabet, weights=weights, k=length))


def get_degenerate_segment(options, alphabet, weights, curr_depth, avg_settings, max_length=None):
    # print('----- RECURSIVE CALL ------')
    if curr_depth == options.max_reds_depth:
        return [get_random_string(alphabet, weights, avg_settings.get_length(max_length))]

    degenerate_segment = []

    number_of_items = avg_settings.get_number()
    for _ in range(number_of_items):
        segment_length = avg_settings.get_length(max_length)

        segment_str = ''
        i = 0
        while i < segment_length:
            if random.random() < options.reds_prob:
                data = get_degenerate_segment(options, alphabet, weights, curr_depth + 1, avg_settings.create_reduced_copy(), segment_length)

                if len(data) == 1:
                    segment_str += data[0]
                else:
                    segment_str += '{' + ','.join(data) + '}'

                i += len(data[0])
            else:
                segment_str += random.choices(alphabet, weights=weights, k=1)[0]
                i += 1

        degenerate_segment.append(segment_str)

    if len(degenerate_segment) == 0:
        degenerate_segment.append(get_random_string(alphabet, weights, avg_settings.get_length(max_length)))
    
    # print('----- END ------')
    return degenerate_segment


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option('-a', '--alphabet', dest='alphabet',
                    help="Alphabet type (valid values: D for DNA, P for Protein")
    parser.add_option('-f', '--freq-file', dest='freqfile',
                    help="Frequency file of the alphabet symbols (see protein.freq for example)")
    parser.add_option('-p', '--probability', dest="probability", default=0.15, metavar='NUMBER',
                    type='float', help="Probability of degenerate segment on each position")
    parser.add_option('-e', '--average', dest="average", default=15, metavar='NUMBER',
                    type='int', help="Gaussian distribution mean of the length of degenerate segment")
    parser.add_option('-d', '--deviation', dest="deviation", default=0.5, metavar='NUMBER',
                    type='float', help="Gaussian distribution standard deviation of the length of degenerate segment")
    parser.add_option('-n', '--number', dest="number", default=4, metavar='NUMBER',
                    type='int', help="Gaussian distribution mean of the number of string in degenerate segment")
    parser.add_option('-m', '--ndeviation', dest="number_deviation", default=0.5, metavar='NUMBER',
                    type='float', help="Gaussian distribution standard deviation of the number of string in degenerate segment")
    parser.add_option('-r', '--recursive', dest="reds_prob", default=0.1, metavar='NUMBER',
                    type='float', help="Probability of starting recursive segment")
    parser.add_option('-b', '--maxdepth', dest="max_reds_depth", default=5, metavar='NUMBER',
                    type='int', help="Maximum depth of recursive segments")
    parser.add_option('-l', '--length', dest="length", default=1000, metavar='NUMBER',
                    type='int', help="Length of EDS - degenerate segments add to length with its first element")
    parser.add_option('--decorate-output', dest="decorateoutput", default=False, action='store_true',
                    help='Append parameters and file type to the output filename')

    (options, args) = parser.parse_args()

    print('Random EDS string generator 0.5')

    alphabet = DNA_ALPHABET
    if options.alphabet == 'P':
        alphabet = PROTEIN_ALPHABET

    weights = None
    if options.freqfile:
        with open(options.freqfile) as f:
            reader = csv.reader(f, delimiter=" ", skipinitialspace=True)
            freqs = [r for r in reader]
            alphabet = [f[1] for f in freqs]
            sum_weights = sum([int(f[0]) for f in freqs])
            weights = [float(f[0]) / sum_weights for f in freqs]

    if options.probability < 0 or options.probability > 1:
        raise Exception('Probability can be only from interval [0, 1] including 0 and 1.')

    output_fname = args[0]
    if options.decorateoutput:
        output_fname = f"{output_fname}_p={options.probability:0.2f}_e={options.average:06.2f}_n={options.number:06.2f}_r={options.reds_prob:0.2f}_l={options.length:010}.eds"

    print('  alphabet: {} [{}]'.format(options.alphabet, alphabet))
    print('  weights: {} {}'.format(options.freqfile, weights))
    print('  Degenerate segment pp: {}'.format(options.probability))
    print('  Gauss avg. length: {}'.format(options.average))
    print('  Gauss standard deviation length: {}'.format(options.deviation))
    print('  Gauss avg. number: {}'.format(options.number))
    print('  Gauss standard deviation number: {}'.format(options.number_deviation))
    print('  REDS segment pp: {}'.format(options.reds_prob))
    print('  Maximum REDS depth: {}'.format(options.max_reds_depth))
    print('  EDS length: {}'.format(options.length))
    print('  EDS output file: {}'.format(args[0]))
    print('  EDS output file decorated: {}'.format(output_fname))

    eds_str = ''
    i = 0
    while i < options.length:
        if random.random() < options.probability:
            # degenerate segment
            avg_settings = AverageSettings(options.average, options.deviation, options.number, options.number_deviation)
            degenerate_segment = get_degenerate_segment(options, alphabet, weights, 0, avg_settings)

            if len(degenerate_segment) == 1:
                eds_str += degenerate_segment[0]
            else:
                eds_str += '{' + ','.join(degenerate_segment) + '}'

            i += len(degenerate_segment[0])
        else:
            eds_str += random.choices(alphabet, weights=weights, k=1)[0]
            i += 1


    with open(output_fname, "w") as f:
        f.write(eds_str)
