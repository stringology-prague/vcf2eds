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
    
    def generate_segment_size(self):
        number = int(math.fabs(random.gauss(self.number, self.number_dev)))
        # print('number - {}'.format(number))
        return number
    
    def generate_element_length(self, max_length=None):
        rnd_length = int(math.fabs(random.gauss(self.length, self.length_dev)))
        if max_length is not None and rnd_length > max_length:
            return max_length

        return rnd_length


def get_degenerate_segment(options, alphabet, weights, curr_depth, avg_settings, max_length=None):
    segment = []
    segment_size = max(2 if options.segment_size_force_min else 1, avg_settings.generate_segment_size())

    for _ in range(segment_size):
        element_len = avg_settings.generate_element_length(max_length)

        # Recursive EDS mode
        if options.reds_prob > 0 and curr_depth < options.max_reds_depth:
            element = ''
            i = 0
            while i < element_len:
                if random.random() < options.reds_prob:
                    data = get_degenerate_segment(options, alphabet, weights, curr_depth + 1,
                                                  avg_settings.create_reduced_copy(), element_len)

                    if len(data) == 1:
                        element += data[0]
                    else:
                        element += '{' + ','.join(data) + '}'

                    i += len(data[0])
                else:
                    element += random.choices(alphabet, weights=weights, k=1)[0]
                    i += 1
        # EDS mode
        else:
            element = ''.join(random.choices(alphabet, weights=weights, k=element_len))

        # Add new element to the segment
        segment.append(element)
    return segment


if __name__ == "__main__":
    parser = OptionParser()

    parser.add_option('-l', '--length', dest="length",
                      default=1000, metavar='NUMBER', type='int',
                      help="Length of EDS - degenerate segments add to length with its first element")

    parser.add_option('--decorate-output', dest="decorateoutput",
                      default=False, action='store_true',
                      help='Append parameters and file type to the output filename')

    # Alphabet options
    parser.add_option('-a', '--alphabet', dest='alphabet',
                      help="Alphabet type (valid values: D for DNA, P for Protein")
    parser.add_option('-f', '--freq-file', dest='freqfile',
                      help="Frequency file of the alphabet symbols (see protein.freq for example)")

    # Structural options
    parser.add_option('-p', '--deg-prob', dest="deg_prob",
                      default=0.15, metavar='NUMBER', type='float',
                      help="Probability of degenerate segment on each position")
    parser.add_option('-r', '--rec-prob', dest="reds_prob",
                      default=0.0, metavar='NUMBER', type='float',
                      help="Probability of starting recursive segment")
    parser.add_option('-b', '--rec-max-depth', dest="max_reds_depth",
                      default=5, metavar='NUMBER', type='int',
                      help="Maximum depth of recursive segments")

    # Segment size options
    parser.add_option('-n', '--segment-size-avg', dest="segment_size_avg",
                      default=4, metavar='NUMBER', type='int',
                      help="Gaussian distribution mean of the number of string in degenerate segment")
    parser.add_option('-m', '--segment-size-stdev', dest="segment_size_stdev",
                      default=0.5, metavar='NUMBER', type='float',
                      help="Gaussian distribution standard deviation of the number of string in degenerate segment")
    parser.add_option('--segment-size-no-force-min', dest="segment_size_force_min",
                      default=True, action='store_false',
                      help="Disables forced minimum segment size of 2 elements. This allows for segments of size 1 to "
                           "be embedded into the parent string, effectively reducing degenerate probability.")

    # Element length options
    parser.add_option('-e', '--element-len-avg', dest="element_len_avg",
                      default=15, metavar='NUMBER', type='int',
                      help="Gaussian distribution mean of the length of degenerate element")
    parser.add_option('-d', '--element-len-stdev', dest="element_len_stdev",
                      default=0.5, metavar='NUMBER', type='float',
                      help="Gaussian distribution standard deviation of the length of degenerate element")

    (options, args) = parser.parse_args()

    print('Random EDS string generator 0.5')

    if options.alphabet and options.freqfile:
        raise Exception('Cannot specify alphabet and weights file at the same time; weights file defines alphabet.')

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

    if options.deg_prob < 0 or options.deg_prob > 1:
        raise Exception('Probability can be only from interval [0, 1] including 0 and 1.')

    output_fname = args[0]
    if options.decorateoutput:
        output_fname = f"{output_fname}_" \
                       f"p={options.deg_prob:0.2f}_" \
                       f"r={options.reds_prob:0.2f}_" \
                       f"s={options.segment_size_avg:06.2f}_" \
                       f"e={options.element_len_avg:06.2f}_" \
                       f"l={options.length:010}.eds"

    print('  EDS length: {}'.format(options.length))
    print('  Alphabet: {} [{}]'.format(options.alphabet, alphabet))
    print('  Weights: {} {}'.format(options.freqfile, weights))
    print('  Degenerate segment probability: {}'.format(options.deg_prob))
    print('  Recursive segment probability: {}'.format(options.reds_prob))
    print('  Maximum recursive depth: {}'.format(options.max_reds_depth))
    print('  Segment size - gauss distribution average: {}'.format(options.segment_size_avg))
    print('  Segment size - gauss distribution standard deviation: {}'.format(options.element_len_stdev))
    print('  Segment size - force minimum segment size of 2: {}'.format(options.segment_size_force_min))
    print('  Element length - gauss distribution average: {}'.format(options.element_len_avg))
    print('  Element length - gauss distribution standard deviation: {}'.format(options.element_len_stdev))
    print('  EDS output file: {}'.format(args[0]))
    print('  EDS output file decorated: {}'.format(output_fname))

    eds_str = ''
    i = 0
    while i < options.length:
        if random.random() < options.deg_prob:
            # degenerate segment
            avg_settings = AverageSettings(options.element_len_avg,
                                           options.element_len_stdev,
                                           options.segment_size_avg,
                                           options.element_len_stdev)
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
