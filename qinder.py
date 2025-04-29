import math

from Bio import SeqIO


def read_file(file_in):
    list_seq, l_header = [], []
    for record in SeqIO.parse(file_in, "fasta"):
        l_header.append(record.id)
        list_seq.append(record.seq)
    return l_header, list_seq


def sequence_score(line, angle):
    g_c_sum_value, a_t_sum_value = 0, 0
    r = 3
    g_ls = []
    c_ls = []
    sample_len = len(line)
    for i, gene in enumerate(line):
        if gene in ('G', 'g'):
            g_val = 0
            if i == 0 or line[i - 1] not in ('G', 'g'):
                for x in range(i, len(line) + 1):
                    if x < len(line) and line[x] in ('G', 'g'):
                        g_val += 1
                    else:
                        g_ls.append(g_val)
                        break
            if i + 1 < len(line):
                if line[i + 1] in ('G', 'g'):
                    g_c_sum_value += r
                if line[i + 1] in ('A', 'a', 'T', 't'):
                    g_c_sum_value += r / 2
                if line[i + 1] in ('A', 'a'):
                    a_t_sum_value += r / 2
                if line[i + 1] in ('T', 't'):
                    a_t_sum_value -= r / 2
        if gene in ('C', 'c'):
            c_val = 0
            if i == 0 or line[i - 1] not in ('C', 'c'):
                for x in range(i, len(line) + 1):
                    if x < len(line) and line[x] in ('C', 'c'):
                        c_val += 1
                    else:
                        c_ls.append(c_val)
                        break
            if i + 1 < len(line):
                if line[i + 1] in ('C', 'c'):
                    g_c_sum_value -= r
                if line[i + 1] in ('A', 'a', 'T', 't'):
                    g_c_sum_value -= r / 2
                if line[i + 1] in ('A', 'a'):
                    a_t_sum_value += r / 2
                if line[i + 1] in ('T', 't'):
                    a_t_sum_value -= r / 2
        if gene in ('A', 'a'):
            if i + 1 < len(line):
                if line[i + 1] in ('A', 'a'):
                    a_t_sum_value += r
                if line[i + 1] in ('G', 'g', 'C', 'c'):
                    a_t_sum_value += r / 2
                if line[i + 1] in ('G', 'g'):
                    g_c_sum_value += r / 2
                if line[i + 1] in ('C', 'c'):
                    g_c_sum_value -= r / 2
        if gene in ('T', 't'):
            if i + 1 < len(line):
                if line[i + 1] in ('T', 't'):
                    a_t_sum_value -= r
                if line[i + 1] in ('G', 'g', 'C', 'c'):
                    a_t_sum_value -= r / 2
                if line[i + 1] in ('G', 'g'):
                    g_c_sum_value += r / 2
                if line[i + 1] in ('C', 'c'):
                    g_c_sum_value -= r / 2

    g_c_value = g_c_sum_value / sample_len
    a_t_value = a_t_sum_value / sample_len

    g4_score_g = math.prod(g_ls) / sample_len
    g4_score_c = math.prod(c_ls) / sample_len

    sin_of_a_t = math.sin(math.radians(angle)) * a_t_value
    final_score = g_c_value + sin_of_a_t

    g4_final_score_g = final_score * g4_score_g
    g4_final_score_c = final_score * g4_score_c

    return line, final_score, g4_final_score_g, g4_final_score_c


class Qinder:

    @staticmethod
    def qinder_app(input_file, window, score, offset, angle, is_all_score):

        limit = window
        offset = offset

        l_header, list_seq = read_file(input_file)

        result_dict_g = dict()
        result_dict_c = dict()

        whole_sequence = list_seq[0]
        for _offset in range(0, int(len(whole_sequence)), offset):
            _sequence, _score, _g4_score_g, _g4_score_c = sequence_score(whole_sequence[_offset:(limit + _offset)],
                                                                         angle)

            if _score > score:
                result_dict_g[str(_sequence)] = _offset, _score, _g4_score_g
            if is_all_score and -score > _score:
                result_dict_c[str(_sequence)] = _offset, _score, _g4_score_c

        return result_dict_g, result_dict_c
