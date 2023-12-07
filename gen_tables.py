from itertools import combinations
import sys

sys.setrecursionlimit(1000000)
import math
import os

try:
    import cPickle as pickle
except ImportError:
    import pickle


@staticmethod
def dotProduct(x, y, size):
    """
    comput the dot product of two binary strings:
    dotProduct(0101, 1100, 4) = 0*1 + 1*1 + 0*0 + 1*0 = 1

    Example:
        >>> bin(15)
        '0b1111'
        >>> bin(11)
        '0b1011'
        >>> dotProduct(15, 11, 4)
        1
    """
    a = bin(x)[2:].zfill(size)
    b = bin(y)[2:].zfill(size)
    t = [int(a[i]) & int(b[i]) for i in range(0, size)]
    return reduce((lambda x, y: x ^ y), t) @ staticmethod


class Sbox:
    @staticmethod
    def basic_operate(d, in_var, area_bound, T=0):  #
        '''
        Apply d-depth operations on operands in in_var. d is limited to [0.0, 0.5, 1.0, 1.5, 2.0]
        :param d: depth
        :param in_var: a list of [[expression, depth, area, truth_table], ...]
                        e.g. [['a', 0, 0, '0000000011111111'], ['b', 0, 0, '0000111100001111'], ['c', 0, 0, '0011001100110011'], ['d', 0, 0, '0101010101010101']]
        :param area_bound: area bound
        :param T: the depth that the operand should reach. Default is 0
        :return: a list [[expression, depth, area, truth_table], ...], where expressions are obtained by applying d-depth operations on input expressions.

        Example:
        >>> basic_in_var=[['a', 0, 0, '0000000011111111'], ['b', 0, 0, '0000111100001111'], ['c', 0, 0, '0011001100110011'], ['d', 0, 0, '0101010101010101']]
        >>> Sbox.basic_operate(2, basic_in_var, 50)
        [['( a ) XOR ( b )', 2, 2, '0000111111110000'], \
        ['( a ) XNOR ( b )', 2, 2, '1111000000001111'], \
        ...\
        ['( c ) XOR ( d )', 2, 2, '0110011001100110'], \
        ['( c ) XNOR ( d )', 2, 2, '1001100110011001']]
        '''

        res = []
        assert (d in [0, 0.5, 1, 1.5, 2])
        if d == 0:
            return in_var
        if d == 0.5:  # NOT
            var_ind = list(combinations(range(len(in_var)), 1))
            for i in range(len(var_ind)):
                each = in_var[var_ind[i][0]]
                assert len(each) == 4  # the 4 items are expression, depth, area and value
                if each[1] < T:
                    continue
                new_area = each[2] + 0.5
                if new_area > area_bound:
                    continue
                new_one = []
                new_expr = 'NOT ( ' + each[0] + ' )'
                new_depth = 0.5 + each[1]
                # new_valu = int(not each[3])
                in_valu = each[3]
                new_valu = ''  # new truth table after composition
                for j in range(len(in_valu)):
                    new_valu = new_valu + str(int(not int(in_valu[j])))

                new_one.append(new_expr)
                new_one.append(new_depth)
                new_one.append(new_area)
                new_one.append(new_valu)
                res.append(new_one)
            return res
        if d == 1:  # NAND, NOR operations
            var_ind = list(combinations(range(len(in_var)), 2))  # number of operands is 2
            for i in range(len(var_ind)):
                each0 = in_var[var_ind[i][0]]
                assert len(each0) == 4
                each1 = in_var[var_ind[i][1]]
                assert len(each1) == 4
                if each0[1] < T and each1[1] < T:
                    continue
                new_area = 1 + each0[2] + each1[2]
                if new_area > area_bound:
                    continue

                # NAND
                new_one = []
                new_expr = '( ' + each0[0] + ' ) NAND ( ' + each1[0] + ' )'
                if each0[1] > each1[1]:
                    new_depth = 1 + each0[1]
                else:
                    new_depth = 1 + each1[1]
                # new_valu = int(not (each0[3] & each1[3]))
                in_valu0 = each0[3]
                in_valu1 = each1[3]
                new_valu = ''  # new truth table after composition
                for j in range(len(in_valu0)):
                    new_valu = new_valu + str(int(not int(int(in_valu0[j]) & int(in_valu1[j]))))
                new_one.append(new_expr)
                new_one.append(new_depth)
                new_one.append(new_area)
                new_one.append(new_valu)
                res.append(new_one)

                # NOR
                new_one = []
                new_expr = '( ' + each0[0] + ' ) NOR ( ' + each1[0] + ' )'
                # new_valu = int(not (each0[3] | each1[3]))
                in_valu0 = each0[3]
                in_valu1 = each1[3]
                new_valu = ''  # new truth table after composition
                for j in range(len(in_valu0)):
                    new_valu = new_valu + str(int(not int(int(in_valu0[j]) | int(in_valu1[j]))))
                new_one.append(new_expr)
                new_one.append(new_depth)
                new_one.append(new_area)
                new_one.append(new_valu)
                res.append(new_one)
            return res
        if d == 1.5:  # AND, OR
            var_ind = list(combinations(range(len(in_var)), 2))  # number of operands is 2
            for i in range(len(var_ind)):
                each0 = in_var[var_ind[i][0]]
                assert len(each0) == 4
                each1 = in_var[var_ind[i][1]]
                assert len(each1) == 4
                if each0[1] < T and each1[1] < T:
                    continue
                new_area = 1.5 + each0[2] + each1[2]
                if new_area > area_bound:
                    continue
                # AND
                new_one = []
                new_expr = '( ' + each0[0] + ' ) AND ( ' + each1[0] + ' )'
                if each0[1] > each1[1]:
                    new_depth = 1.5 + each0[1]
                else:
                    new_depth = 1.5 + each1[1]
                # new_valu = int(each0[3] & each1[3])
                in_valu0 = each0[3]
                in_valu1 = each1[3]
                new_valu = ''  # new truth table after composition
                for j in range(len(in_valu0)):
                    new_valu = new_valu + str(int(int(in_valu0[j]) & int(in_valu1[j])))
                new_one.append(new_expr)
                new_one.append(new_depth)
                new_one.append(new_area)
                new_one.append(new_valu)
                res.append(new_one)

                # OR
                new_one = []
                new_expr = '( ' + each0[0] + ' ) OR ( ' + each1[0] + ' )'
                # new_valu = int(each0[3] | each1[3])
                in_valu0 = each0[3]
                in_valu1 = each1[3]
                new_valu = ''  # new truth table after composition
                for j in range(len(in_valu0)):
                    new_valu = new_valu + str(int(int(in_valu0[j]) | int(in_valu1[j])))
                new_one.append(new_expr)
                new_one.append(new_depth)
                new_one.append(new_area)
                new_one.append(new_valu)
                res.append(new_one)
            return res
        if d == 2:  # XOR, XNOR
            var_ind = list(combinations(range(len(in_var)), 2))  # number of operands is 2
            for i in range(len(var_ind)):
                each0 = in_var[var_ind[i][0]]
                assert len(each0) == 4
                each1 = in_var[var_ind[i][1]]
                assert len(each1) == 4
                if each0[1] < T and each1[1] < T:
                    continue
                new_area = 2 + each0[2] + each1[2]
                if new_area > area_bound:
                    continue

                # XOR
                new_one = []
                new_expr = '( ' + each0[0] + ' ) XOR ( ' + each1[0] + ' )'
                if each0[1] > each1[1]:
                    new_depth = 2 + each0[1]
                else:
                    new_depth = 2 + each1[1]
                # new_valu = int(each0[3] ^ each1[3])
                in_valu0 = each0[3]
                in_valu1 = each1[3]
                new_valu = ''  # new truth table after composition
                for j in range(len(in_valu0)):
                    new_valu = new_valu + str( int(in_valu0[j]) ^ int(in_valu1[j]))
                new_one.append(new_expr)
                new_one.append(new_depth)
                new_one.append(new_area)
                new_one.append(new_valu)
                res.append(new_one)

                # XNOR
                new_one = []
                new_expr = '( ' + each0[0] + ' ) XNOR ( ' + each1[0] + ' )'
                # new_valu = int(not (each0[3] ^ each1[3]))
                in_valu0 = each0[3]
                in_valu1 = each1[3]
                new_valu = ''  # new truth table after composition
                for j in range(len(in_valu0)):
                    new_valu = new_valu + str( int(not (int(in_valu0[j]) ^ int(in_valu1[j]))))
                new_one.append(new_expr)
                new_one.append(new_depth)
                new_one.append(new_area)
                new_one.append(new_valu)
                res.append(new_one)
            return res

    @staticmethod
    def full_basic_depth(d, area_bound):
        '''
        Return all expressions with exact depth d. d is limited to [0, 0.5, 1, 1.5, 2]
        :param d: depth
        :param area_bound: area bound of depth d
        :return: a list of expressions with exact depth d.

        Example:
        >>> Sbox.full_basic_depth(2, 50)
        [['( a ) XOR ( b )', 2, 2, '0000111111110000'], \
        ['( a ) XNOR ( b )', 2, 2, '1111000000001111'], \
        ...\
        ['NOT ( NOT ( ( c ) NOR ( d ) ) )', 2.0, 2.0, '1000100010001000']]  ---> not good. There are repeats, so need to reduce.
        '''
        assert d in [0, 0.5, 1, 1.5, 2]
        basic_in_var = [['a', 0, 0, '0000000011111111'], ['b', 0, 0, '0000111100001111'], \
                        ['c', 0, 0, '0011001100110011'], ['d', 0, 0, '0101010101010101']]
        if d == 0.0:
            return Sbox.basic_operate(d, basic_in_var, area_bound, 0)
        elif d == 0.5 or d == 1.0:  # NOT / NAND / NOR operate on basic_in_var
            return Sbox.basic_operate(d, basic_in_var, area_bound, 0)
        elif d == 1.5:
            res = []
            # 1.5-depth operation on basic_in_var
            res = res + Sbox.basic_operate(d, basic_in_var, area_bound)
            # 1-depth operation work on 0.5-depth operands
            res = res + Sbox.basic_operate(1.0, basic_in_var + Sbox.full_basic_depth(0.5, area_bound - 1), area_bound, 0.5)
            # 0.5-depth operation work on 1-depth operands
            res = res + Sbox.basic_operate(0.5, Sbox.full_basic_depth(1, area_bound - 0.5), area_bound, 1)
            return res
        elif d == 2.0:
            res = []
            res = res + Sbox.basic_operate(d, basic_in_var, area_bound)
            res = res + Sbox.basic_operate(1.0,
                                         basic_in_var + Sbox.full_basic_depth(0.5, area_bound - 1) + Sbox.full_basic_depth(
                                             1, area_bound - 1), area_bound, 1)
            res = res + Sbox.basic_operate(1.5, basic_in_var + Sbox.full_basic_depth(0.5, area_bound - 1.5),
                                         area_bound, 0.5)
            res = res + Sbox.basic_operate(0.5, Sbox.full_basic_depth(1.5, area_bound - 0.5), area_bound, 1.5)
            return res


    @staticmethod
    def reduced_depth(depth, area):
        '''
        Return the area-optimal boolean function with the exact given depth within the given area bound
        :param depth:
        :param area:
        :return: a dictionary with truth table being the key, the optimal [expression, depth, area] as the value

        Example:
        >>> Sbox.reduced_depth(2, 50)
        {'0000111111110000': ['( a ) XOR ( b )', 2, 2], \
        '1111000000001111': ['( a ) XNOR ( b )', 2, 2], \
        ... \
        '1110111011101110': ['NOT ( ( c ) AND ( d ) )', 2.0, 2.0]}
        '''
        full_list = Sbox.full_basic_depth(depth, area)

        l = len(full_list)  # there are l expressions with given depth
        res = dict()  # a dictionary to store expressions. The keys are truth table of boolean function.


        for i in range(l):
            ind = full_list[i][3]
            if ind in res:
                if full_list[i][2] < res[ind][2]:  # store the one whose area is lower
                    res[ind][0] = full_list[i][0]  # expression
                    assert res[ind][1] == full_list[i][1]  # depth
                    res[ind][2] = full_list[i][2]  # area
            else:
                res[ind] = [full_list[i][0], full_list[i][1], full_list[i][2]]

        return res


    @staticmethod
    def merge_dict(a, c):
        '''
        merge two dictionaries. For each truth table, keep the one with the lowest depth and area
        :param a:
        :param c:
        :return: merged dictionary
        '''
        b = c.copy()
        for k in a:
            if k in b:
                if a[k][1] < b[k][1]:  # depth
                    b[k] = a[k]
                elif a[k][1] == b[k][1]:
                    if a[k][2] < b[k][2]:  # area
                        b[k] = a[k]
            else:
                b[k] = a[k]

        return b

    @staticmethod
    def dict_operate(d, in_dict, T=0):
        '''
        Apply d-depth operations on operands in in_var in dictionary fashion.
        :param d: depth of operation
        :param in_dict: a dictionary with truth table being the key, [expression, depth, area] being the value
                        e.g. {'0000111111110000': ['( a ) XOR ( b )', 2, 2], \
                            '1111000000001111': ['( a ) XNOR ( b )', 2, 2], \
                            ...,\
                            '1110111011101110': ['NOT ( ( c ) AND ( d ) )', 2.0, 2.0]}
        :param T:
        :return: a dictionary with truth table being the key, [expression, depth, area] being the value

        Example:
        >>> in_dict = Sbox.reduced_depth(2, 50)
        >>> Sbox.dict_operate(0.5, in_dict)
        {'1111000000001111': ['NOT ( ( a ) XOR ( b ) )', 2.5, 2.5], \
        '0000111111110000': ['NOT ( ( a ) XNOR ( b ) )', 2.5, 2.5], \
        ... \
        '0001000100010001': ['NOT ( NOT ( ( c ) AND ( d ) ) )', 2.5, 2.5]}
        '''
        res = dict()
        key = list(in_dict.keys())  # all truth tables of given expressions
        if d == 0:
            return in_dict
        elif d == 0.5:
            var_ind = list(combinations(range(len(in_dict)), 1))

            for i in range(len(var_ind)):
                each = in_dict[key[var_ind[i][0]]]  # the operand
                if each[1] < T:
                    continue
                new_area = each[2] + 0.5

                new_one = []
                new_expr = 'NOT ( ' + each[0] + ' )'
                new_depth = 0.5 + each[1]
                in_str = key[var_ind[i][0]]
                out_str = ''  # new truth table after composition
                for j in range(len(in_str)):
                    out_str = out_str + str(int(not int(in_str[j])))
                if out_str in res:  # store the one with lower area
                    assert new_depth == res[out_str][1]
                    if new_area < res[out_str][2]:
                        res[out_str] = [new_expr, new_depth, new_area]
                else:
                    res[out_str] = [new_expr, new_depth, new_area]
            return res
        elif d == 1:
            var_ind = list(combinations(range(len(in_dict)), 2))  # number of inputs is 2
            for i in range(len(var_ind)):
                each0 = in_dict[key[var_ind[i][0]]]  # the operands
                each1 = in_dict[key[var_ind[i][1]]]
                if each0[1] < T and each1[1] < T:
                    continue
                new_area = 1 + each0[2] + each1[2]

                if each0[1] > each1[1]:
                    new_depth = 1 + each0[1]
                else:
                    new_depth = 1 + each1[1]

                in_str0 = key[var_ind[i][0]]
                in_str1 = key[var_ind[i][1]]
                assert len(in_str0) == len(in_str1)

                new_expr = '( ' + each0[0] + ' ) NAND ( ' + each1[0] + ' )'
                out_str = ''
                for j in range(len(in_str0)):
                    out_str = out_str + str(int(not (int(in_str1[j]) & int(in_str0[j]))))
                if out_str in res:
                    assert new_depth == res[out_str][1]
                    if new_area < res[out_str][2]:
                        res[out_str] = [new_expr, new_depth, new_area]
                else:
                    res[out_str] = [new_expr, new_depth, new_area]

                new_expr = '( ' + each0[0] + ' ) NOR ( ' + each1[0] + ' )'
                out_str = ''
                for j in range(len(in_str0)):
                    out_str = out_str + str(int(not (int(in_str1[j]) | int(in_str0[j]))))
                if out_str in res:
                    assert new_depth == res[out_str][1]
                    if new_area < res[out_str][2]:
                        res[out_str] = [new_expr, new_depth, new_area]
                else:
                    res[out_str] = [new_expr, new_depth, new_area]

            return res
        elif d == 1.5:
            var_ind = list(combinations(range(len(in_dict)), 2))  # number of inputs is 2
            for i in range(len(var_ind)):
                each0 = in_dict[key[var_ind[i][0]]]  # 0:expression,1:depth,2:area
                each1 = in_dict[key[var_ind[i][1]]]
                if each0[1] < T and each1[1] < T:
                    continue
                new_area = 1.5 + each0[2] + each1[2]

                if each0[1] > each1[1]:
                    new_depth = 1.5 + each0[1]
                else:
                    new_depth = 1.5 + each1[1]

                in_str0 = key[var_ind[i][0]]
                in_str1 = key[var_ind[i][1]]
                assert len(in_str0) == len(in_str1)

                new_expr = '( ' + each0[0] + ' ) AND ( ' + each1[0] + ' )'
                out_str = ''
                for j in range(len(in_str0)):
                    out_str = out_str + str(int(int(in_str1[j]) & int(in_str0[j])))
                if out_str in res:
                    assert new_depth == res[out_str][1]
                    if new_area < res[out_str][2]:
                        res[out_str] = [new_expr, new_depth, new_area]
                else:
                    res[out_str] = [new_expr, new_depth, new_area]

                new_expr = '( ' + each0[0] + ' ) OR ( ' + each1[0] + ' )'
                out_str = ''
                for j in range(len(in_str0)):
                    out_str = out_str + str(int(int(in_str1[j]) | int(in_str0[j])))
                if out_str in res:
                    assert new_depth == res[out_str][1]
                    if new_area < res[out_str][2]:
                        res[out_str] = [new_expr, new_depth, new_area]
                else:
                    res[out_str] = [new_expr, new_depth, new_area]

            return res
        elif d == 2:
            var_ind = list(combinations(range(len(in_dict)), 2))  # number of inputs is 2
            for i in range(len(var_ind)):
                each0 = in_dict[key[var_ind[i][0]]]  # 0:expression,1:depth,2:area
                each1 = in_dict[key[var_ind[i][1]]]
                if each0[1] < T and each1[1] < T:
                    continue
                new_area = 2 + each0[2] + each1[2]

                if each0[1] > each1[1]:
                    new_depth = 2 + each0[1]
                else:
                    new_depth = 2 + each1[1]

                in_str0 = key[var_ind[i][0]]
                in_str1 = key[var_ind[i][1]]
                assert len(in_str0) == len(in_str1)

                new_expr = '( ' + each0[0] + ' ) XOR ( ' + each1[0] + ' )'
                out_str = ''
                for j in range(len(in_str0)):
                    out_str = out_str + str(int(int(in_str1[j]) ^ int(in_str0[j])))
                if out_str in res:
                    assert new_depth == res[out_str][1]
                    if new_area < res[out_str][2]:
                        res[out_str] = [new_expr, new_depth, new_area]
                else:
                    res[out_str] = [new_expr, new_depth, new_area]

                new_expr = '( ' + each0[0] + ' ) XNOR ( ' + each1[0] + ' )'
                out_str = ''
                for j in range(len(in_str0)):
                    out_str = out_str + str(int(not (int(in_str1[j]) ^ int(in_str0[j]))))
                if out_str in res:
                    assert new_depth == res[out_str][1]
                    if new_area < res[out_str][2]:
                        res[out_str] = [new_expr, new_depth, new_area]
                else:
                    res[out_str] = [new_expr, new_depth, new_area]

            return res
        else:
            print('wrong')
            return 0
    @staticmethod
    def all_simple_depth(d):
        '''
        Return boolean functions that can be implemented with exact depth d.
        :param d: depth
        :return: a dictionary with truth table being the key, [expression, depth, area] being the value

        Example:
        >>> Sbox.all_simple_depth(1.0)
        {'1111111111110000': ['( a ) NAND ( b )', 1, 1], \
        '1111000000000000': ['( a ) NOR ( b )', 1, 1], \
        ...\
        '1110111011101110': ['( c ) NAND ( d )', 1, 1], \
        '1000100010001000': ['( c ) NOR ( d )', 1, 1]}
        '''
        assert (d in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5])
        file_name = 'new_depth{}.txt'.format(d)
        if os.path.exists(file_name):
            f = open(file_name, 'rb')
            r = pickle.load(f)
            f.close()
            return r

        if d in [0.0, 0.5, 1.0, 1.5, 2.0]:
            res = Sbox.reduced_depth(d, 50)
            f = open(file_name, 'wb')
            pickle.dump(res, f)
            f.close()
            return res

        else:
            # construct by NOT operating on (d-0.5)-depth operands
            res = Sbox.dict_operate(0.5, Sbox.all_simple_depth(d - 0.5), d - 0.5)

            # construct by 1-depth operation operating on (d-1)-depth operands
            r = 2 * d - 2
            in_var = dict()
            for i in [j / 2.0 for j in range(int(r + 1))]:
                in_var = Sbox.merge_dict(in_var, Sbox.all_simple_depth(i))
            res = Sbox.merge_dict(res, Sbox.dict_operate(1.0, in_var, d - 1))

            # construct by 1.5-depth operation operating on (d-1.5)-depth operands
            r = 2 * d - 3
            in_var = dict()
            for i in [j / 2.0 for j in range(int(r + 1))]:
                in_var = Sbox.merge_dict(in_var, Sbox.all_simple_depth(i))
            res = Sbox.merge_dict(res, Sbox.dict_operate(1.5, in_var, d - 1.5))

            # construct by 2-depth operation operating on (d-2)-depth operands
            r = 2 * d - 4
            in_var = dict()
            for i in [j / 2.0 for j in range(int(r + 1))]:
                in_var = Sbox.merge_dict(in_var, Sbox.all_simple_depth(i))
            res = Sbox.merge_dict(res, Sbox.dict_operate(2.0, in_var, d - 2))

            f = open(file_name, 'wb')
            pickle.dump(res, f)
            f.close()

            return res


    @staticmethod
    def str_hw(s):
        res = 0
        assert len(s) == 16
        for i in range(16):
            if s[i] == '1':
                res = res + 1
        return res

    @staticmethod
    def all_simple_depth_hw8(d):
        '''
        Give balanced boolean functions that can be implemented with exact depth d
        :param d: depth
        :return: a dictionary of boolean functions with the truth table being the keys, [expression, depth, area] being the value
        Example:
        >>> Sbox.all_simple_depth_hw8(2)
        {'0000111111110000': ['( a ) XOR ( b )', 2, 2], \
        '1111000000001111': ['( a ) XNOR ( b )', 2, 2], \
        ... \
        '0011001100110011': ['( NOT ( c ) ) NAND ( ( a ) NAND ( c ) )', 2, 2.5], \
        '0101010101010101': ['( NOT ( d ) ) NAND ( ( a ) NAND ( d ) )', 2, 2.5]}
        '''
        res = dict()
        a = Sbox.all_simple_depth(d)
        for k in a:
            if Sbox.str_hw(k) == 8:
                res[k] = a[k]

        return res

    @staticmethod
    def write_bool_func_hw8(d):
        '''
        Generate the file of all balanced boolean functions that can be implemented with exact depth d
        :param d: depth
        :return: a dictionary of boolean functions with the truth table being the keys, [expression, depth, area] being the value
                eg.
        '''
        res = Sbox.all_simple_depth_hw8(d)
        myfile = open('bool_func_hw8_depth' + str(d) + '.txt', 'wb')
        pickle.dump(res, myfile)
        myfile.close()
        return res

    @staticmethod
    def load_expr_hw8_depth(d):
        assert d in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
        filename = 'bool_func_hw8_depth' + str(d) + '.txt'
        if os.path.exists(filename):
            myfile = open('bool_func_hw8_depth' + str(d) + '.txt', 'rb')
            res = pickle.load(myfile)
            myfile.close()
            return res
        else:
            print(filename + ' file has not been generated. Generating...')
            res = Sbox.write_bool_func_hw8(d)
            return res

    @staticmethod
    def write_all_expr():
        '''
        Generate the final file used for Sbox depth evaluation.
        :return: none
        '''
        all_expr = dict()
        all_expr = Sbox.merge_dict(Sbox.load_expr_hw8_depth(0.0), Sbox.load_expr_hw8_depth(0.5))
        all_expr = Sbox.merge_dict(all_expr, Sbox.load_expr_hw8_depth(2.0))
        all_expr = Sbox.merge_dict(all_expr, Sbox.load_expr_hw8_depth(2.5))
        all_expr = Sbox.merge_dict(all_expr, Sbox.load_expr_hw8_depth(3.0))
        all_expr = Sbox.merge_dict(all_expr, Sbox.load_expr_hw8_depth(3.5))
        all_expr = Sbox.merge_dict(all_expr, Sbox.load_expr_hw8_depth(4.0))
        all_expr = Sbox.merge_dict(all_expr, Sbox.load_expr_hw8_depth(4.5))
        myfile = open('all_expr_hw8_within_depth4.5.txt', 'wb')
        pickle.dump(all_expr, myfile)
        myfile.close()
        return

    @staticmethod
    def load_all_expr():
        filename = 'all_expr_hw8_within_depth4.5.txt'
        if os.path.exists(filename):
            myfile = open(filename, 'rb')
            res = pickle.load(myfile)
            myfile.close()
            return res
        else:
            print(filename + ' File has not been generated. Generating...')
            res = Sbox.write_all_expr()
            return res

    @staticmethod
    def find_expression(sbox):
        all_expr = Sbox.load_all_expr()
        key = list(all_expr.keys())
        for i in range(4):
            bl = ''
            for k in range(16):
                bl = bl + str((sbox[k] >> (3 - i)) & 0x1)
            if bl in all_expr:
                if i == 0:
                    print('a\' = ')
                if i == 1:
                    print('b\' = ')
                if i == 2:
                    print('c\' = ')
                if i == 3:
                    print('d\' = ')

                print(all_expr[bl])
            else:
                print('no appear')

        return

    @staticmethod
    def number_of_balanced():
        '''
        Obtain the number of balanced 4-bit boolean functions that can be implemented with each exact depth
        '''
        for d in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]:
            n = len(list(Sbox.all_simple_depth_hw8(d)))
            if n==0:
                print('Number of balanced boolean functions with exact depth {} is {}'.format(d, n))
            else:
                p = math.log(n, 2)
                print('Number of balanced boolean functions with exact depth {} is 2^{}'.format(d, p))
        return



if __name__ == "__main__":
    Sbox.write_all_expr()
    Sbox.find_expression([0xc,0xa,0xd,0x3,0xe,0xb,0xf,0x7,0x8,0x9,0x1,0x5,0x0,0x2,0x4,0x6])
    #Sbox.number_of_balanced()
