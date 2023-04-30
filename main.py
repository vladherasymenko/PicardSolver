import re
import numpy as np
import scipy.special
from scipy.misc import derivative
from math import sin, cos, e


def custom_func(func_str):
    # Функція приймає функцію-рядок
    # І повертає функцію пайтона зі строки.
    # Наприклад для "sin(x)" вона поверне функцію, що приймає 1 аргумент і обчислює значення синуса
    def result_func(x):
        # Це функція, що буде повернута
        nonlocal func_str
        # Обчисоюємо значення "функції", заданої рядком
        func_str = func_str.replace('^', '**')
        try:
            return eval(func_str.replace('x', '(' + str(float(x)) + ')'))
        except OverflowError:
            # Якщо отримане значення завелике - повернемо, як результат, 10 в степені 25
            return 10**+25
    return result_func


def integrate_linear(coefs, y_0_i = 0):
    res = [y_0_i]
    for i in range(len(coefs)):
        res.append(coefs[i]/(i+1))
    return res


def list_to_formula(coefs, precision=2):
    res = str(coefs[0])
    for i in range(1, len(coefs)):
        if not np.isclose(coefs[i], 0):
            res+=f"+{round(coefs[i], precision)}*x^{i}"
    return res.replace("+-", "-")


def multiply_polys(poly_list1, poly_list2):
    p1 = np.array(poly_list1)
    term_mult = np.zeros(len(poly_list1)-1)
    for i in range(len(poly_list2)):
        zeros = np.zeros(i)
        term_mult = np.concatenate((term_mult, np.zeros(1)))
        term = p1 * poly_list2[i]
        term_mult += np.concatenate((zeros, term))
    return np.array(term_mult)


def poly_to_power(poly_list, power):
    base = poly_list
    for i in range(int(power)-1):
        poly_list = multiply_polys(poly_list, base)
    return np.array(poly_list)


def add_arr(a, b):
    if len(a) < len(b):
        c = b.copy()
        c[:len(a)] += a
    else:
        c = a.copy()
        c[:len(b)] += b
    return c


def factorial(n):
    if n == 0:
        return 1
    return n*factorial(n-1)


def to_taylor_series(func_str, y_i_estimations, x0, n): # TODO : Add precision
    for key in y_i_estimations.keys():
        func_str = func_str.replace(key, list_to_formula(y_i_estimations[key]))
    f_i = custom_func(func_str)
    result_poly = np.array([f_i(x0)])
    for i in range(1, n+1):
        coef = derivative(f_i, x0, n=i, order=2*i+1)/factorial(i)  # TODO : MB symbolic diffentiation
        result_poly = add_arr(result_poly, coef*poly_to_power([-x0, 1], i))
    return result_poly


def picard_general(func_strings, x0, init_cond, taylor_order=5, num_iters=3):
    # func_strings - list of strings, which contain right hand sides of equations
    # x0 - float - point around which to find solutions
    # init_cond - dict of initial conditions like {"y1" : [1], "y2" : [2]}
    # current_y_estim - dict like init_cond - current estimation of solution
    current_y_estim = init_cond
    for iteration in range(num_iters):
        equation_inx = 0
        prev_y_estim_buffer = current_y_estim
        for y_i in current_y_estim.keys():
            # taylor_poly - list, which represents the polynomial approximation of the current solution
            taylor_poly = to_taylor_series(func_strings[equation_inx], current_y_estim, x0, taylor_order)
            # integrate the polynomial to get next approximation
            y_0_i=init_cond[y_i][0]
            y_i_next = integrate_linear(taylor_poly, y_0_i=y_0_i)
            equation_inx += 1
            prev_y_estim_buffer[y_i] = y_i_next
        current_y_estim = prev_y_estim_buffer
    return current_y_estim # return dictionary - last approximation of solutions - format like init_cond


"""
def string_to_coefs(right_hand_side_strings): # ! ФОРМАТ : -2 * y1^2 + 4 * y2^1 - всі коефіцієнти і всі степені явно виражені !
    coefs = []
    powers = []
    for text in right_hand_side_strings:
        text = text.replace("-", "+-") #TODO : function pre-processing

        coefs_dict = {}
        powers_dict = {}
        pattern = r'([\+\%\*\^])'
        elements = re.split(pattern, text)
        elements.remove('')
        for i in range(len(elements)):
            if elements[i] == '*':
                coefs_dict[elements[i+1]] = float(elements[i-1])
            if elements[i] == '^':
                powers_dict[elements[i-1]] = float(elements[i+1])
        coefs.append(coefs_dict)
        powers.append(powers_dict)
    return coefs, powers

def init_picard(coef_dicts, powers_dict, init_cond, y_i_0 = 1, num_iters = 3, is_non_lin = True):
    y_i_prev = {y_i : [y_i_0] for y_i in init_cond.keys()} # np.array([[y_i_0] for eq in num_equations])
    for iteration in range(num_iters):
        i = 0  # index of current equation
        y_i_prev_buffer = {y_i : [y_i_0] for y_i in init_cond.keys()}
        for y_i in init_cond.keys():  # в теорії це можна паралелізувати - або викорстовувати знайдені наближення
            y_i_next = np.array([0.0])
            for token in coef_dicts[i].keys():
                if token != "x":
                    term = coef_dicts[i][token] * poly_to_power(y_i_prev[token], powers_dict[i][token])
                    y_i_next = add_arr(y_i_next, term)
                else:
                    term = coef_dicts[i][token] * poly_to_power(np.array([0, 1]), powers_dict[i][token])
                    y_i_next = add_arr(y_i_next, term)
            y_i_next = integrate_linear(y_i_next, y_0_i = init_cond[y_i])
            y_i_prev_buffer[y_i] = y_i_next
            print(f"{y_i}_({iteration+1}) = ", list_to_formula(y_i_next))
            i += 1
        y_i_prev = y_i_prev_buffer
        
# example of usage

y1 = 'e^x'
y2 = 'y1'
equations = [y1, y2]
init_cond = {"y1":[1], "y2":[1]}
x0 = 0

picard_general(equations, x0, init_cond, num_iters = 2)
"""

