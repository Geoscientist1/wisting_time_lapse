import numpy as np
from itertools import product


def four_thought():
    operators = [' * ', ' + ', ' - ', ' // ']
    output = {}

    for itr in product(operators, repeat=3):
        exp1 = f"4 {itr[0]} 4 {itr[1]} 4 {itr[2]} 4"
        exp2 = f"(4 {itr[0]} 4) {itr[1]} (4 {itr[2]} 4)"
        exp3 = f"(4 {itr[0]} 4 {itr[1]} 4) {itr[2]} 4"
        exp4 = f"4 {itr[0]} (4 {itr[1]} 4 {itr[2]} 4)"
        exp5 = f"4 {itr[0]} 4 {itr[1]} (4 {itr[2]} 4)"
        exp6 = f"4 {itr[0]} (4 {itr[1]} 4) {itr[2]} 4"
        exp7 = f"((4 {itr[0]} 4) {itr[1]} 4) {itr[2]} 4"

        all_express = [exp1, exp2, exp3, exp4, exp5, exp6, exp7]

        for i in all_express:
            try:
                result = eval(i)
                if result not in output:
                        output[result] = i.replace('//', '/')
            except ZeroDivisionError:
                continue

    return output


def main():
    with open('4thought.txt', "r") as file:
        read_data = []
        for line in file:
            read_data.append(line)

    num_of_tests = read_data[0]

    integer_val = [int(read_data[i]) for i in range(1, int(num_of_tests) + 1)]

    expressions = four_thought()

    results = []
    for case in integer_val:
        if case in expressions:
            results.append(f"{expressions[case]} = {case}")
        else:
            results.append("no solution")

    for result in results:
        print(result)


if __name__ == "__main__":
    main()
