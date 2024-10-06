
def delete_fourth_fifth_col(file):
    line_count = 0

    with open(file, "w") as file:
        for line in file:
            if line != "\n":
                line_count += 1

            arr = line.split()
            for k in range(len(arr)):
                arr[3] = " "
                arr[4] = " "
    file.close()


atr_csem = 'atr_csem_inversion.txt'
delete_fourth_fifth_col(atr_csem)

