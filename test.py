import norm_Sean_v20181029 as normalization
import utility_functions
import time
import difflib

normalization.normalize_spectra()

expected_results = "correct_output.txt"
current_results = "log.txt"

with open(expected_results, 'r') as file1:
    with open(current_results, 'r') as file2:
        diff = difflib.unified_diff(
            file1.readlines(),
            file2.readlines(),
            fromfile=expected_results,
            tofile=current_results,
        )
        _empty = object()
        if next(diff, _empty) == _empty:
            print("Test passed, outputs are identical!")
        else:
            for line in diff:
                print(line)


