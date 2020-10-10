import norm_Sean_v20181029 as normalization
import utility_functions
import time
import difflib

normalization.normilize_spectra()

expected_results = "correct_output.txt"
current_results = "log.txt"

with open(expected_results, 'r') as hosts0:
    with open(current_results, 'r') as hosts1:
        diff = difflib.unified_diff(
            hosts0.readlines(),
            hosts1.readlines(),
            fromfile=expected_results,
            tofile=current_results,
        )
        _empty = object()
        if next(diff, _empty) == _empty:
            print("Test passed, outputs are identical!")
        else:
            for line in diff:
                print(line)


