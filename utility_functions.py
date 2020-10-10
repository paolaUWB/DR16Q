def print_to_file(text: str, file_name: str):
    print(text, file = open(file_name, 'a'))

def file_length(file_name):
    with open(file_name) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def clear_file(file_name: str):
    open(file_name, 'w').close()

def open_file(file_name: str):
    return open(file_name, 'r')