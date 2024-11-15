def write_dict_to_file(dictionary, filename):
    with open(filename, 'w') as file:
        for key, value in dictionary.items():
            file.write(f"{key}:{value}\n")


def read_dict_from_file(filename):
    dictionary = {}
    with open(filename, 'r') as file:
        for line in file:
            key, value = line.strip().split(':', 1)  # Split on first colon
            dictionary[key] = value
    return dictionary


def file_name_to_unix(f):
    return f.replace('/mnt/g', 'G:').replace("/", chr(92))
