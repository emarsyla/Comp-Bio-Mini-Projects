# Sample lists of strings
list_of_lists = [
    ['G,G,G,G'],
    ['A,B,C,D'],
    ['E,F,G,H']
]

# Initialize a list of empty lists to collect characters for each position
max_length = max(len(lst) for lst in list_of_lists)  # Find the maximum number of strings in any list
columns = [[] for _ in range(max_length)]

# Iterate through each list of strings
for list_of_strings in list_of_lists:
    # Iterate through each string in the current list
    for index, string in enumerate(list_of_strings):
        # Iterate through each character in the current string
        for char in string:
            # Skip commas and add characters to the appropriate column
            if char != ',':
                columns[index].append(char)

# Convert lists in columns to tuples
result_tuples = [tuple(col) for col in columns]

print(result_tuples)
