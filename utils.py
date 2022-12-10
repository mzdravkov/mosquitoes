from prettytable import PrettyTable

def print_table(data, columns=None):
    """
    Prints a table with the given data.
    """
    table = PrettyTable()
    if columns:
        table.field_names = columns
    else:
        table.field_names = data[0]
    rows = data if columns else data[1:]
    for row in rows:
        table.add_row(row)
    print(table) 