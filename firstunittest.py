"""Example functions"""


def test1():
    return 1

def calculate_fibonacci_number(n):
    """Calculates the fibonacci number on index n. Index starts at 0 with value 0."""
    
    result = 0
    if n == 0:
        return 0
    elif n == 1:
        return 1 
    else:
        result = calculate_fibonacci_number(n-1) + calculate_fibonacci_number(n-2)
    return result