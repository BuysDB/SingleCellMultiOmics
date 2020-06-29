import re

# From https://github.com/django/django/blob/master/django/utils/text.py, then converted to use a precomiled pattern
file_replace_pattern = re.compile(r'(?u)[^-\w.]')
def get_valid_filename(s):
    """
    Return the given string converted to a string that can be used for a clean
    filename. Remove leading and trailing spaces; convert other spaces to
    underscores; and remove anything that is not an alphanumeric, dash,
    underscore, or dot.
    >>> get_valid_filename("john's portrait in 2004.jpg")
    'johns_portrait_in_2004.jpg'
    """
    s = str(s).strip().replace(' ', '_')
    return file_replace_pattern.sub( '', s)
