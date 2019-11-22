#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def style_str(s, color='black', weight=300):
    """
    Style the supplied string with HTML tags

    Args:
        s (str):
            string to format

        color( str ):
            color to show the string in

        weight( int ):
            how thick the string will be displayed

    Returns:
        html(string) : html representation of the string

    """

    return f'<text style="color:{color}; font-weight:{weight}" >{s}</text>'
