import unicodedata
import re

def normalize(s):
    s = replace_some_greek_letters(s)
    s = replace_non_alnum_intelligent(s)
    return s

def replace_some_greek_letters(s):
    s = s.replace('alpha', 'a')
    s = s.replace('beta', 'b')
    s = s.replace('gamma', 'g')
    s = s.replace('α', 'a')
    s = s.replace('𝛼', 'a')
    s = s.replace('β', 'b')
    s = s.replace('γ', 'g')
    return ''.join(c for c in unicodedata.normalize('NFD', s)
                   if unicodedata.category(c) != 'Mn')

def replace_non_alnum_intelligent(line):
    ret = []
    for i in range(len(line)):
        char = line[i]
        if re.search(r'[(){}\[\]!.;,*/<>]', char):
            before = line[i - 1:i + 1]
            if before != 'c.' and before != 'p.':
                ret.append(' ')
                ret.append(char)
                ret.append(' ')
            else:
                ret.append(char)
        else:
            ret.append(char)
    ret = ''.join(ret)
    oldret = ret
    while oldret != ret:
        oldret = ret
        ret = ret.replace('  ', ' ')
    return ret
