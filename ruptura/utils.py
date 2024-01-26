isothermMeta = {
    'Langmuir': {
        'typeid': 0,
        'labels': ['q_sat', 'b']
    },
    'Anti-Langmuir': {
        'typeid': 1,
        'labels': ['a', 'b']
    },
    'BET': {
        'typeid': 2,
        'labels': ['q_sat', 'b', 'c']
    },
    'Henry': {
        'typeid': 3,
        'labels': ['a']
    },
    'Freundlich': {
        'typeid': 4,
        'labels': ['a', 'nu']
    },
    'Sips': {
        'typeid': 5,
        'labels': ['q_sat', 'b', 'nu']
    },
    'Langmuir-Freundlich': {
        'typeid': 6,
        'labels': ['q_sat', 'b', 'nu']
    },
    'Redlich-Peterson': {
        'typeid': 7,
        'labels': ['a', 'b', 'nu']
    },
    'Toth': {
        'typeid': 8,
        'labels': ['q_sat', 'b', 'nu']
    },
    'Unilan': {
        'typeid': 9,
        'labels': ['q_sat', 'b', 'eta']
    },
    'OBrien-Myers': {
        'typeid': 10,
        'labels': ['q_sat', 'b', 'sigma']
    },
    'Quadratic': {
        'typeid': 11,
        'labels': ['q_sat', 'b', 'c']
    },
    'Temkin': {
        'typeid': 12,
        'labels': ['q_sat', 'b', 'c']
    },
    'BingelWalton': {
        'typeid': 13,
        'labels': ['q_sat', 'a', 'c']
    }
}

pressureScales = {"log": 0, "linear": 1}
markers = ["o", "+", "^", "D", "x", "*", "p", "s", "v"]
getMarker = lambda i: markers[i%len(markers)]
