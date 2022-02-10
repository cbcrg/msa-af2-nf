#!/usr/bin/env python3

import json
import sys

f = open(sys.argv[1],)
data = json.load(f)
print(max(list(data['plddts'].values())))
f.close()
