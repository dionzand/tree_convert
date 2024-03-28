# check isogg_tree.json to see if each child only belongs to one parent

import json

isogg_tree = json.load(open("isogg_tree.json"))

child_to_parent = {}

for parent, children in isogg_tree.items():
    for child in children:
        if child in child_to_parent:
            print(f"{child} belongs to {child_to_parent[child]} and {parent}")
        child_to_parent[child] = parent

# check if each parent is also a child of another parent
for parent, children in isogg_tree.items():
    if parent not in child_to_parent:
        print(f"{parent} is not a child of any parent")