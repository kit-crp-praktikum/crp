# **Setup**

First create a virtual environment outside of this project.

We don't won't to push the dependencies to git and so you would have to install them everytime you switch branches.
```
python3 -m venv /path/to/new/virtual/environment
```

Activate the environment:
```
source /path/to/new/virtual/environment/bin/activate
```

Install the dependencies.
```
pip3 install -r requirements.txt
```

# **Usage**

Remember to activate the environment before you run the script.
```
source /path/to/new/virtual/environment/bin/activate
```

Example usage for partition visualization:
```
python3 vis.py -g ./karlsruhe -p ka_partition -l 4 -c 4 -m 5
```

Example usage for search space visualization:
```
python3 vis.py -g ./germany/ -p ger_c8_l4 -l 4 -c 8 -e query_ger_c8_l4 -w 0.05
```

If you print help, run:
```
python3 vis.py -h
```