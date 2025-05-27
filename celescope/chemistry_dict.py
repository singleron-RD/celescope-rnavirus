from importlib import resources

chemistry_dict = {
    "auto": {
        "pattern": ""  # auto detect
    },
    "mobiu_5p-1": {
        "pattern": "U8C9L4C9L4C9",
    },
    "mobiu_3p-1": {
        "pattern": "C9L4C9L4C9U8",
    },
    "mobiu-1": {
        "pattern": "C9C9C9U8",  # converted
        "bc": ["bc1.txt", "bc2.txt", "bc3.txt"],
    },
    "mobiu_5p-2": {
        "pattern": "C8L4U9C8L4C8",
    },
    "mobiu_3p-2": {
        "pattern": "C8L4C8U9L4C8",
    },
    "mobiu-2": {
        "pattern": "C8C8C8U9",  # converted
        "bc": ["bc1.txt", "bc2.txt", "bc3.txt"],
    },
    "mobiu_5p-3": {
        "pattern": "C9L6U9C9L4C9",
    },
    "mobiu_3p-3": {
        "pattern": "C9L4C9U9L6C9",
    },
    "mobiu-3": {
        "pattern": "C9C9C9U9",  # converted
        "bc": ["bc1.txt", "bc2.txt", "bc3.txt"],
    },
}

chemistry_dir = str(resources.files("celescope.data.chemistry"))
