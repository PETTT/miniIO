# Readme

Instructions and notes when using `libautotuner.so`.

## Vagrant

The terminal app doesn't work out of the box in the vagrant environment.
See https://bbs.archlinux.org/viewtopic.php?id=180103 for a fix.

The VM should be restarted after initial creation.

## Links

* https://github.com/HDFGroup/H5Tuner
* https://github.com/HDFGroup/H5Tuner/wiki


## Environmental Variables

* `LD_PRELOAD` - Set to libautotuner.
* `H5TUNER_VERBOSE` - Set 0 to 4 (least to most) for output when using libautotuner.
* `H5TUNER_CONFIG_FILE` - Location of the h5tuner config file.
