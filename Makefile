# Subdirectories to recurse in build order
SUBDIR = rnascoring gtmfe iB4e bin

.PHONY: $(SUBDIR) recurse

$(MAKECMDGOALS) recurse: $(SUBDIR)

$(SUBDIR):
	$(MAKE) -C $@ $(MAKECMDGOALS)
