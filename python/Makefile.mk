$(top_srcdir)/python/dist/gridlabd-$(notdir $(prefix))-.tar.gz: source/build.h
	@echo "building $@"
	@python3 -m pip install build 1>/dev/null
	@( export SRCDIR=$(realpath $(top_srcdir)) ; python3 -m build $(top_srcdir)/python 1>/dev/null )

python-install: $(top_srcdir)/python/dist/gridlabd-$(notdir $(prefix))-.tar.gz
	@echo "installing $@"
	@python3 -m pip install --ignore-installed $(top_srcdir)/python/dist/gridlabd-*.whl

python-clean:
	@rm $(top_srcdir)/python/dist/*.{whl,tar.gz}
	@echo "uninstalling $(top_srcdir)/python"
	@python3 -m pip uninstall gridlabd -y || (echo "Use '. utilities/cleanwc' instead to clean this build."; exit 1)
