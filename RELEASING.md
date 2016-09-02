# Releasing Varcode

This document explains what do Pull Request has been reviewed/finalized and you're
ready to merge your branch into master.

* Assign a version to the release you are preparing. Varcode uses [versioneer], so rather
than explicitly editing `__version__` variable in the code, you must instead tag the branch with your new version number (e.g. `git tag v1.2.3`).

* Once your candidate release branch is tagged you must then run `git push --tags` and merge the branch.

* After the Varcode unit tests complete successfully on Travis then the latest version
of the code (with the version specified above) will be pushed to [PyPI](https://pypi.python.org/pypi) automatically. If you're curious about how automatic deployment is achieved, see our [Travis configuration](https://github.com/hammerlab/varcode/blob/master/.travis.yml#L51).