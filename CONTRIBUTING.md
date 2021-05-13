Contributions are welcomed via [pull requests on GitHub](https://github.com/mosdef-hub/mbuild/pulls). Developers and/or
users will review requested changes and make comments. The rest of this file will serve as a set of general guidelines
for contributors.

# Features

## Implement functionality in a general and flexible fashion

mBuild is designed to be general and flexible, not limited to single chemistries, file formats, simulation engines, or
simulation methods. Additions to core features should attempt to provide something that is applicable to a variety of
use-cases and not targeted at only the focus area of your research. However, some specific features targeted toward
a limited use case may be appropriate. Speak to the developers before writing your code and they will help you make design
choices that allow flexibility.

# Version control

We currently use the "standard" Pull Request model. Contributions should be implemented on feature branches of forks.
Please try to keep the `master` branch of your fork up-to-date with the `master` branch of the main repository.

## Propose a single set of related changes

Small changes are preferred over large changes. A major contribution can often be broken down into smaller PRs. Large PRs that
affect many parts of the codebase can be harder to review and are more likely to cause merge conflicts.

# Source code

## Use a consistent style

It is important to have a consistent style throughout the source code. The following criteria are desired:

* Lines wrapped to 80 characters
* Lines are indented with spaces
* Lines do not end with whitespace
* For other details, refer to [PEP8](https://www.python.org/dev/peps/pep-0008)

We use [pre-commit](https://pre-commit.com/) to automatically check our code style. Pre-commit is included in the dev environment and its git hooks can be installed using:

```bash
pre-commit install
```

## Document code with comments

All public-facing functions should have docstrings using the numpy style. This includes concise paragraph-style description
of what the class or function does, relevant limitations and known issues, and descriptions of arguments. Internal functions
can have simple one-liner docstrings.


# Tests

## Write unit tests

All new functionality in mBuild should be tested with automatic unit tests that execute in a few seconds. These tests
should attempt to cover all options that the user can select. All or most of the added lines of source code should be
covered by unit test(s). We currently use [pytest](https://docs.pytest.org/en/latest/), which can be executed simply by calling
`pytest` from the root directory of the package.
