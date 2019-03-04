Contributions are welcomed via [pull requests on GitHub](https://github.com/mosdef-hub/mbuild/pulls). Developers and/or 
users will review requested changes and make comments. This The rest of this file will serve as a set of general guidelines
for contributors.

# Features

## Implement functionality in a general and flexible fashion

mBuild is designed to be general and flexible, not limited to single chemistries, file formats, simulations engines, or 
simulation methods. Additions to core features shoudl attempt to provide something that is applicable to a variety of 
use-cases and not  just the one thing you might need it to do for your research. But some specific features targetted toward 
a limited use case may be appropriate. Speak to the developers before writing your code and they will help you make design 
choices that allow flexibility.

# Version control

We currently use the "standard" Pull Request model. Contributions should be implemented on feature branches of forks.
Please try to keep the `master` branch of your fork up-to-date with the `master` branch of the main repository.

## Propose a single set of related changes

Small changes are preferred over large changes. A major contribution can be broken down into smaller PRs. Large PRs that 
affect many parts of the codebase can be harder to review and cause merge conflicts if they drag on for long periods of time. 

# Source code

## Use a consistent style

It is important to have a consistent style throughout the source code. The following criteria are recommended (but strict
adherence is not necessary):

* Lines to 80 characters
* General consistent with [PEP8](https://www.python.org/dev/peps/pep-0008)
* The use of linters such as [Black](https://github.com/ambv/black) or [pycodestyle](https://github.com/PyCQA/pycodestyle) 

## Document code with comments

All public-facing functions should have doscstrings using the numpy style. This includes concise paragraph-style description
of what the class or function does, relevant limitations and known issues, and descriptions of arguments. Internal functions
can have simple one-liner docstrings.


# Tests

## Write unit tests

All new functionality in mBuild should be tested with automatic unit tests that execute in a few seconds. These tests 
should attempt to cover all options that the user can select. All or most of the added lines of source code should be 
covered by unit test(s).
