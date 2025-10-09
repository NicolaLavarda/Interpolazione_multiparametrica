# Interpolazione_multiparametrica

**Generic multi-parameter interpolation from file data**

## Overview

Interpolazione_multiparametrica is a C++ library for performing multi-parameter interpolation using data loaded from files. It is designed to be flexible and efficient, allowing users to interpolate values based on multiple input variables and supporting various use cases in engineering, scientific computing, and data analysis.

## Features

- **Multi-parameter interpolation:** Supports interpolation across any number of input dimensions.
- **File-based data input:** Easily load datasets from files for interpolation.
- **Generic design:** Adaptable to different file formats and data structures.
- **High performance:** Implemented in modern C++ for speed and reliability (also through multi-threading).
- **Simple integration:** Minimal dependencies for easy inclusion in other projects.

## Prerequisites

- **CERN ROOT:** You must have [ROOT](https://root.cern/) (the data analysis framework developed at CERN) installed and configured on your system to compile and use this project.
- A C++ compiler supporting C++11 or newer.
- `make` utility for building (optional, if using the provided Makefile).

## Building

Clone the repository:
```bash
git clone https://github.com/NicolaLavarda/Interpolazione_multiparametrica.git
cd Interpolazione_multiparametrica
```

Build using `make`:
```bash
make
```

> **Note:** The build process requires the ROOT environment to be set up. Make sure you've sourced the ROOT setup script, e.g.:
> ```bash
> source /path/to/root/bin/thisroot.sh
> ```
> Adjust the path to your specific ROOT installation.

## Usage

```
 Usage: ./interpolate <file.txt> <param1> <param2> ... <interpolating_function> [options]

Example:
  ./interpolate file.txt 1.5 a a "a*sin(x*b)+c" plot improve

Arguments:
  <file.txt>               - Input data file containing values to interpolate.
  <param1> <param2> ...    - Initial values for parameters in the interpolation function.
                             Use 'a' to indicate parameters that should be determined automatically.
  <interpolating_function> - Mathematical function to be used for interpolation.
                             This must be enclosed in double quotes ("").

Optional Flags:
  improve  - Uses previously computed parameters from the input file as starting values.
  faster   - Speeds up the interpolation process.
  approx   - Rounds results to the appropriate significant figures.
  complex  - Displays intermediate steps in the chi-squared minimization.
  save     - Saves results to the input file only if the chi-squared improves or is not present.
  save!    - Forces saving of results to the input file regardless of improvement.
  plot     - Generates a graphical representation of the interpolation.

Notes:
- Ensure the function syntax follows mathematical conventions and includes parameters defined in input.
- The names of the parameters used in the interpolating function must be written in alphabetical order (e.g., if using three parameters, they must be 'a', 'b', 'c').
- Flags should be placed at the end of the command.
- Use 'help' as the only argument to display this message.
```

### Data File Format

Accepted formats for the input text file:
- **Four columns:** `x`, `sigma_x`, `y`, `sigma_y`
- **Three columns:** `x`, `y`, `sigma_y`

Columns should be separated by spaces or tabs.

**Data File Examples:**
```
0.208   0.0032   2.42   0.0352
0.306   0.0041   3.62   0.0478
...
```
or
```
0.208   2.42   0.0352
0.306   3.62   0.0478
...
```
If your file has only three columns, the third is interpreted as the statistic error on y.

### Other Usage Examples

The key classes and functions typically used in the main program and same exemples can be consulted in the file [Usage_Examples.md](https://github.com/NicolaLavarda/Interpolazione_multiparametrica/blob/main/Usage_Examples.md)

## Language

This project is primarily written in **C++** (99.9%), with a minimal amount of **Makefile** (0.1%) for build automation.

## Contributing

Feel free to open issues or submit pull requests for improvements and bug fixes.

## License

This project is open source. See the [LICENSE](LICENSE) file for details.

## Author

Nicola Lavarda  
For questions or suggestions, open an issue on GitHub.
