## Usage Examples

Here are example code snippets demonstrating the primary usage of the main classes and methods from your project—based on their direct use in the main program, including setup, optimization, file output, and plotting.

---

### 1. Initialize Input Data from Command Line Arguments

```cpp
#include "input.h"

// Parse input data and options from command-line arguments
int main(int argc, char* argv[]) {
    input Input(argc, argv);

    // Typical workflow might continue with computation and access to parsed data:
    // std::map<std::string, bool> options;
    // Input.compute(options);
    // auto data = Input.GetData();
    // auto interpolation = Input.GetInterpolation();
    return 0;
}
```
The `input` class parses user input and data files as described in your README.

---

### 2. Automatic Research for Best Parameters

```cpp
#include "AutomaticResearch.h"
#include <vector>
#include <map>

// Suppose par_best is a vector of optimal parameters and options["complex"] is a boolean
void run_auto_research(const std::vector<double>& par_best, const std::map<std::string, bool>& options) {
    AutomaticResearch Auto(par_best, options.at("complex"));
    // Auto runs parameter refinement automatically, using multi-threading if enabled.
}
```

---

### 3. Chi-Square Minimization

```cpp
#include "ChiSquareMinimizer.h"
#include <map>

// Suppose options is a map containing various optimization flags
void optimize_chi_square(const std::map<std::string, bool>& options) {
    ChiSquareMinimizer Optimizer(options);
    // Optimizer.minimize(...); // See class for available methods
}
```
This class helps perform parameter minimization based on the chi-square metric.

---

### 4. Write Results to File

```cpp
#include "file.h"
#include "input.h"

// Given input::Data data, input::Interpolation interpolation, and an 'approx' flag
void save_results(const input::Data& data, const input::Interpolation& interpolation, bool approx) {
    writeFile(data, interpolation, approx);
}
```
This function writes the results and fitting parameters to an output file.

---

### 5. Generate Plots

```cpp
#include "PlotGenerator.h"
#include "input.h"

// Given input::Data and input::Interpolation
void generate_plot(const input::Data& data, const input::Interpolation& interpolation) {
    PlotGenerator p_generator(data, interpolation);
    // p_generator.generate(); // Or another method, depending on implementation
}
```
This class handles the visualization of your interpolation results.

---

For more advanced examples, see main program sources and header files such as:
- [input.h](https://github.com/NicolaLavarda/Interpolazione_multiparametrica/blob/main/Interpolazione_multiparametrica/include/input.h)
- [file.h](https://github.com/NicolaLavarda/Interpolazione_multiparametrica/blob/main/Interpolazione_multiparametrica/include/file.h)

These examples show the core workflow and can be adapted for your main analysis script.
