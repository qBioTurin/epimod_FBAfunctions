#include <Rcpp.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>

using namespace Rcpp;
using namespace std;

/**
 * @brief Writes optimization model data to a text file.
 * 
 * This function exports a stoichiometry matrix and associated vectors into a text file,
 * formatted specifically for use by optimization software. The file includes reaction identifiers,
 * the stoichiometry matrix, objective coefficients, variable bounds, and additional metadata
 * such as maximum and average biomass estimates for the species being modeled.
 *
 * @param S NumericMatrix, the stoichiometry matrix of the system, where rows represent constraints
 *          and columns represent variables.
 * @param react_id CharacterVector, the identifiers for the reactions corresponding to the columns of S.
 * @param obj_coef NumericVector, the coefficients for the objective function.
 * @param lowbnd NumericVector, the lower bounds for the variables.
 * @param uppbnd NumericVector, the upper bounds for the variables.
 * @param rb NumericVector, the right-hand side coefficients of the constraints, assuming all are equalities.
 * @param model_name string, the base name for the output file, which will be saved as <model_name>.txt.
 * @param write bool, controls whether the file is actually written; if false, the function returns immediately.
 * @param wd string, the working directory where the file will be saved, specifically in a subdirectory 'data'.
 * @param bioMax double, the maximum biomass this model can achieve (default -1 if unspecified).
 * @param bioMean double, the average biomass for this species (default -1 if unspecified).
 *
 * @details
 * The output format is as follows:
 * - Line 1: Space-separated reaction IDs.
 * - Line 2: Two integers representing the dimensions of matrix S and the string "GLP_MAX" indicating a maximization problem.
 * - Line 3: The total number of non-zero entries in the matrix.
 * - Line 4: Biomax value.
 * - Line 5: Biomean value.
 * - Subsequent lines: Each non-zero entry in S formatted as "row_index ; column_index ; value".
 * - Followed by objective coefficients, one per line.
 * - Then, bounds for each constraint (row) formatted as "GLP_FX ; lower_bound ; upper_bound".
 * - Finally, bounds for each variable (column) formatted as "GLP_DB ; lower_bound ; upper_bound".
 *
 * Example:
 * writeModelCpp(S, react_id, obj_coef, lowbnd, uppbnd, rb, "example_model", true, "/my/directory", 0.75, 0.45);
 */

// [[Rcpp::export]]
void writeModelCpp(NumericMatrix S, CharacterVector react_id, NumericVector obj_coef, NumericVector lowbnd, NumericVector uppbnd, NumericVector rb, string model_name, bool write, string wd, double bioMax = -1, double bioMean = -1, double bioMin = -1) {
    if (!write) return;

    int nrow = S.nrow();
    int ncol = S.ncol();

    // Prepare the file path and open the file
    string fileName = wd + model_name + ".txt";
    ofstream file(fileName);
    if (!file.is_open()) {
        Rcout << "Failed to open file: " << fileName << endl;
        return;
    }

    // Buffer for all file content
    stringstream buffer;

    // Write the reaction IDs as a single line
    for (int i = 0; i < react_id.size(); i++) {
        buffer << as<string>(react_id[i]);
        if (i < react_id.size() - 1) buffer << " "; // Add spaces between IDs, not after the last
    }
    buffer << '\n';  // New line after IDs

    // Write the matrix dimensions and problem type
    buffer << nrow << " ; " << ncol << " ; GLP_MAX" << '\n';

    // Initialize the non-zero element counter
    int nnz = 0;
    stringstream matrixData;

    // Use scientific notation for floating point numbers
    matrixData << scientific;

    // Process the matrix for non-zero elements
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            if (S(i, j) != 0) {
                matrixData << (i + 1) << " ; " << (j + 1) << " ; " << S(i, j) << '\n';
                nnz++;
           }
        }
    }

    // Write the number of non-zero elements
    buffer << nnz << '\n';
    
    // Write BioMax, BioMin, BioMean
    buffer << bioMax << '\n';
    buffer << bioMean << '\n';
    buffer << bioMin << '\n';
    

    // Write the objective coefficients
    for (int i = 0; i < obj_coef.size(); i++) {
        buffer << obj_coef[i];
        if (i < obj_coef.size() - 1) buffer << " ";
    }
    buffer << '\n';

    // Write row bounds (constraints)
    for (int i = 0; i < nrow; i++) {
        buffer << "GLP_FX ; " << rb[i] << " ; " << rb[i] << '\n';
    }

    // Write column bounds (variables)
    for (int j = 0; j < ncol; j++) {
        buffer << "GLP_DB ; " << lowbnd[j] << " ; " << uppbnd[j] << '\n';
    }

    // Append matrix data
    buffer << matrixData.str();

    // Write everything from buffer to file
    file << buffer.str();
    file.close();
    Rcout << "File written to: " << fileName << endl;
}


