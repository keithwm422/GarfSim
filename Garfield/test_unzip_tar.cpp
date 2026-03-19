#include <iostream>
#include <string>
#include <cstdlib> // For system()
#include <fstream>
void extractFileFromTarGz(const std::string& archivePath, const std::string& fileNameToExtract, const std::string& outputPath) {
    // Construct the command to extract the specific file
    std::string command = "tar -xzvf " + archivePath + " " + " -C " + outputPath + " " + fileNameToExtract; //tar -xzvf HEATModelGarfieldFiles.tar.gz -C /home/kmcbride/HEATModel_files_for>    std::cout << "string command is " << command << std::endl;
    // Execute the command
    int result = std::system(command.c_str());

    if (result == 0) {
        std::cout << "Successfully extracted '" << fileNameToExtract << "' to '" << outputPath << "'." << std::endl;
    } else {
        std::cerr << "Error extracting file: " << result << std::endl;
    }
}

int main() {
    std::string archive = "HEATModelGarfieldFiles.tar.gz";
    std::string fileToExtract = "HEATModelForGarfield/HEATModel_xslice_0_column_1.csv"; // Path within the archive // looks like HEATModelForGarfield/HEATModel_xslice_0_column_0.csv
    std::string outputDir = "/home/kmcbride/garfield/isoHEATB_codes/GarfSim/Garfield/HEATModel_files_for_garfield";

    extractFileFromTarGz(archive, fileToExtract, outputDir);

    // Now you can read the extracted file:
    std::string extractedFilePath = outputDir + "/" + fileToExtract; // Adjust if file path within archive differs from output path

    return 0;
}
