// memory_usage_macro.C
void memory_usage_macro(const char* filename = "your_file.root") {
    TFile* file = TFile::Open(filename, "READ");
    if (!file || file->IsZombie()) {
        printf("Error: Could not open file %s\n", filename);
        return;
    }

    printf("Analyzing objects in file: %s\n", filename);
    printf("--------------------------------------------------\n");
    printf("%-40s %-15s\n", "Object Name", "Size (bytes)");
    printf("--------------------------------------------------\n");

    // Get the list of keys in the file
    TList* keys = file->GetListOfKeys();
    if (!keys) {
        printf("No objects found in the file.\n");
        file->Close();
        delete file;
        return;
    }

    // Iterate through the keys and print information
    TIter next(keys);
    TKey* key;
    while ((key = (TKey*)next())) {
        printf("%-40s %-15lld\n", key->GetName(), key->GetObjlen());
    }

    printf("--------------------------------------------------\n");
    file->Close();
    delete file;
}
