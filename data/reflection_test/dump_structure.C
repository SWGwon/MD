void dump_structure(const char* filename) {
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        printf("Error opening file\n");
        return;
    }
    f->ls();
    TIter next(f->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())) {
        TClass *cl = TClass::GetClass(key->GetClassName());
        if (cl->InheritsFrom("TTree")) {
            TTree *tree = (TTree*)key->ReadObj();
            tree->Print();
        }
    }
}
