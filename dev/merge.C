#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <Riostream.h>
#include <TFileMerger.h>

void merge(const char* files, const char* out, const char* skip)
{
  TFileMerger *merger = new TFileMerger(kFALSE); 
  merger->OutputFile(out);
  merger->SetFastMethod(kTRUE);

  TString strFiles(files);
  TObjArray *list = strFiles.Tokenize(" ");
  for (Int_t i = 0; i < list->GetEntriesFast(); i++) {
    TObjString *objStr = dynamic_cast<TObjString*>(list->At(i));
    if (!objStr) continue;
    merger->AddFile(objStr->GetString());
  }

  merger->AddObjectNames(skip);

  if (merger->PartialMerge(TFileMerger::kAll | TFileMerger::kIncremental | TFileMerger::kSkipListed)) {
    Printf("Merging success.");
  }
  else {
    Printf("Merging failed!");
  }

  delete list;
  list = 0;

  delete merger;
  merger = 0;
}
