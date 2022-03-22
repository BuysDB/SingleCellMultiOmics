import json

def mark_autodetect(libname):
    lib = libname.lower()
    if 'k9me3' in lib:
        return 'H3K9me3'
    if 'k4me1' in lib:
        return 'H3K4me1'
    if 'k4me3' in lib:
        return 'H3K4me3'
    if 'k27me3' in lib:
        return 'H3K27me3'
    if 'k36me3' in lib:
        return 'H3K36me3'

    return 'NA'


class SampleSheet:
    """
    Class for reading and writing sample sheet information
    """
    def __init__(self, path=None):
        self.data = {}
        if path is not None:
            self.load_json(path)

    def load_json(self,sample_sheet_location,merge=False):
        """
        Load the data from the specified json file.
        When merge=True the data will be added to the current data (used for merging sheets)
        Duplicates will be overwritten by the last load file
        """
        with open(sample_sheet_location) as f:
            if merge:
                for k,v in json.load(f).items():
                    if k not in self.data:
                        self.data[k] = v
                    else:
                        self.data[k].update(v)
            else:
                self.data = json.load(f)
    def get(self,key,default):
        return self.data.get(key,default)

    def __getitem__(self,k):
        return self.data[k]

    def layout_well2index(self, layout_name):
        """ Get the well to index mapping for the supplied layout """
        return {k:tuple(v) for k,v in self['well2coord'][self['layout_format'][layout_name]].items()}

    def index_to_condition(self, layout_name):
        well_to_condition = self['layouts'][layout_name]
        return {self['well2index'][self['layout_format'][layout_name]][k]:v for k,v in well_to_condition.items()}


    def index2well(self, layout_name: str) -> dict:
        """
        Returns a dictionary index->well
        {1: 'A1',
         2: 'A2',
         3: 'A3',
         4: 'A4',
         5: 'A5',
        """
        return {v:k for k,v in self['well2index'][self['layout_format'][layout_name]].items() }

    @property
    def libraries(self):

        libraries = set()
        for k in ['libtype','marks']:
            libraries.update( set(self.get(k,dict()).keys()))
        return libraries


    def drop_library(self, library):
        if type(library) is list:
            for l in library:
                self.drop_library(l)
            return
        for main_key in ['marks','condition','libtype','library_layout']:
            if main_key not in self.data:
                continue
            if library in self[main_key]:
                del self[main_key][library]

    def drop_mark(self, mark):
        """ Remove all libraries with the given mark """
        if type(mark) is list:
            for m in mark:
                self.drop_mark(m)
            return
        to_prune = []
        for sample, _mark in self['marks'].items():
            if mark==_mark:
                to_prune.append(sample)

        self.drop_library(to_prune)
