from .statistic import StatisticHistogram
import singlecellmultiomics.pyutils as pyutils
import collections
import pandas as pd
import matplotlib.pyplot as plt


class MappingQualityHistogram(StatisticHistogram):
    def __init__(self,args):
        StatisticHistogram.__init__(self, args)
        self.histogram = collections.Counter()
    def processRead(self,read):
        self.histogram[read.mapping_quality]+=1
    def __repr__(self):
        return f'The average mapping quality is {pyutils.meanOfCounter(self.histogram)}, SD:{pyutils.varianceOfCounter(self.histogram)}'


    def get_df(self):
        return pd.DataFrame.from_dict({'mq':self.histogram})

    def to_csv(self, path):
        self.get_df().to_csv(path)

    def plot(self, target_path, title=None):
        df = self.get_df() #,'UnmappedReads']]

        df['mq'].plot.bar(figsize=(10,4))
        ax = plt.gca()
        ax.set_xlabel('Mapping quality')
        ax.set_ylabel('Frequency (reads)')


        if title is not None:
            plt.title(title)
        plt.tight_layout()
        plt.savefig(target_path)

        ax.set_yscale('log')
        plt.tight_layout()
        plt.savefig(target_path.replace('.png','.log.png'))

        plt.close()
