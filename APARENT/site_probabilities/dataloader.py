from kipoi.data import SampleIterator
from kipoiseq.transforms import OneHot
from kipoiseq.extractors import SingleSeqVCFSeqExtractor, VariantSeqExtractor, MultiVariantsMatcher, \
    FastaStringExtractor
from kipoi.metadata import GenomicRanges
import pyranges


class FixedSeqPolyaDl(SampleIterator):
    seq_length = 205

    def __init__(
            self,
            gtf_file,
            fasta_file,
            # disable_infer_transcripts=True,
            # disable_infer_genes=True
    ):
        self.gtf_file = gtf_file
        self.fasta_file = fasta_file
        self.fasta_extractor = FastaStringExtractor(fasta_file, use_strand=True)

        rngs = pyranges.read_gtf(self.gtf_file)
        rngs = rngs.subset(lambda df: df.Feature == "transcript")

        rngs.Start -= 1000

        self.intervals = rngs
        self.one_hot = OneHot()
        self.n_upstream = self.seq_length // 2 + self.seq_length % 2
        self.n_downstream = self.seq_length // 2

    def __next__(self):
        sequence = self.one_hot(self.fasta_extractor.extract(interval))
        return {
            "inputs": {
                "sequence": sequence,
            },
            "metadata": {
                "ranges": GenomicRanges.from_interval(interval),
                "gene_id": interval.attrs.get('gene_id', ''),
                "transcript_id": interval.attrs.get('transcript_id', ''),
                "gene_biotype": interval.attrs.get('gene_biotype', '')
            }
        }

    def __iter__(self):
        return self
