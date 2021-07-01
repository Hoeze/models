from kipoi.data import SampleIterator
from kipoiseq.transforms import OneHot
from kipoiseq.extractors import VariantSeqExtractor, SingleVariantMatcher
from kipoi.metadata import GenomicRanges

import pyranges


class FixedSeqPolyaDl(SampleIterator):
    seq_length = 205

    def __init__(
            self,
            gtf_file,
            fasta_file,
            vcf_file,
            # disable_infer_transcripts=True,
            # disable_infer_genes=True
            # **kwargs,
    ):
        self.gtf_file = gtf_file
        self.fasta_file = fasta_file
        self.variant_seq = VariantSeqExtractor(fasta_file)
        self.fasta_extractor = self.variant_seq.fasta
        self.vcf_file = vcf_file

        self.n_upstream = self.seq_length // 2 + self.seq_length % 2
        self.n_downstream = self.seq_length // 2

        rngs = pyranges.read_gtf(self.gtf_file)
        rngs = rngs.subset(lambda df: df.Feature == "transcript")

        rngs.Start = rngs.End - self.n_upstream + 1
        rngs.End = rngs.End + self.n_downstream

        self.matcher = SingleVariantMatcher(
            vcf_file,
            pranges=rngs,
            interval_attrs=['gene_id', 'gene_name', 'transcript_id', 'gene_biotype']
        )

        self.pairs = iter(self.matcher)
        self.one_hot = OneHot()

    def __next__(self):
        interval, variant = next(self.pairs)
        return {
            "inputs": {
                "ref_seq": self.one_hot(self.fasta_extractor.extract(interval)),
                "alt_seq": self.one_hot(self.variant_seq.extract(interval, variant, anchor=self.seq_length // 2))
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
