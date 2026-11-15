from geniml.bbclient import BBClient
from txdb import TxDbRegistry

registry = TxDbRegistry()

# List available databases
print(registry.list_txdb())

# Load a specific database (downloads and caches it automatically)
# Example: hg38 knownGene
txdb = registry.load_db("TxDb.Hsapiens.UCSC.hg38.knownGene")

# Access features
transcripts = txdb.transcripts()
exons = txdb.exons()

tss = transcripts.resize(width=1, fix="start")
promoters = tss.promoters(upstream=2000, downstream=200)

# Load ChIP-seq BED file from BEDbase as GenomicRanges
bbclient = BBClient(cache_folder="cache", bedbase_api="https://api.bedbase.org")
bedfile_id = "be4054acf6e3feeb4dc490e6430e358e"
bedfile = bbclient.load_bed(bedfile_id)
chipseq_dataset = bedfile.to_granges()

peaks_by_tss = chipseq_dataset.subset_by_overlaps(tss)
peaks_by_promoters = chipseq_dataset.subset_by_overlaps(promoters)

