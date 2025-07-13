import json
import logging
from typing import Any, Dict

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    repaq is a tool that generates compressed FASTQ files from FastQ files and
    provides summary statistics in JSON format.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="repaq",
            anchor="repaq",
            href="https://github.com/OpenGene/repaq",
            info="A tool for lossless compression of FASTQ files with quality control",
            extra="""
            repaq performs lossless compression of FASTQ files and provides quality control
            metrics including read counts, base counts, and compression results.
            """,
            doi="",
        )

        data_by_sample = dict()
        for f in self.find_log_files("repaq", filehandles=True):
            s_name, parsed_data = self.parse_repaq_json(f)
            if not s_name:
                continue
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            data_by_sample[s_name] = parsed_data

        # Ignore if module was called with no files
        if len(data_by_sample) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(data_by_sample)} repaq reports")

        # Write parsed data to file
        self.write_data_file(data_by_sample, "multiqc_repaq")

        # General stats table
        self.repaq_general_stats_table(data_by_sample)

        # Main detailed table
        self.add_section(
            name="Summary Statistics",
            anchor="repaq-stats",
            description="Statistics from repaq compression and QC analysis",
            plot=self.repaq_stats_table(data_by_sample),
        )

    def parse_repaq_json(self, f):
        """Parse a repaq JSON output file"""
        try:
            data = json.load(f["f"])
        except (json.JSONDecodeError, KeyError, TypeError):
            log.warning(f"Could not parse repaq JSON: {f['fn']}")
            return None, None

        # Extract sample name from filename
        s_name = self.clean_s_name(f["fn"], f)

        # Validate required fields
        required_fields = ['result', 'fastq_reads', 'rfq_reads', 'fastq_bases', 'rfq_bases']
        if not all(field in data for field in required_fields):
            log.warning(f"Missing required fields in repaq JSON: {f['fn']}")
            return None, None

        return s_name, data

    def repaq_general_stats_table(self, data_by_sample):
        """Add key stats to the general statistics table"""
        headers = {
            'result': {
                'title': 'Result',
                'description': 'Compression result status',
                'scale': 'RdYlGn',
                'format': lambda x: 'PASS' if x == 'passed' else 'FAIL'
            },
            'fastq_reads': {
                'title': 'Reads',
                'description': 'Number of reads in input FASTQ',
                'scale': 'Blues',
                'format': '{:,.0f}',
                'shared_key': 'read_count'
            },
            'fastq_bases': {
                'title': 'Bases',
                'description': 'Number of bases in input FASTQ',
                'scale': 'Greens',
                'format': '{:,.0f}',
                'shared_key': 'base_count'
            }
        }

        self.general_stats_addcols(data_by_sample, headers)

    def repaq_stats_table(self, data_by_sample):
        """Create a detailed statistics table"""
        headers = {
            'result': {
                'title': 'Result',
                'description': 'Compression result status',
                'scale': 'RdYlGn'
            },
            'msg': {
                'title': 'Message',
                'description': 'Status message'
            },
            'fastq_reads': {
                'title': 'Input Reads',
                'description': 'Number of reads in input FASTQ file',
                'scale': 'Blues',
                'format': '{:,.0f}'
            },
            'rfq_reads': {
                'title': 'Output Reads',
                'description': 'Number of reads in output RFQ file',
                'scale': 'Blues',
                'format': '{:,.0f}'
            },
            'fastq_bases': {
                'title': 'Input Bases',
                'description': 'Number of bases in input FASTQ file',
                'scale': 'Greens',
                'format': '{:,.0f}'
            },
            'rfq_bases': {
                'title': 'Output Bases',
                'description': 'Number of bases in output RFQ file',
                'scale': 'Greens',
                'format': '{:,.0f}'
            }
        }

        table_config = {
            'namespace': 'repaq',
            'id': 'repaq_stats_table',
            'title': 'repaq: Summary Statistics',
            'col1_header': 'Sample Name'
        }

        return table.plot(data_by_sample, headers, table_config)

