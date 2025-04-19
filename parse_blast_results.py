from Bio.Blast import NCBIXML
import glob
import sys

# 从命令行参数读取模式（blastn 或 blastx）
if len(sys.argv) != 2 or sys.argv[1] not in ["blastn", "blastx"]:
    print("Usage: python parse_blast_results.py [blastn|blastx]")
    sys.exit(1)

mode = sys.argv[1]
xml_files = glob.glob(f"{mode}_*.xml")

print(f"Found {len(xml_files)} {mode.upper()} result files.\n")

for xml_file in xml_files:
    print(f"📄 Parsing {xml_file}")
    with open(xml_file) as handle:
        blast_record = NCBIXML.read(handle)

        top_hits = blast_record.alignments[:3]
        print(f"Top {len(top_hits)} hits:\n")

        for alignment in top_hits:
            hsp = alignment.hsps[0]
            print(f"  ➤ Title: {alignment.title[:80]}...")
            print(f"     Length: {alignment.length}")
            print(f"     E-value: {hsp.expect}")
            print(f"     Score: {hsp.score}")
            print(f"     Identity: {hsp.identities}/{hsp.align_length} ({100*hsp.identities/hsp.align_length:.1f}%)\n")

    print("-" * 60)
