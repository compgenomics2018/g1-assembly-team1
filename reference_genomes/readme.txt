awk '/Klebsiella/ && $13 == "Complete" {print }' assembly_summary.txt
