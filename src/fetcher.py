
from Bio import Entrez, Medline
import pandas as pd
from typing import Optional
from src.utils import is_non_academic_author, extract_email
Entrez.email = "your-email@example.com"

def fetch_papers(query: str, filename: Optional[str], debug: bool):
    if debug:
        print(f"Searching PubMed for query: '{query}'")

    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=10)
        record = Entrez.read(handle)
        ids = record.get("IdList", [])
        handle.close()
    except Exception as e:
        print(" Error during PubMed search:", e)
        return

    if not ids:
        print("⚠️ No papers found for this query.")
        return

    if debug:
        print(f" Found {len(ids)} papers: {ids}")

    try:
        handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
        records = list(Medline.parse(handle))
        handle.close()
    except Exception as e:
        print(" Error fetching paper details:", e)
        return

    results = []
    for rec in records:
        pubmed_id = rec.get("PMID", "")
        title = rec.get("TI", "").replace("\n", " ")
        pub_date = rec.get("DP", "")
        authors = rec.get("FAU", [])
        affiliations = rec.get("AD", [])
        email = extract_email(affiliations)
        aff_list = affiliations if isinstance(affiliations, list) else [affiliations] * len(authors)

        non_acad_authors = []
        company_affiliations = []

        for author, aff in zip(authors, aff_list):
            if is_non_academic_author(aff):
                non_acad_authors.append(author)
                company_affiliations.append(aff)

        if non_acad_authors:
            results.append({
                "PubmedID": pubmed_id,
                "Title": title,
                "Publication Date": pub_date,
                "Non-academic Author(s)": "; ".join(non_acad_authors),
                "Company Affiliation(s)": "; ".join(company_affiliations),
                "Corresponding Author Email": email or "Not found"
            })

    if not results:
        print("⚠️ No papers found with non-academic authors.")
        return

    df = pd.DataFrame(results)

    if filename:
        df.to_csv(filename, index=False)
        print(f"✅ Results saved to {filename}")
    else:
        print(" Results:")
        print(df.to_string(index=False))
