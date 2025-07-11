from typer import Typer, Option, Argument
from src.fetcher import fetch_papers
app = Typer(help="Fetch PubMed papers based on a query and export to CSV or console.")

@app.command()
def main(
    query: str = Argument(..., help="Search query to use with PubMed"),
    file: str = Option(None, "-f", "--file", help="Filename to save CSV results. If not provided, prints to console."),
    debug: bool = Option(False, "-d", "--debug", help="Print debug information during execution."),
):
    """
    Fetch research papers from PubMed matching the query.
    Filters for non-academic authors and outputs results.
    """
    fetch_papers(query=query, filename=file, debug=debug)
