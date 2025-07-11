
# PubMed Paper Fetcher CLI

This is a command-line Python tool that helps you find biomedical research papers from PubMed. It looks for papers that have at least one author working at a pharmaceutical or biotech company, and lets you save those results as a CSV file or print them to the terminal.

---

## How the Code is Organized

```
pubmed-paper-fetcher/
├── src/
│   ├── cli.py         # Main CLI interface using Typer
│   ├── fetcher.py     # Code to search, fetch, and process PubMed data
│   ├── utils.py       # Helper functions (e.g., email extraction, heuristics)
├── tests/
│   └── test_utils.py  # Unit tests
├── pyproject.toml     # Poetry config with dependencies and script entry point
├── README.md          # Project documentation
├── results.csv        # Example CSV output
```

* `cli.py`: Starts the command-line tool using [Typer](https://typer.tiangolo.com)
* `fetcher.py`: Talks to the PubMed API and filters the results
* `utils.py`: Has helper logic for checking company affiliations and getting emails
* `test_utils.py`: Contains test cases for utility functions

---

## Features

* Search PubMed with flexible queries
* Identify papers with non-academic (company) authors
* Extract email addresses of corresponding authors when available
* Output the result as a CSV file or to the console
* Built with modern tools: **Python 3.10+**, **Typer**, **Poetry**, **Biopython**, **Pandas**
* Includes test coverage using `pytest`

---

## How I Detect Non-Academic Authors

I use a simple rule-based method:

* If the affiliation contains company-related words like `Inc`, `LLC`, `Pharma`, `Biotech`, etc., it is considered **non-academic**.
* If it also includes words like `University`, `Hospital`, `Institute`, etc., we assume it's academic and skip it.
* Emails are pulled using a regular expression from the affiliation text.

---

## Installation

Make sure Python 3.10 or higher is installed.

Clone the project and install dependencies using Poetry:

```bash
git clone https://github.com/your-username/pubmed-paper-fetcher.git
cd pubmed-paper-fetcher
poetry install
```

---

## How to Run the Program

Search and print to the terminal:

```bash
poetry run get-papers-list "covid vaccine"
```

Save results to a CSV file:

```bash
poetry run get-papers-list "cancer immunotherapy" -f results.csv
```

Enable debug mode to see more details:

```bash
poetry run get-papers-list "covid" --debug
```

---

## Running the Tests

First install testing dependencies:

```bash
poetry add --group dev pytest
```

Then run:

```bash
poetry run pytest
```

Tests cover the academic detection and email extraction functions.

---

## Tools and Libraries Used

| Tool                                           | Description                                                         |
| ---------------------------------------------- | ------------------------------------------------------------------- |
| [Python 3.10+](https://www.python.org/)        | Programming language                                                |
| [Typer](https://typer.tiangolo.com/)           | CLI framework                                                       |
| [Biopython](https://biopython.org/)            | Interface to the PubMed API                                         |
| [Pandas](https://pandas.pydata.org/)           | CSV and data handling                                               |
| [Poetry](https://python-poetry.org/)           | Dependency and script manager                                       |
| [pytest](https://docs.pytest.org/)             | Testing framework                                                   |
| [ChatGPT (OpenAI)](https://openai.com/chatgpt) | Used to assist in debugging code, allowed by assignment |

---

## Assignment Checklist


| Requirement                           | Status    |
| ------------------------------------- | --------- |
| Matches task description              | Success   |
| Outputs all required CSV fields       | Success   |
| CLI with `--debug` and `--file` flags | Success   |
| Uses type hints                       | Success   |
| Structured code                       | Success   |
| Handles API errors and edge cases     | Success   |
| Has unit tests                        | Success   |
| Clearly documents LLM use             | Success   |

---

## Author

**Vidya Vandana**

---

## License
MIT License

Copyright (c) 2025 Vidya Vandana

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the “Software”), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

=======
