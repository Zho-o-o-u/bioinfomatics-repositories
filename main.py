import json
import time
from collections import Counter, defaultdict
from datetime import date, datetime, timedelta

import fire
from requests_html import HTMLSession

MIN_SCORE = 4

OUT = "README.md"

url_one = "https://github.com/{user}?tab=repositories"
url_team = "https://github.com/{user}"

TITLE = "# Bioinfomatics Repositories\n\n"
DESCRIPTION = "Collection of bioinformatician and their repositories.\n\n"
UPDATE = f"**Update {date.today()}**\n\n"
ITEM1 = "- [Top50](https://github.com/xhqsm/bioinfomatics-repositories#Top50)\n"
ITEM2 = "- [Repositories](https://github.com/xhqsm/bioinfomatics-repositories#repositories)\n\n"
HEADER1 = "## Top50\n\n"
HEADER2 = "## Repositories\n\n"
TABLE_HEADER = "| ID | User | Type |Name | Stars |\n| :-: | :-: | :-: | :-: | :-: |\n"
template = """##### {user}/{repository}

ID = {id}, Stars = {star}, Language = {language}, Update = {update}

{description}
"""


def get_repositories(user: str, url: str = None) -> dict:
    rec = False
    cell = defaultdict(dict)
    if not ("?after=" in url and "?page=" in url):
        url = url.format(user=user)
    print(url)
    r = HTMLSession().get(url)
    h1_text = r.html.find('h1')[0].text
    if "504" in h1_text or "Whoa there!" in h1_text:
        time.sleep(7)
        yield from get_repositories(user, url)
    else:
        next_elements = None
        if url.endswith("repositories"):
            username = "One |" + r.html.find('h1>span')[0].text
            parents = r.html.find(
                'div[class="col-10 col-lg-9 d-inline-block"]')
            next_elements = r.html.find(
                'a[class="btn btn-outline BtnGroup-item"]')
        else:
            username = "Team |" + h1_text
            parents = r.html.find(
                'li[class="public source d-block py-4 border-bottom"]')
            next_elements = r.html.find('a[class="next_page"]')
        for e in parents:
            stars = e.find('a[href$="stargazers"]')
            if stars:
                star = int(stars[0].text.replace(",", ""))
                update = e.find('relative-time')[0].attrs["datetime"].split(
                    "T")[0]
                years = (datetime.now() -
                         datetime.strptime(update, "%Y-%m-%d")).days / 365
                # score filter
                score = star / (years + 1)
                if score >= MIN_SCORE:
                    language = e.find('span[itemprop="programmingLanguage"]')
                    language = language[0].text if language else None
                    if language:
                        r = e.find('a[itemprop="name codeRepository"]')[0]
                        repository = f"[{r.text}]({list(r.absolute_links)[0]})"
                        description = e.find('p[itemprop="description"]')
                        description = description[
                            0].text if description else "None"
                        user_link = f"[{user}](https://github.com/{user})"
                        cell[repository] = {
                            "star": star,
                            "username": username,
                            "user": user_link,
                            "description": description.strip(),
                            "language": language,
                            "update": update,
                        }
        if next_elements:
            for e in next_elements:
                if e.text == "Next":
                    rec = True
                    url = list(e.absolute_links)[0]
                    break
        yield cell
        if rec:
            yield from get_repositories(user, url)


def get_users(config: str) -> tuple:
    with open(config, encoding="utf-8") as f:
        users = json.load(f)
        return (set(x.lower() for x in users["ones"]),
                set(x.lower() for x in users["teams"]))


def get_content(config: str) -> dict:
    content = defaultdict(dict)
    ones, teams = get_users(config)
    for user in ones:
        for cell in get_repositories(user, url_one):
            content.update(cell)
    for user in teams:
        for cell in get_repositories(user, url_team):
            content.update(cell)
    return content


def get_exist() -> list:
    repositories = []
    with open(OUT, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#####"):
                row = line.strip().split()
                repositories.append(row[1])
    return repositories


def main(config: str = "default.json"):
    exist = get_exist()
    usernames = {}
    counter = Counter()
    content = get_content(config)
    with open(OUT, "w", encoding="utf-8") as out:
        for line in [
                TITLE, DESCRIPTION, ITEM1, ITEM2, UPDATE, HEADER1, TABLE_HEADER
        ]:
            out.write(line)
        for repository, cell in content.items():
            user = cell["user"]
            usernames[user] = cell["username"]
            counter.update({user: cell["star"]})
        for i, (user, stars) in enumerate(counter.most_common(), start=1):
            line = f"| {i} | {user} | {usernames[user]} | {stars} |\n"
            out.write(line)
            if i == 50:
                break
        out.write(HEADER2)
        for i, (r, c) in enumerate(sorted(content.items(),
                                          key=lambda t: t[1]["star"],
                                          reverse=True),
                                   start=1):
            if f'{c["user"]}/{r}' not in exist:
                r = r + " *new+*"
            c.update({"id": i, "repository": r})
            out.write(template.format_map(c))


if __name__ == "__main__":
    fire.Fire(main)
