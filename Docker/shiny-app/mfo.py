import os
from playwright.sync_api import sync_playwright


def render_html_to_image() -> str:
    """
    Render the HTML file to an image using Playwright

    Returns:
        str: The path to the downloaded image
    """

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()

        try:
            path = os.path.abspath("output.html")
            page.goto(f"file://{path}")
            page.wait_for_selector("svg")

            with page.expect_download() as download_info:
                page.click("text=Download SVG")
            download = download_info.value
            download.save_as(download.suggested_filename)

        finally:
            browser.close()
            return download.suggested_filename
