library(curl)

# page API de skitour
# https://skitour.fr/api/

# Notice package curl
# https://cran.r-project.org/web/packages/curl/vignettes/intro.html#Setting_handle_options

# Notice packahe httr
# https://cran.r-project.org/web/packages/httr/vignettes/quickstart.html 


# API pour les sorties
# cle: yYOlbMcps8rOuUCdVDTwbawLXO26IHLM

req <- curl_fetch_memory("https://skitour.fr/api/sorties?a=2022&s=84")
str(req)
parse_headers(req$headers)
jsonlite::prettify(rawToChar(req$content))
