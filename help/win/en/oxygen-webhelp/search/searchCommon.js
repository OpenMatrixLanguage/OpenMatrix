/**
 * An object containing the search result for a single topic/HTML page.
 * Contains pointer to the topicID, title, short description and the list of words that were found.
 *
 * @param {string} topicID The ID of the topic. Can be used to identify unique a document in the search result.
 * @param {string} relativePath The relative path to the topic.
 * @param {string} title The topic title.
 * @param {string} shortDescription The topic short description.
 * @param {[string]} words The array with words contained by this topic.
 * @param {int} scoring The search scoring computed for this document.
 * @param {int} startsWith The number used to display 5 stars ranking.
 * @param {int} resultID The search result ID.
 * @param {int} linkID The search link ID.
 * @constructor
 */
function SearchResultInfo(topicID, relativePath, title, shortDescription, words, scoring, starWidth, resultID, linkID) {
    this.topicID = topicID;
    this.relativePath = relativePath;
    this.title = title;
    this.shortDescription = shortDescription;
    this.words = words;
    this.scoring = scoring;
    this.starWidth = starWidth;
    this.resultID = resultID;
    this.linkID = linkID;
    this.similarResults = [];
}

/**
 * Last displayed search results items.
 * @type {Array}
 */
var lastSearchResultItems = [];

/**
 * Last displayed search result.
 * @type {Array}
 */
var lastSearchResult;

/**
 * Pre process search result to compute similar results and scoring.
 * The lastSearchResultItems variable will be updated.
 *
 * @param searchResult The seach result to process.
 * @param whDistribution The WebHelp distribution.
 */
function preprocessSearchResult(searchResult, whDistribution) {
    lastSearchResult = searchResult;
    lastSearchResultItems = [];

    var wh_mobile =
        (typeof whDistribution != 'undefined') && whDistribution == 'wh-mobile';
    var wh_Classic =
        (typeof whDistribution != 'undefined') && whDistribution == 'wh-classic';


    if (searchResult.documents !== undefined && searchResult.documents.length > 0) {
        var allPages = searchResult.documents;

        // The score for fist item
        var ttScore_first = 1;
        if (allPages.length > 0) {
            ttScore_first = allPages[0].scoring;
        }

        var currentSimilarPage={};
        for (var page = 0; page < allPages.length; page++) {
            /*debug("Page number: " + page);*/

            if (allPages[page].relativePath == 'toc.html') {
                continue;
            }

            var starWidth = 0;
            if (typeof webhelpSearchRanking != "undefined" && webhelpSearchRanking) {
                var hundredPercent = allPages[page].scoring + 100 * allPages[page].words.length;
                var numberOfWords = allPages[page].words.length;
                /*debug("hundredPercent: " + hundredPercent + "; ttScore_first: " + ttScore_first + "; numberOfWords: " + numberOfWords);*/
                var ttScore = allPages[page].scoring;

                // Fake value
                var maxNumberOfWords = allPages[page].words.length;
                starWidth = (ttScore * 100 / hundredPercent) / (ttScore_first / hundredPercent) * (numberOfWords / maxNumberOfWords);
                starWidth = starWidth < 10 ? (starWidth + 5) : starWidth;
                // Keep the 5 stars format
                if (starWidth > 85) {
                    starWidth = 85;
                }
            }

            var idLink = 'foundLink' + page;
            var idResult = 'foundResult' + page;

            // topicID, relativePath, title, shortDescription, words, scoring, starWidth, resultID, linkID, similarResults
            var csri = new SearchResultInfo(
                allPages[page].topicID,
                allPages[page].relativePath,
                allPages[page].title,
                allPages[page].shortDescription,
                allPages[page].words,
                allPages[page].scoring,
                starWidth,
                idResult,
                idLink
            );

            // Similar pages
            var similarPages = !wh_mobile && similarPage(allPages[page], allPages[page - 1]);
            if (!similarPages) {
                currentSimilarPage = csri;
                lastSearchResultItems.push(csri);
            } else {
                currentSimilarPage.similarResults.push(csri);
            }

        }
    }
}

/**
 * Compute the HTML to be displayed in the search results page.
 *
 * @param whDistribution The string with WebHelp distribution. One of wh-classic, wh-mobile or wh-responsive.
 * @param pageNumber The page number to display.
 * @param totalPageNumber The total page number.
 * @param itemsPerPage The number of items to display on a page.
 * @returns {string} The HTML to be displayed as search result.
 */
function computeHTMLResult(whDistribution, pageNumber, totalPageNumber, itemsPerPage) {
    // Empty jQuery element
    var results = $();

    var $wh_search_results_items = $();

    if (lastSearchResult.searchExpression.length > 0) {
        if (lastSearchResultItems.length > 0) {
            $wh_search_results_items = $('<div/>', {
                class: 'wh_search_results_items'
            });

            // Start and end index depending on the current presented page
            var s = 0;
            var e = lastSearchResultItems.length;

            if (typeof pageNumber != "undefined" && typeof itemsPerPage != "undefined") {
                s = (pageNumber - 1) * itemsPerPage;
                e = Math.min(s + itemsPerPage, lastSearchResultItems.length);
            }

            // Result for: word1 word2
            var txt_results_for = "Results for:";
            var $headerHTML = $('<div/>', {
                class: 'wh_search_results_header'
            });

            var $whSearchResultsHeaderDocs = $('<div/>', {
                class: 'wh_search_results_header_docs'
            }).text(
                lastSearchResultItems.length +
                ' ' +
                getLocalization(txt_results_for) + ' '
            );

            var $span = $('<span/>', {
                class: 'wh_search_expression'
            }).text(lastSearchResult.originalSearchExpression);

            $whSearchResultsHeaderDocs.append($span);
            $headerHTML.append($whSearchResultsHeaderDocs);

            if (typeof pageNumber != "undefined" && typeof totalPageNumber != "undefined" && totalPageNumber > 1) {
                var $wh_search_results_header_pages = $('<div/>', {
                    class: 'wh_search_results_header_pages'
                }).text(getLocalization('Page') + ' ' + pageNumber + '/' + totalPageNumber);
                $headerHTML.append($wh_search_results_header_pages);
            }

            $wh_search_results_items.append($headerHTML);

			// EXM-38967 Start numbering
            var start = (pageNumber - 1) * 10 + 1;
            var $ol = $('<ol/>', {
                class: 'searchresult',
                start: start
            });

            for (var page = s; page < e; page++) {
                var csri = lastSearchResultItems[page];

                var hasSimilarPages =
                    csri.similarResults != null &&
                    csri.similarResults.length > 0;

                var siHTML = computeSearchItemHTML(
                    csri,
                    whDistribution,
                    hasSimilarPages,
                    null);
                $ol.append(siHTML);

                if (hasSimilarPages) {
                    // Add HTML for similar pages
                    for (var smPage = 0; smPage < csri.similarResults.length; smPage++) {
                        var simHTML = computeSearchItemHTML(
                            csri.similarResults[smPage],
                            whDistribution,
                            false,
                            csri.resultID);

                        $ol.append(simHTML);
                    }
                }
            }

            $wh_search_results_items.append($ol);

            if ($wh_search_results_items.find('li').length == 0) {
                $wh_search_results_items = $('<div/>', {
                    class: 'wh_search_results_for'
                });
                var $span = $('<span/>', {
                    class: 'wh_search_expression'
                }).text(lastSearchResult.originalSearchExpression);

                $wh_search_results_items.append($span);
                }
            } else {
            $wh_search_results_items = $('<div/>', {
                class: 'wh_search_results_for'
            }).text(getLocalization('Search no results') + ' ');
            var $span = $('<span/>', {
                class: 'wh_search_expression'
            }).text(lastSearchResult.originalSearchExpression);
            $wh_search_results_items.append($span);
            }
        } else {
        // Search expression is empty. If there are stop words, display a message accordingly
        if (lastSearchResult.excluded.length > 0) {
            $wh_search_results_items = $();
            var $p = $('<p/>', {
                class: 'wh_search_results_for'
            }).text(getLocalization("no_results_only_stop_words1"));
            $wh_search_results_items.append($p);

            $p.text(getLocalization('no_results_only_stop_words2'));
            $wh_search_results_items.append($p);
        }
    }

    return $wh_search_results_items;
}

function computeSearchItemHTML(searchItem, whDistribution, hasSimilarPages, similarPageID) {
    // New empty jQuery element
    var htmlResult = $();

    var wh_mobile =
        (typeof whDistribution != 'undefined') && whDistribution == 'wh-mobile';
    var wh_Classic =
        (typeof whDistribution != 'undefined') && whDistribution == 'wh-classic';

    var allSearchWords = lastSearchResult.searchExpression.split(" ");

    var tempPath = searchItem.relativePath;

    // EXM-27709 START
    // Display words between '<' and '>' in title of search results.
    var tempTitle = searchItem.title;
    // EXM-27709 END
    var tempShortDesc = searchItem.shortDescription;
    var starWidth = searchItem.starWidth;
    var rankingHTML = $();

    if (!wh_mobile && (typeof webhelpSearchRanking != 'undefined') && webhelpSearchRanking) {
        // Add rating values for scoring at the list of matches
        rankingHTML =  $("<div/>", {
            id: 'rightDiv'
        });
        if (displayScore) {
            rankingHTML.attr('title', 'Score: ' + searchItem.scoring);
        }

        var rankingStar =
            $('<div/>', {
                id: 'star'
            }).append(
                $('<div/>', {
                    id: 'star0',
                    class: 'star'
                }).append(
                    $('<div/>', {
                        id: 'starCur0',
                        class: 'curr',
                        style: 'width: ' + starWidth + 'px'
                    }).append(
                        $('<br/>', {
                            style: 'clear: both;'
                        })
                    )
                )
            );
        rankingHTML.append(rankingStar);
    }

    var finalArray = searchItem.words;
    var arrayStringAux = [];
    var arrayString = '';

    for (var x in finalArray) {
        if (finalArray[x].length >= 2 || useCJKTokenizing || (indexerLanguage == "ja" && finalArray[x].length >= 1)) {
            arrayStringAux.push(finalArray[x]);
        }
    }
    arrayString = arrayStringAux.toString();

    // Add highlight param
    if (!wh_Classic && !wh_mobile) {
        tempPath += '?hl=' + encodeURIComponent(arrayString);
    }

    var idLink = searchItem.linkID;
    var idResult = searchItem.resultID;

    var link = '';

    if (wh_Classic) {
        link = 'return openAndHighlight(\'' + tempPath + '\', Array(';
        for (var i in arrayStringAux) {
            link +='\'' + arrayStringAux[i] + '\', ';
        }
        link = link.substr(0, link.length-2) + '))';
    }

    // Similar pages
    if (similarPageID == null) {
        htmlResult = $('<li/>', {
            id: idResult
        });

        var $a = $('<a/>', {
            id: idLink,
            href: tempPath,
            class: 'foundResult'
        }).html(tempTitle);
        if (wh_Classic) {
            $a.attr('onclick', link);
        }

        htmlResult.append($a);
    } else {
        htmlResult = $('<li/>', {
            id: idResult,
            class: 'similarResult',
            'data-similarTo': similarPageID
        });

        var $a = $('<a/>', {
            id: idLink,
            href: tempPath,
            class: 'foundResult'
        }).html(tempTitle);
        if (wh_Classic) {
            $a.attr('onclick', link);
        }

        htmlResult.append($a);
    }

    // Also check if we have a valid description
    if ((tempShortDesc != "null" && tempShortDesc != '...')) {
        var $shortDescriptionDIV = $('<div/>', {
            class: 'shortdesclink'
        }).html(tempShortDesc);

        // Highlight the search words in short description
        for (var si = 0; si < allSearchWords.length; si++) {
            var sw = allSearchWords[si];
            $shortDescriptionDIV.highlight(sw, 'search-shortdescription-highlight');
        }

        htmlResult.append($shortDescriptionDIV);
    }

    // Empty jQuery element
    var searchItemInfo = $('<div/>', {
        class: 'missingAndSimilar'
    });

    // Relative Path
    $a = $('<a/>', {
        href: tempPath
    }).html(searchItem.relativePath);
    if (wh_Classic) {
        $a.attr('onclick', link);
    }

    var relPathStr = $('<div/>', {
        class: 'relativePath'
    }).append($a);

    searchItemInfo.append(relPathStr);

    // Missing words
    if (!wh_mobile && allSearchWords.length != searchItem.words.length) {
        var missingWords = [];
        allSearchWords.forEach(function (word) {
            if (searchItem.words.indexOf(word) == -1) {
                missingWords.push(word);
            }
        });

        var missingHTML = $('<div/>', {
            class: 'wh_missing_words'
        });
        missingHTML.html(getLocalization('missing') + ' : ');


        for (var widx = 0; widx < missingWords.length; widx++) {
            var $span = $('<span/>', {
                class: 'wh_missing_word'
            }).html(missingWords[widx]);
            missingHTML.append($span).append(' ');
        }

        searchItemInfo.append(missingHTML);
    }

    if (!wh_mobile && hasSimilarPages) {
        var $similarHTML = $('<a/>', {
            class: 'showSimilarPages',
            onclick: 'showSimilarResults(this)'
        }).html(getLocalization('Similar results') +  '...');

        searchItemInfo.append($similarHTML);
    }

    if (rankingHTML.html() != '' && searchItemInfo.html() != '') {
        var $searchItemAdditionalData = $('<div/>', {
            class: 'searchItemAdditionalData'
        }).append(searchItemInfo).append(rankingHTML);

        htmlResult.append($searchItemAdditionalData);
    } else if (searchItemInfo.html() != '') {
        htmlResult.append(searchItemInfo);
    } else if (rankingHTML.html() != '') {
        htmlResult.append(rankingHTML);
    }

    return htmlResult;
}

/**
 * @description Compare two result pages to see if there are similar
 * @param result1 Result page
 * @param result2 Result page
 * @returns {boolean} true - result pages are similar
 *                    false - result pages are not similar
 */
function similarPage(result1, result2) {
    var toReturn = false;

    if (result1 === undefined || result2 === undefined) {
        return toReturn;
    }

    var pageTitle1 = result1.title;
    var pageShortDesc1 = result1.shortDescription;

    var pageTitle2 = result2.title;
    var pageShortDesc2 = result2.shortDescription;

    if (pageTitle1.trim() == pageTitle2.trim() && pageShortDesc1.trim() == pageShortDesc2.trim()) {
        toReturn = true;
    }

    return toReturn;
}

/**
 * When it is true, then the score is displayed as tooltip.
 *
 * @type {boolean}
 */
var displayScore = false;

/**
 * @description Show similar results that are hidden by default
 * @param link Link clicked to show similar results
 */
function showSimilarResults(link) {
    var parentLiElement = $(link).parents('li[id]');
    var currentResultId = parentLiElement.attr('id');

    $('[data-similarTo="' + currentResultId + '"]').toggle();
    $(link).toggleClass('expanded');
}