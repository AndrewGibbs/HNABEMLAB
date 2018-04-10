(function ($) {
    function getIntrnllinks() {
        var locationPath = filterPath(location.pathname);
        if (locationPath.charAt(0) === "/") {
            locationPath = locationPath.substr(1);
        }
        var anchors = $('a:not([data-toggle]):not([data-slide])');
        var intrnllnks = [];
        $.each(anchors, function (index, anchor) {
            var anchorPath = anchor.pathname.substring(anchor.pathname.indexOf("#"));
            anchorPath = filterPath(anchorPath);
            if (anchorPath.charAt(0) === "/") {
                anchorPath = anchorPath.substr(1);
            }
            if (locationPath == anchorPath && anchor.hash !== '') {
                intrnllnks.push(anchor);
            }
        });
        return intrnllnks;
    }

    function filterPath(string) {
        return string.replace(/(index.html)$/, '')
    }

    function addSmoothScroll() {
        $.each(getIntrnllinks(), function(index, anchor) {
            var hash = this.hash;
            var target = getInternalLinkTarget(hash);
            if (target.length > 0) {
                $(this).on('click', function (evt) {
                    evt.preventDefault();
                    var nextSibling = getNextSiblingForAnchorTarget(target);
                    var scrollParameter = getScrollParameter();
                    var scrollTop = target.offset().top - getScrollTopAdjustment();
                    $(scrollParameter).animate({scrollTop: scrollTop}, 700, function () {
                        target.trigger('afterScroll');
                    });
                 
                    if (history.pushState) {
                        history.pushState(null, null, hash);
                    }
                    else {
                        location.hash = hash;
                    }

                })
            }
        });
    }

    function getScrollParameter() {
        //On IE and FF, the slow scroll parameter is the HTML dom element. On webkit, it is the body.
        return $('html, body');
    }


    function getInternalLinkTarget(hash) {
        //search for anchor with given hash as "name" atrribute value;
        var target = [];

        //Remove the first '#' character from the name attribute. Escape any special character from the name/id.
        var escapedHash = hash.substring(1).replace(/([;&,.+*~':"!^#$%@\[\]\(\)=>\|])/g, '\\$1');

        target = $("#" + escapedHash);
        if(target.length === 0) {
            target = $('[name=' + escapedHash + ']');
        }
        return target;
    }

    function getNextSiblingForAnchorTarget(target) {
        var nextSibling;
        nextSibling = target.next();
        return nextSibling;
    }

    function getScrollTopAdjustment() {
        var scrollTop = 0;
        scrollTop = $('.sticky_header_container').height() + 10;
        return scrollTop;
    }

    function scrollToHash() {
        if (location.hash.length > 0) {
            var target = getInternalLinkTarget(location.hash);
            if (target.length > 0) {
                var nextSibling = getNextSiblingForAnchorTarget(target);
                var scrollParameter = getScrollParameter();
                var scrollTop = target.offset().top - getScrollTopAdjustment();
                $(scrollParameter).scrollTop(scrollTop);
                target.trigger('afterScroll');
            }
        }
    }
    
    $(function() {
        addSmoothScroll();
    });

    $(window).on('load', function() {
        scrollToHash();
    });
    
    

})(jQuery);
