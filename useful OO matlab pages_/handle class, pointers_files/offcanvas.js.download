(function() {

$(document).ready(function () {
  $('[data-toggle=offcanvas]').click(function () {
    $('.row-offcanvas').toggleClass('active');
    $('#nav_toggle').toggleClass('active');
    $('#responsive_offcanvas').removeClass('no_animate');
      setCookie('MW_toc_visible', $('.row-offcanvas').hasClass('active'), getRootPath());
  });
});

if (!isMobileWidth() && isTocOpen()) {
    $('#responsive_offcanvas').addClass('no_animate');
    $('.row-offcanvas').addClass('active');
    $('#nav_toggle').addClass('active');
}

 function isTocOpen() {
        var tocCookie = getCookie('MW_toc_visible');
        return tocCookie === null || tocCookie === 'true';
    }

    function getCookie(name) {
        var cookies = document.cookie.split(';');
        for (var i = 0; i < cookies.length; i++) {
            var cookie = cookies[i].replace(/^\s+|\s+$/g, '');
            if (cookie.indexOf(name) === 0) {
                return cookie.substring(name.length + 1, cookie.length);
            }

        }
        return null;
    }

    function setCookie(name, value, path) {
        var date = new Date();
        date.setTime(date.getTime() + (7 * 24 * 60 * 60 * 1000));
        var expiresDate = date.toGMTString();
        document.cookie = name + "=" + value
            + "; expires=" + expiresDate
            + "; path=" + path ;

    }

    // Get the root path. This function returns the top level by default. 
    function getRootPath() {
        var pathname = $(location).attr('pathname');
        var pathArray = pathname.split("/");
        var rootPath = "/" + pathArray[1];
        return rootPath;
    }

    function isMobileWidth () {
        return $(window).width() < 768;
    }
})();

