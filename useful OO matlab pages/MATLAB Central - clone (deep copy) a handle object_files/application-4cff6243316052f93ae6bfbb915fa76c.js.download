function loadPopoverClickHandlers(selector) {

   $(selector).click(function(event){
    event.preventDefault();
    var link = event.currentTarget;
    var id = $(link).data("cp-link-id")
    var target = $("div[data-cp-popup-id=" + id + "]");

    // prevent clipping of profile link via ancestor panels' overflow value
    target.parents(".content_container .panel-body").css("overflow-y", "visible");

    var windowWidth = $(window).width();
    var windowHeight = $(window).height();
    var windowScrollHeight = $(window).scrollTop();

    var profileLinkLeft = $(link).offset().left;
    var profileLinkRight = $(link).offset().left + target.outerWidth();
    var profileLinkBottom = $(link).offset().top + target.outerHeight();

    var leftPosition = $(link).position().left - 20;
    var topPosition = $(link).position().top - 25;

    // adjust horizontal placement so profile link is clipped as little as possible on the right
    // this assumes that there is no horizontal scroll due to the resposive design
    if (profileLinkRight > windowWidth) {

        // find maximum left adjustment
        if ((target.outerWidth() + 6) > windowWidth) {
          leftPosition -= profileLinkLeft - 24;
        }
        else {
          leftPosition -= (profileLinkRight - windowWidth - 17);
        }
    }

    // adjust vertical placement so profile link isn't clipped on the bottom
    if (profileLinkBottom > (windowHeight + windowScrollHeight)) {
        topPosition -= ((profileLinkBottom - windowHeight - windowScrollHeight) - 21);
    }

    target.css({top: topPosition, left: leftPosition, position:'absolute'});
    target.show();
  });
}

function loadPopoverMouseleaveHandlers(selector) {

  $(selector).mouseleave(function(event){

    // set any ancestor panels' overflow back to original value
    $(this).parents(".content_container .panel-body").css("overflow-y", "auto");

    $(event.currentTarget).hide();
  });
}

(function($){

  $( document ).ready(function() {

    loadPopoverClickHandlers("a[data-cp-link-id]");
    loadPopoverMouseleaveHandlers("div[data-cp-popup-id]");

    $('body').on('touchstart', function(e){
      $.each($("div[data-cp-popup-id]:visible"), function(index, profile_div) {
        if (!($(e.target).data("cp-popup-id") == $(profile_div).data("cp-popup-id")) &&
            !($(e.target).parents("[data-cp-popup-id]").size() > 0)) {
          $(profile_div).hide();
        };
      });
    });

  });

})(jQuery);
// This is a manifest file that'll be compiled into application.js, which will include all the files
// listed below.
//
// Any JavaScript/Coffee file within this directory, lib/assets/javascripts, vendor/assets/javascripts,
// or vendor/assets/javascripts of plugins, if any, can be referenced here using a relative path.
//
// It's not advisable to add code directly here, but if you do, it'll appear at the bottom of the
// the compiled file.
//
// WARNING: THE FIRST BLANK LINE MARKS THE END OF WHAT'S TO BE PROCESSED, ANY BLANK LINE SHOULD
// GO AFTER THE REQUIRES BELOW.
//
// require jquery
// require jquery_ujs

;
