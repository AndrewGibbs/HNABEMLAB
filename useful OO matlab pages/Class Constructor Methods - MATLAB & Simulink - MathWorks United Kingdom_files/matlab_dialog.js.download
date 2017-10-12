$(document).ready(function () {
    registerMatlabCommandDialogAction();
    registerOpenExampleAction();
});

function registerOpenExampleAction() {
    var openExampleButtons = $('.examples_short_list a.btn[href^="matlab:"]');
    $.each(openExampleButtons, function() {
        var href = $(this).attr('href');
        addOpenExampleLinkClickHandler($(this));
        var matlabCommand = getMatlabCommand(href);
        var openWithCommand = getOpenWithCommand(matlabCommand);
        if (window.ow !== undefined &&  openWithCommand) {
            ow.doesExampleExist(openWithCommand, function (status) {
                var matlabLink = $('a[href="matlab:openExample(\'' + openWithCommand+ '\')"]');
                if (status === 'true') {
                    matlabLink.html("Try it in MATLAB");
                    matlabLink.addClass('visible-xs').css('display', 'inline-block');
                    var dropDown = $('<div class="pull-right hidden-xs"><div class="btn-group"><button class="btn btn_color_blue dropdown-toggle dropdown-toggle" data-toggle="dropdown" href="#">Try this Example <span class="caret"></span></button><ul class="dropdown-menu dropdown-menu-right"><li><a href="#" class="analyticsOpenWith">Try it in your browser</a></li><li><a href="#">Try it in MATLAB</a></li></ul></div></div>');
                    $(matlabLink).parent().append(dropDown);
                    $(dropDown).find("li:first a").on('click', function(e) {
                        e.preventDefault();
                        ow.startOpenWith(openWithCommand);
                        ow.loadExample(openWithCommand);
                    });
                    $(dropDown).find("li:last a").on('click', function(e) {
                        e.preventDefault();
                        showMatlabDialog(matlabCommand);
                    });
                } else {
                    matlabLink.html("Try it in MATLAB");
                    matlabLink.css('display', 'inline-block');
                }
            });
        } else {
            if (window.ow !== undefined) {
                $(this).html("Try it in MATLAB");
                $(this).css('display', 'inline-block');
            }
        }
    });
}

function addOpenExampleLinkClickHandler(link) {
    $(link).on('click', function(e) {
        e.preventDefault();
        var href = $(this).attr('href');
        var matlabCommand = getMatlabCommand(href);
        showMatlabDialog(matlabCommand);
    });

}

$(window).bind('examples_cards_added', function(e) {
    $('.card_container a[href^="matlab:"]').hide();
});

function registerMatlabCommandDialogAction() {
    $('a[href^="matlab:"]').not('.card_container a[href^="matlab:"], .examples_short_list a.btn[href^="matlab:"]').on('click', function (e) {
        e.preventDefault();
        var href = $(this).attr('href');
        var matlabCommand = getMatlabCommand(href);
        showMatlabDialog(matlabCommand);
    });
}

function getMatlabCommand(href) {
    var matlabCommand = null;
    var match = href.match(/matlab:(.*)/);
    if (match) {
        matlabCommand = match[1];
    }
    return matlabCommand;
}

function getOpenWithCommand(matlabCommand) {
    var openWithCommand = null;
    var match = matlabCommand.match(/openExample\('(.*)'\)/);
    if (match) {
        openWithCommand = match[1];
    }
    return openWithCommand;
}

function showMatlabDialog(matlabCommand) {
    if (matlabCommand) {
        $("#matlab-command-dialog #dialog-body #dialog-matlab-command").text(matlabCommand);
    } else {
        $("#matlab-command-dialog #dialog-body #dialog-matlab-command").hide();
    }
    $("#matlab-command-dialog").modal();
}