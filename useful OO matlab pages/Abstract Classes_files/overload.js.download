$(document).ready( function() {
    /* Windows load */
    $(window).on("load", function() {
        if (typeof callOverloadCommand === 'function') {
            callOverloadCommand(); 
        }
    });
});

function callOverloadCommand() {
    if(window.location.href.indexOf("overload") > -1) {
        var overloadCommand = $.getParameterByName("overload");
        window.location = 'overload:' + overloadCommand;
    }
}

//===============================================================================
////==================overload page templates======================================

window.JST = window.JST || {};

// Overload function message bar
JST['messagebar'] = _.template(
    '<div id="message_bar" class="alert alert-info alert-dismissible fade in" role="alert" style="display:none;">' +
        '<div class="messagebar_info"></div>' +
        '<button type="button" class="close" data-dismiss="alert" aria-label="Close">' +
            '<span aria-hidden="true">×</span>' +
        '</button>' +
        '<a href="#" id="overload_message_bar_title" data-toggle="modal" data-target="#helpcenter-overload-function-dialog"><%= linktext %></a>' +
    '</div>'
);

// Modal Dialog with overload function list body
JST['overloaddialog'] = _.template(
    '<div class="modal fade" id="helpcenter-overload-function-dialog" tabindex="-1" role="dialog" aria-labelledby="helpcenterOverloadFunctionDialog">' +
        '<div class="modal-dialog modal-lg" role="document">' +
            '<div class="modal-content">' +
                '<div class="modal-header">' +
                    '<button type="button" class="close" data-dismiss="modal" aria-label="Close"><span aria-hidden="true">&times;</span></button>' +
                    '<h3 class="overload-dialog-title modal-title" id="helpcenterOverloadFunctionDialog" style="color:#C45400"><%= overloadjson.overloaddialogtitle %></h3>' +
                '</div>' +
                '<div class="modal-body">' +
                    '<pre id="dialog-overload">' +
                        '<% _.each(overloadjson.list, function(product) { %>' +
                            '<p style="font-weight:bold; font-size:11px;"><%= product.displayname %></p>' +
                            '<ul class="list-unstyled" style="margin-bottom:0px;">' +
                                '<% _.each(product.functionlist, function(overloadFunc, i) { %>' +
                                '<li class="overload_func" style="margin-left: 25px; font-size:12px;">' +
                                    '<a class="overload_func_link" href="<%= overloadjson.docroot%>/<%= overloadFunc.path %>?overload=<%=overloadjson.term%>%20<%=overloadjson.methodorproperty%>">' +
                                      '<%= overloadFunc.entity_name %>' +
                                    '</a>' +
                                    '<span> - <%= overloadFunc.summary %></span>' +
                                '</li>' +
                                '<% }) %>' +
                            '</ul>' +
                        '<% }) %>' +
                    '</pre>' +
                '</div>' +
                '<div class="modal-footer">' +
                    '<button type="button" class="btn btn-default" data-dismiss="modal">Close</button>' +
                '</div>' +
            '</div>' +
        '</div>' +
    '</div>'
);


function getOverloadFunctionsListHtml(overloadFunctionsJson) {
    var jsonData = {overloadjson: overloadFunctionsJson};
    return JST['overloaddialog'](jsonData);

}

function getOverloadMessageBarInfo(overloadFunctionsJson) {
    var jsonData = {linktext: overloadFunctionsJson.overloadmessagebarinfo};
    return JST['messagebar'](jsonData);
}

function storeMessagebarStatus() {
    sessionStorage.setItem('overloadmessagebar', 'hide');
}

function loadMessagebarStatus() {
    return sessionStorage.getItem('overloadmessagebar');
}

function populateOverloadStrings(overloadjson) {
    var messagebarVisible = loadMessagebarStatus();
    if (overloadjson !== 'undefined' && overloadjson !== null && messagebarVisible !== 'hide') {        
        var messagebar = getOverloadMessageBarInfo(overloadjson);
        if ($('#message_bar').length) {
                
        } else {
            $('.sticky_header_container').addClass("messagebar_active"); // fix the left nav breadcrumb scroll fixed top
            $('.sticky_header_container').prepend(messagebar);
            $("#message_bar").show();   
            $(window).trigger('content_resize'); // fix the scroll
            $('div#message_bar').on('close.bs.alert', function () {
                //
                storeMessagebarStatus();
                $('.sticky_header_container').removeClass("messagebar_active");  // fix the left nav breadcrumb scroll fixed top
                $(window).trigger('content_resize'); // fix the scroll
            });
        }

        var  html = getOverloadFunctionsListHtml(overloadjson);
        $('#responsive_offcanvas').prepend(html);
    }   
}