$(document).ready(function () {
    registerMatlabCommandDialogAction();
}); 

$(window).bind('examples_cards_added', function(e) {
  $('.card_container a[href^="matlab:"]').off('click');
  registerMatlabCommandDialogAction();
});

function registerMatlabCommandDialogAction() {
  $('a[href^="matlab:"]').on('click', function (e) {
      e.preventDefault();
      var href = $(this).attr('href'),
      match = href.match(/matlab:(.*)/),
      matlabCommand = null;

      if (match) {
      matlabCommand = match[1];
      }

      if (matlabCommand) {
      $("#matlab-command-dialog #dialog-body #dialog-matlab-command").text(matlabCommand);
      } else {
      $("#matlab-command-dialog #dialog-body #dialog-matlab-command").hide();
      }
      $("#matlab-command-dialog").modal();
  });
};        