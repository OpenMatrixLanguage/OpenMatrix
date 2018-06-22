'use babel';

var connection = require('./oml-connection');
var omlCompletiondata    = require('./oml-completiondata');
(function() {
  this.provider = {
    selector: '.source.oml',
    disableForSelector: '.source.oml .comment',
    inclusionPriority: 1,
    excludeLowerPriority: false,
    suggestionPriority: 2,
    getSuggestions: function(arg) {
      if(connection.isOmlAlive()) {
        return new Promise(function(resolve) {
          return resolve((omlCompletiondata.getAutoCompleteArr(arg.prefix)));
        });
      }
      else {
        return [];
      }
    }
  };
}).call(this);
