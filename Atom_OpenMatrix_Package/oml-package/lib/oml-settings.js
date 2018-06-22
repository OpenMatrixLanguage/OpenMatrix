'use babel';

var globals = require('./oml-globals');
(function() {
  module.exports = {
    openOmlSettingsWindow () {
      atom.workspace.open(globals.omlSettingsUri);
    }
  };
}).call(this);
