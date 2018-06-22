'use babel';

var globals       = require('./oml-globals');

(function() {
  module.exports = {
    openOmlCommandWindow() {
      atom.workspace.open(globals.commandWindowUri);
    },
    clearOmlCommandWindow() {
      if (null !== globals.commandWindowObj) {
        globals.commandWindowObj.clear();
        globals.commandWindowObj.insertBanner();
        globals.commandWindowObj.insertPrompt();
      }
    },
    resetOmlCommandWindowZoom() {
      globals.commandWindowObj.resetZoom();
    }
  };
}).call(this);
