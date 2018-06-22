'use babel';
(function() {
  module.exports = {
    openOmlHelp () {
      atom.applicationDelegate.openExternal('https://www.openmatrix.org');
    }
  };
}).call(this);