'use babel';

var globals = require('./oml-globals');
var statusbarItem = document.createElement('PROGRESS');

(function() {
  module.exports = {
    statusbar: '',
    addProgress() {
      this.statusbar = globals.statusBar.addLeftTile({
        item: statusbarItem,
        priority: -2
      });
    },
    removeProgress() {
      if (this.statusbar !== '') {
        this.statusbar.destroy();
        this.statusbar.status = '';
      }
    }
  };
}).call(this);
