'use babel';

var connection = require('./oml-connection');

(function() {
  module.exports = {
    runOmlFile: function() {
      let editor = atom.workspace.getActiveTextEditor();
      if (typeof editor !== 'undefined') {
        let filepath = editor.getPath();
        if (typeof filepath !== 'undefined') {
          let fileExtension = filepath.substr(filepath.length - 4, 4).toLowerCase();
          let isOmlFile = true;
          if ('.oml' !== fileExtension) {
            isOmlFile = false;
          }
          if(isOmlFile) {
            if(connection.isOmlAlive()) {
              connection.isOmlReady().then(function () {
                connection.execute('evalfile,'+filepath);
                connection.startProgress();
              }, function () {
                atom.notifications.addInfo('OML interpreter status:', {
                  dismissable: true,
                  description: 'OML interpreter is busy.'
                });
              });
            }
          } else {
            console.error('Not a valid OML file.File extension must be .oml');
          }
        } else {
          console.info('File not exist.');
        }
      } else {
        console.info('No active editor found.');
      }
    },
    runOmlScript: function() {
      editor = atom.workspace.getActiveTextEditor();
      if (typeof editor !== 'undefined') {
        let selection = editor.getLastSelection();
        if (typeof selection !== 'undefined') {
          let text = selection.getText();
          if (text !== '' && connection.isOmlAlive()) {
            connection.isOmlReady().then(function () {
              connection.execute('evalstring,' + text);
              connection.startProgress();
            },function () {
              atom.notifications.addInfo('OML interpreter status:', {
                dismissable: true,
                description: 'OML interpreter is busy.'
              });
            });
          }
        }
      } else {
        console.info('No active editor found.');
      }
    },
  };
}).call(this);
