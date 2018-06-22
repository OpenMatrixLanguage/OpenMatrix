'use babel';

(function() {
  var CompositeDisposable;

  CompositeDisposable = require('atom').CompositeDisposable;

  module.exports = {
    activate: function() {
      var menu;
      this.subs = new CompositeDisposable;
      this.subs.add(atom.menu.add([
        {
          label: 'Packages',
          submenu: this.menu
        }
      ]));

      this.subs.add = atom.menu.add(this.menu);
      menu = atom.menu.template.pop();
      atom.menu.template.splice(3, 0, menu);

    },
    deactivate: function() {
      this.subs.dispose();
    },
    menu: [
      {
        'label': 'OML',
        'submenu': [
          {
            'label': 'Start OML',
            'command': 'oml-connection:startOml'
          },{
            'label': 'Interrupt OML',
            'command': 'oml-connection:interruptOml'
          },{
            'label': 'Stop OML',
            'command': 'oml-connection:stopOml'
          },{
            'label': 'Run File',
            'command': 'oml-execution:runOmlFile'
          },{
            'label': 'Run Selection',
            'command': 'oml-execution:runOmlScript'
          },{
            'label': 'Open command window',
            'command': 'oml-commandwindow:openOmlCommandWindow'
          },{
            'label': 'Clear command window',
            'command': 'oml-commandwindow:clearOmlCommandWindow'
          },{
            'label': 'Settings',
            'command': 'oml-settingswindow:openOmlSettingsWindow'
          },{
            'label': 'Help',
            'command': 'oml-help:openOmlHelp'
          }
        ]
      }
    ]
  };
}).call(this);
