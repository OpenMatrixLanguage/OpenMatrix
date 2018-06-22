'use babel';

import OmlCommandWindowView from './oml-commandwindow-view';
import OmlSettingsView from './oml-settings-view';
import { CompositeDisposable } from 'atom';
var globals       = require('./oml-globals');
var connection    = require('./oml-connection');
var execution     = require('./oml-execution');
var menus         = require('./oml-menu');
var commandWindow = require('./oml-command-window');
var settings      = require('./oml-settings');
var help          = require('./oml-help');
var omlCompletionsProvider = require('./oml-completions-provider');
var commandHistory = require('./oml-command-history');

export default {
  subscriptions: null,
  activate(state) {
    if(state.omlexepath){
      globals.omlExePath = state.omlexepath;
    }
    if(state.omlToolboxPath){
      globals.omlToolboxPath = state.omlToolboxPath;
    }
    if(state.commandHistory){
      commandHistory.cmdHistory = state.commandHistory;
    }    
    // Events subscribed to in atom's system can be easily cleaned up with a CompositeDisposable
    this.subscriptions = new CompositeDisposable();
    this.subscriptions.add(atom.commands.add('atom-workspace', {
      'oml-connection:startOml': () => connection.startOml()
    }));
    this.subscriptions.add(atom.commands.add('atom-workspace', {
      'oml-connection:interruptOml': () => connection.interrputOml()
    }));
    this.subscriptions.add(atom.commands.add('atom-workspace', {
      'oml-connection:stopOml': () => connection.stopOml()
    }));
    this.subscriptions.add(atom.commands.add('atom-workspace', {
      'oml-execution:runOmlFile': () => execution.runOmlFile()
    }));
    this.subscriptions.add(atom.commands.add('atom-workspace', {
      'oml-execution:runOmlScript': () => execution.runOmlScript()
    }));
    this.subscriptions.add(atom.commands.add('atom-workspace', {
      'oml-commandwindow:openOmlCommandWindow': () => commandWindow.openOmlCommandWindow()
    }));
    this.subscriptions.add(atom.commands.add('atom-workspace', {
      'oml-commandwindow:clearOmlCommandWindow': () => commandWindow.clearOmlCommandWindow()
    }));
    this.subscriptions.add(atom.commands.add('atom-workspace', {
      'oml-settingswindow:openOmlSettingsWindow': () => settings.openOmlSettingsWindow()
    }));
    this.subscriptions.add(atom.commands.add('atom-workspace', {
      'oml-commandwindow:resetOmlCommandWindowZoom': () => commandWindow.resetOmlCommandWindowZoom()
    }));
    this.subscriptions.add(atom.commands.add('atom-workspace', {
      'oml-help:openOmlHelp': () => help.openOmlHelp()
    }));
    
    atom.workspace.addOpener(uri => {
      if (globals.commandWindowUri === uri) {
        if (null === globals.commandWindowObj) {
          globals.commandWindowObj = new OmlCommandWindowView();
        }
        return globals.commandWindowObj;
      }
    });

    atom.workspace.addOpener(uri => {
      if (globals.omlSettingsUri === uri) {
        if ('' === globals.omlSettingsObj) {
          globals.omlSettingsObj = new OmlSettingsView();
        }
        return globals.omlSettingsObj;
      }
    });
    menus.activate();
  },

  deactivate() {
    connection.deactivate();
    this.subscriptions.dispose();
    menus.deactivate();
  },
  serialize () {
    let omlPath = globals.omlExePath;
    let omlToolboxPath = globals.omlToolboxPath;
    let cmdHistory = commandHistory.cmdHistory;
    if(null !== document.getElementById('omlsettings-exepath')) {
      omlPath = document.getElementById('omlsettings-exepath').value;
    }
    if(null !== document.getElementById('omlsettings-toolboxpath')) {
      omlToolboxPath = document.getElementById('omlsettings-toolboxpath').value;
    }
    return {
      omlexepath: omlPath,
      commandHistory: cmdHistory,
      omlToolboxPath: omlToolboxPath
    };
  },
  consumeStatusBar: function(statusBar) {
    globals.statusBar = statusBar;
  },
  omlCompletions: function() {
    return omlCompletionsProvider.provider;
  }
};
