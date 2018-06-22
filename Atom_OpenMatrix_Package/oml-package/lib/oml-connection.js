'use babel';

var globals      = require('./oml-globals');
var progress     = require('./oml-progress');
var completionData = require('./oml-completiondata.js');
var fs           = require('fs');
var os           = require('os');
var path         = require('path');
var childProcess = require('child_process');
var net          = require('net');
var DEFAULTPORT  = 9099;
var port         = this.DEFAULTPORT;
var isOmlRunning = false;
var userInputInProgress = false;
var isPause = false;
var childObject = {};
var Emitter      = require('atom').Emitter;
var connectionRetryCount = 0;
var isSocketConnected = false;
//possibe states: notstarted ready executing
var omlInterPreterStatus = 'notstarted';

(function() {
  module.exports = {
    emitter: new Emitter,
    clientSocket: '',
    getFreePort: function(self) {
      return new Promise(function(resolve) {
        let server = net.createServer();
        return server.listen(0, '127.0.0.1', function() {
          var port;
          port = server.address().port;
          server.close();
          return resolve([port,self]);
        });
      });
    },
    startOml: function() {
      if (!this.isOmlAlive())  {
        (this.getFreePort(this)).then(function([fPort, self]){
          self.clientSocket = new net.Socket();
          if(isNaN(fPort)){
            port = DEFAULTPORT;
          } else {
            port = fPort;
          }
          let omlPath = '';
          if(null !== document.getElementById('omlsettings-exepath')) {
            omlPath = path.normalize(document.getElementById('omlsettings-exepath').value);
          } else {
            omlPath = globals.omlExePath;
          }
          let omlToolboxPath = '';
          if(null !== document.getElementById('omlsettings-exepath')) {
            omlToolboxPath = path.normalize(document.getElementById('omlsettings-toolboxpath').value);
          } else {
            omlToolboxPath = globals.omlToolboxPath;
          } 
          fs.stat(omlPath, function (err, stats) {
            if (err) {
              atom.notifications.addError('Can not start OML', { dismissable: true,
                description: 'Not a valid OML executable path. It can be changed from OML settings.'});
            } else if (stats.isFile()) {
              let omlrootdir =  path.dirname(omlPath);
              if ('linux' ===  os.platform()) {
                process.env['LD_LIBRARY_PATH'] = process.env['LD_LIBRARY_PATH']+":"+omlrootdir;
              }
              let args =  ['-atom','-port',port]
              fs.stat(omlToolboxPath, function (err, stats) {
                if (err) {                                                                
                  atom.notifications.addError('Can not load tool boxes', { dismissable: true,
                    description: 'Not a valid OML file path. It can be changed from OML settings.'});
                } else if (stats.isFile()) {
                    args = ['-atom','-port',port,'-f',omlToolboxPath]
                }
                
                childObject = childProcess.spawn(omlPath, args);              
                if (null === globals.commandWindowObj) {
                  atom.workspace.open(globals.commandWindowUri);
                }
                self.startProgress();
                if (childObject.pid) {
                  childObject.stdin.setEncoding('utf8');
                  childObject.stdout.setEncoding('utf8');
                  childObject.stderr.setEncoding('utf8');
                  childObject.stdout.on('data', function(data){
                    globals.commandWindowObj.stdOut(data);
                  });
  
                  childObject.stderr.on('data', function(data){
                    globals.commandWindowObj.stdErr(data);
                  });
  
                  childObject.on('exit', function (code) {
                    console.log('OML interpreter exited with code ' + code);
                    globals.commandWindowObj.stdOut('OML has stopped.',globals.messageColor);
                    isOmlRunning = false;
                    userInputInProgress = false;
                    isPause = false;
                    omlInterPreterStatus = 'notstarted';
                    self.stopProgress();
                  });
                  isOmlRunning = true;
                  connectionRetryCount = 1;
                  isSocketConnected = false;
                  globals.commandWindowObj.stdOut('OML started.',globals.messageColor);
  
                  setTimeout(function () {self.connect(self);},2000);
  
  
                } else {
                  console.log('Failed to creare OML interpreter.');
                }
              });
            } else {
              atom.notifications.addError('Can not start OML', { dismissable: true,
                description: 'Not a valid OML executable path. It can be changed from OML settings.'});
            }
          });
        });
      }  else {
        console.log('OML interpreter already running. Can not start again.');
      }
    },
    interrputOml: function() {
      let self = this;
      if (this.isOmlAlive())  {
        self.isOmlReady().then(function () {
          self.stopProgress();
        }, function () {
          self.execute('interrupt,');
          if(userInputInProgress || isPause) {
            self.sendStdin('');
          }
        });
      }
    },
    stopOml: function() {
      if (!this.isOmlAlive())  {
        console.log('OML interpreter is not running.');
      }  else {
        childObject.kill();         
        isOmlRunning = false;
        userInputInProgress = false;
        isSocketConnected = false;
        isPause = false;
        omlInterPreterStatus = 'notstarted';
      }
    },
    execute: function(cmd) {
      if(this.isOmlAlive()) {
        this.clientSocket.write(cmd);
        this.startProgress();
      } else {
        console.info('OML interpreter is not running.');
      }
    },
    isOmlAlive: function() {
      return isOmlRunning;
    },
    isOmlReady: function() {
      return new Promise(function(resolve,reject) {
        if('ready' === omlInterPreterStatus){
          return resolve();
        } else {
          return reject();
        }
      });
    },
    isWaitingForUserInput: function() {
      return new Promise(function(resolve,reject) {
        if(userInputInProgress){
          return resolve();
        } else {
          return reject();
        }
      });
    },
    isExecutionPaused: function () {
      return isPause;
    },
    isPartialExpression(expression) {
      if(this.isOmlAlive()) {
        this.clientSocket.write('ispartialexpression,'+expression);
        let self = this;
        return new Promise(function(resolve,reject) {
          self.onPartialExpressionCheck(function(isPartialExpression) {
            if(isPartialExpression) {
              resolve();
            } else {
              reject();
            }
          });
        });
      } else {
        return new Promise(function(resolve,reject) {resolve();});
      }
    },
    connect: function (self) {
      if(self.isOmlAlive()) {
        self.clientSocket.connect({'port':port, host:'127.0.0.1'});

        self.clientSocket.on('connect', function () {
           isSocketConnected = true;
           console.log('OML server socket connected.');
        });
        
        self.clientSocket.on('data', function (data) {
          let jObjs = JSON.parse('[' + (data.toString()).replace(/}{/g, '},{') + ']');
          for (var i = 0; i < jObjs.length; i++) {
            if('' === jObjs[i]) {
              continue;
            }
            switch (jObjs[i].type) {
              case 'status': {
                omlInterPreterStatus = jObjs[i].data;
                if('ready' === omlInterPreterStatus){
                  self.stopProgress();
                  isSocketConnected = true;
                }
                break;
              } 
              case 'keywordslist': {
                completionData.setKeywords(jObjs[i].data);
                break;
              }
              case 'functionslist': {
                completionData.setFunctions(jObjs[i].data);
                break;
              }
              case 'userinput': {
                userInputInProgress = jObjs[i].data;
                self.emitter.emit('onWaitngForUserInput',userInputInProgress);
                break;
              }
              case 'pause': {
                isPause = jObjs[i].data;
                break;
              }
              case 'partialexpression': {
                self.emitter.emit('onPartialExpressionCheck',jObjs[i].data);
                break;            
              }
              case 'control': {
                switch (jObjs[i].data) {
                  case 'clearcmdwindow':
                    if (null !== globals.commandWindowObj) {
                      globals.commandWindowObj.clear();
                      globals.commandWindowObj.insertBanner();
                      globals.commandWindowObj.insertPrompt();
                    }
                    break;
                }
                break;              
              }
            }
          }
        });
  
        self.clientSocket.on('error', function (error) {
          console.log('OML server socket connection error.',error);
          if (self.isOmlAlive() && !isSocketConnected) {
            if(connectionRetryCount < globals.connectionRetryLimit) {
              self.clientSocket.end();
              self.clientSocket = new net.Socket();
              console.log('Retrying to connect OML server.');
              setTimeout(function () {self.connect(self);},1000);
              connectionRetryCount = connectionRetryCount + 1;
            } else {
              console.log('Exceeded maximum number of retries. Stopping OML server.');
              isSocketConnected = false;
              childObject.kill();            
            }
          }            
        });
        
        self.clientSocket.on('end', function () {
          isSocketConnected = false;
          console.log('OML server socket closed.');
        });
      }
    },
    updateProgress: function (self) {
      if(self.isOmlAlive()) {
        self.isOmlReady().then(function () {
          self.stopProgress();
        },function () {
          setTimeout(function() {self.updateProgress(self);},1);
        });
      } else {
        self.stopProgress();
      }
    },
    startProgress: function () {
        progress.addProgress();
    },
    stopProgress: function () {
      progress.removeProgress();
    },
    onInterpreterReady: function(callback) {
      return this.emitter.on('interpreterReady', callback);
    },
    onWaitngForUserInput: function(callback) {
      return this.emitter.on('onWaitngForUserInput', callback);
    },
    onPartialExpressionCheck: function(callback) {
      return this.emitter.on('onPartialExpressionCheck',callback);
    },
    sendStdin: function(data) {
      childObject.stdin.write(data+'\n');
    },
    deactivate: function () {
      this.emitter.dispose();
    }
  };
}).call(this);
