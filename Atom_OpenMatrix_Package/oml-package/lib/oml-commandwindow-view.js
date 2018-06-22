'use babel';

var globals    = require('./oml-globals');
var connection = require('./oml-connection');
var cmdHistory = require('./oml-command-history');
export default class OmlCommandWindowView {

  constructor(serializedState) {
    this.txtareaObj = '';
    this.zoom = 100;
    // Create root element
    this.element = document.createElement('TABEL');
    this.element.addEventListener('resize', function () { console.log('resize');});

    this.element.id = 'oml-package-cmdwindow';
    this.element.name = 'oml-package-cmdwindow';
    this.insertBanner();
    this.element.style.overflow = 'scroll';
    this.element.style.zoom = '100%';
    this.element.style.width = '100%';
    this.insertPrompt();
    connection.onWaitngForUserInput(this.waitingForUserInput);
    let self = this;
    this.element.addEventListener('contextmenu',this.addContextMenu);
    this.element.addEventListener('mousewheel',function(evt) {
      if(typeof evt !== 'undefined') {
        if(evt.ctrlKey) {
          if (evt.wheelDelta > 0) {
            self.zoom = self.zoom+5;
            if(self.zoom > 500) {
              self.zoom = 500;
            }
            self.element.style.zoom = self.zoom+'%';
          } else {
            self.zoom = self.zoom-5;
            if(self.zoom < 50) {
              self.zoom = 50;
            }
            self.element.style.zoom = self.zoom+'%';
          }
        }
      }
     });
     atom.window.addEventListener('resize',this.adjustElements);
     this.element.addEventListener('click',function(evt) { 
       let activePrompt = document.getElementById('activeomlinput');
       let rect = activePrompt.getBoundingClientRect();
       if(evt.clientY > rect.bottom) {activePrompt.focus()}
      });
     setTimeout(function() {globals.commandWindowObj.applyFontSize(globals.fontSize);} ,10);
  }
  // Returns an object that can be retrieved when package is activated
  serialize() {}
  // Tear down any state and detach
  destroy() {
    this.element.remove();
  }

  getElement() {
    return this.element;
  }
  getTitle() {
    // Used by Atom for tab text
    return 'OML command window';
  }
  getIconName() {
    return 'terminal';
  }
  getURI() {
    return globals.commandWindowUri;
  }
  getDefaultLocation() {
    return 'bottom';
  }
  getAllowedLocations() {
    return ['left', 'right', 'bottom'];
  }
  insertPrompt () {
    let self = this;
    if ('' !== this.txtareaObj) {
      this.txtareaObj.readOnly = true;
      this.txtareaObj.id  = 'inactiveomlinput';
      this.txtareaObj.name  = 'inactiveomlinput';
      this.txtareaObj.removeEventListener('keydown', handleKeyDownEvent);
      this.txtareaObj.removeEventListener('paste',autosize);
      this.txtareaObj.removeEventListener('keyup', handleKeyUpEvent);
      let activePrompt = document.getElementById('activeomlprompt');
      if(null !== activePrompt) {
        activePrompt.id = 'inactiveomlprompt';
        activePrompt.name = 'inactiveomlprompt';
      }      
   }

    let tr = document.createElement('TR');
    this.element.appendChild(tr);

    let td1 = document.createElement('TD');
    tr.appendChild(td1);

    let prompt = document.createTextNode('>');
    td1.appendChild(prompt);
    td1.id = 'activeomlprompt';
    td1.name = 'activeomlprompt';
    prompt.className = 'native-key-bindings';
    prompt.tabIndex = -1;
    td1.style = 'font-size:'+globals.fontSize+'px;color: #0962af;font-weight: bold;vertical-align:top'

    let td2 = document.createElement('TD');
    tr.appendChild(td2);
    td2.style.width = '100%';

    cmdArea = document.createElement('TEXTAREA');
    td2.appendChild(cmdArea);

    this.txtareaObj = cmdArea;
    this.txtareaObj.id = 'activeomlinput';
    this.txtareaObj.name = 'activeomlinput';
    cmdArea.className = 'native-key-bindings';
    cmdArea.tabIndex = -1;
    cmdArea.style = 'font-size:'+globals.fontSize+'px;letter-spacing:2px;padding-left: 2px;' +
                    'float: left;display: inline-block;vertical-align: top;' +
                    'resize: none;overflow: visible;width: 100%;border-radius: 0px;' +
                    'border-width: 0px;' +
                    'background-color: transparent;border-bottom-color: transparent;';
    cmdArea.rows = 1;

    length = this.txtareaObj.value.length;
    this.txtareaObj.setSelectionRange(length, length);
    this.txtareaObj.focus();

    cmdArea.addEventListener('keydown', handleKeyDownEvent);
    cmdArea.addEventListener('paste',autosize);
    cmdArea.addEventListener('keyup', handleKeyUpEvent);

    function autosize(evt) {
      setTimeout(function(){
        evt.target.style.height = '8px';
        height = (evt.target.scrollHeight)+2;
        evt.target.style.height = height+'px';
        let eid = document.getElementById('oml-package-cmdwindow');
        if (null !== eid) {
          eid.scrollTop = eid.scrollHeight;
        }
      },0);
    }

    function handleKeyDownEvent(evt){
        let key = evt.which || evt.keyCode;
        if (connection.isExecutionPaused()) {
          evt.preventDefault();
          connection.sendStdin('');
        } else {
          var ctrl = evt.ctrlKey ? evt.ctrlKey : ((key === 17) ? true : false);
          if ( key == 67 && ctrl ) {
            //Ctrl+c
            connection.isWaitingForUserInput().then(function() {
              connection.interrputOml();
            }, function () {
              connection.isOmlReady().then(function () {
                  //Do nothing
              },function () {
                connection.interrputOml();
              });
            });
          } else if (key == 13) {
            //Enter key
            /*if(connection.isOmlAlive()) {
              connection.isWaitingForUserInput().then(function() {
                evt.preventDefault();
              }, function () {
                connection.isPartialExpression(evt.target.value).then(function () {
                }, function () {
                  evt.preventDefault();
                });
              });
            } else {*/
              evt.preventDefault();
            /*}*/
          } else if (38 == key) {
            //Up arrow
            let activePrompt = document.getElementById('activeomlinput');
            if (null !== activePrompt) {
              let lineNo = activePrompt.value.substr(0, activePrompt.selectionStart).split('\n').length;
              if (lineNo == 1) {
                let cmd = cmdHistory.getPrevious();
                if ('' !== cmd) {
                  activePrompt.value = cmd;
                }
              }
            }
          } else if (40 == key) {
            //Down arrow
            let activePrompt = document.getElementById('activeomlinput');
            if (null !== activePrompt) {
              let lineNo = activePrompt.value.substr(0, activePrompt.selectionStart).split('\n').length;
              let linesCount = activePrompt.value.split('\n').length;
              if(lineNo == linesCount) {
                let cmd = cmdHistory.getNext();
                if ('' !== cmd) {
                  if (null !== activePrompt) {
                    activePrompt.value = cmd;
                  }
                }
              }
            }
          }
        }
        autosize(evt);
    }

    function handleKeyUpEvent(evt) {
      evt = evt || window.event;
      let eid = document.getElementById('oml-package-cmdwindow');
      if (null !== eid) {
        eid.scrollTop = eid.scrollHeight;
      }
      if (evt.keyCode == 13) {
        if (evt.target.value !== '' && connection.isOmlAlive() && evt.target.id == 'activeomlinput') {
            connection.isWaitingForUserInput().then(function () {
              evt.target.readOnly = true;
              self.insertPrompt();
              connection.sendStdin(evt.target.value);
            }, function () {
                connection.isPartialExpression(evt.target.value).then(function () {
                  let startPos = evt.target.selectionStart;
                  let endPos = evt.target.selectionEnd;
                  let selectionBefore = evt.target.value.substring(0, startPos);
                  let selectionAfter = evt.target.value.substring(endPos);
                  evt.target.value = selectionBefore+'\n'+selectionAfter;
                  evt.target.focus();
                  evt.target.selectionEnd = startPos+1;
                  autosize(evt);
                }, function () {
                  evt.preventDefault();
                  connection.isOmlReady().then(function () {
                    evt.target.readOnly = true;
                    self.insertPrompt();
                    connection.execute('evalstring,' + evt.target.value);
                    connection.startProgress();
                    cmdHistory.addCmdHistory(evt.target.value);
                  },function () {
                    atom.notifications.addInfo('OML interpreter status:', {
                    dismissable: true,
                    description: 'OML interpreter is busy.'
                  });
                });
              });
          });
        }
        return false;
      }
      return true;
    };
  }

  stdOut(data,color='') {
    data = data.replace(/(?:\r\n|\r|\n)/g, '<br />');
    let tr = document.createElement('TR');
    let td1 = document.createElement('TD');
    let td2 = document.createElement('TD');
    let tn = document.createElement('DIV');
    tn.innerHTML = data;
    tr.id = 'omloutput';
    tr.name = 'omloutput';
    tr.className = 'native-key-bindings';
    tr.tabIndex = -1;
    td2.style = 'font-size:'+globals.fontSize+'px;word-break:break-all;overflow: auto;word-wrap: break-word;color:'+color;
    this.element.insertBefore(tr,this.element.lastChild);

    tr.appendChild(td1);
    tr.appendChild(td2);
    td2.appendChild(tn);
    let eid = document.getElementById('oml-package-cmdwindow');
    if (null !== eid) {
      eid.scrollTop = eid.scrollHeight;
    }
  }
  stdErr(data) {
    this.stdOut(data,'#ff2525');
  }
  getElement() {
    return this.element;
  }
  getCurrentTextarea () {
    return this.txtareaObj;
  }
  insertBanner() {
    // Create message element
    let tr = document.createElement('TR');
    let td1 = document.createElement('TD');
    let td2 = document.createElement('TD');
    let tn = document.createTextNode('OML command window.');
    td2.style = 'font-size:'+globals.fontSize+'px;color:'+globals.messageColor;
    tr.id = 'omlbanner';
    tr.name = 'omlbanner';
    tr.className = 'native-key-bindings';
    tr.tabIndex = -1;
    this.element.appendChild(tr);

    tr.appendChild(td1);
    tr.appendChild(td2);
    td2.appendChild(tn);
  }
  clear() {
    while(this.element.hasChildNodes()) {
      this.element.removeChild(this.element.lastChild);
    }
  }
  waitingForUserInput (state) {
    if(null !== document.getElementById('activeomlprompt')) {
      if(state) {
        document.getElementById('activeomlprompt').firstChild.textContent = '<';
      } else {
        document.getElementById('activeomlprompt').firstChild.textContent = '>';
      }
    }
  }
  addContextMenu () {
    let selectedTxt = atom.window.getSelection().toString();
    let isZoomed = ((document.getElementById('oml-package-cmdwindow')).style.zoom !== '100%');
    atom.contextMenu.add({
      '#inactiveomlinput, #inactiveomlprompt, #oml-package-cmdwindow, #omlbanner, #omloutput':
      [
        {label: 'Cut', command: 'core:cut',enabled: 'false'},
        {label: 'Copy', command: 'core:copy',enabled: ('' !== selectedTxt)? true:false},
        {label: 'Paste', command: 'core:paste',enabled: 'false'},
        {type: 'separator'},
        {label: 'Reset Zoom', command: 'oml-commandwindow:resetOmlCommandWindowZoom',enabled: isZoomed},
        {label: 'Clear', command: 'oml-commandwindow:clearOmlCommandWindow'},
        {type: 'separator'}
      ]
    });

    let input = document.getElementById('activeomlinput');
    let isSelction  = (input.selectionStart !== input.selectionEnd);
    let clipTxt = atom.clipboard.read();
    atom.contextMenu.add({
      '#activeomlinput, #activeomlprompt':
      [
        {label: 'Cut', command: 'core:cut',enabled: isSelction},
        {label: 'Copy', command: 'core:copy',enabled: isSelction},
        {label: 'Paste', command: 'core:paste',enabled: ('' !== clipTxt)},
        {type: 'separator'},
        {label: 'Reset Zoom', command: 'oml-commandwindow:resetOmlCommandWindowZoom',enabled: isZoomed},
        {label: 'Clear', command: 'oml-commandwindow:clearOmlCommandWindow'},
        {type: 'separator'}
      ]
    });
  }
  resetZoom () {
    this.zoom = 100;
    this.element.style.zoom = '100%';
  }
  adjustElements () {
    let txtWidgets = [];
    let table = document.getElementById('oml-package-cmdwindow');
    txtWidgets.unshift(table.getElementsByTagName('TEXTAREA'));
    if(txtWidgets[0].length) {
      for (let i = 0; i < txtWidgets[0].length; i++) {      
        setTimeout(function(){
          txtWidgets[0][i].style.height = '8px';
          height = (txtWidgets[0][i].scrollHeight)+2;
          txtWidgets[0][i].style.height = height+'px';
        },0);
      }
    }
    let eid = document.getElementById('oml-package-cmdwindow');
    if (null !== eid) {
      eid.scrollTop = eid.scrollHeight;
    }
  }
  applyFontSize (points) {
    if (points === parseInt(points, 10)) {
      globals.fontSize = points;
      let table = document.getElementById('oml-package-cmdwindow');
      if(null !== table) {
        let children = table.getElementsByTagName('td');
  
        for (let i = 0; i < children.length; i++) {
          children[i].style.fontSize = globals.fontSize+'px'
        }
        children = table.getElementsByTagName('textarea')
        for (let i = 0; i < children.length; i++) {
          children[i].style.fontSize = globals.fontSize+'px'
        }
        globals.commandWindowObj.adjustElements();
      }
    }
  }
}
