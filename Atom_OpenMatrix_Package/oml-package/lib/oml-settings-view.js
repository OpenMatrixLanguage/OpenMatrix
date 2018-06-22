'use babel';

var globals = require('./oml-globals');

export default class OmlSettingsView {

  constructor(serializedState) {
    // Create root element
    this.element = document.createElement('div');
    this.element.classList.add('oml-package-settings');
    // Create message element
    this.element.style.overflow = 'scroll';

    h1   = document.createElement('H1');
    txt  = document.createTextNode('OML executable path');
    h1.style = 'margin-left: 20px;resize: none;float: left;display: inline-block;' +
                'vertical-align: top;overflow: visible;width: 95%'
    h1.appendChild(txt);
    this.element.appendChild(h1);

    omlPath = document.createElement('TEXTAREA');
    omlPath.placeholder = 'Example: C:/OpenMatrix/bin/win64/omlconsole.exe'
    omlPath.style = 'margin-left: 20px;letter-spacing:2px;font-size:14px;' +
                    'float: left;display: inline-block;vertical-align: top;' +
                    'resize: none;overflow: visible;width: 95%;border-radius: 1px;' +
                    'border-width: 1px 1px 1px;min-height: 14px;line-height: 14px;' +
                    'background-color: transparent;';
    omlPath.id  = 'omlsettings-exepath';
    omlPath.tabIndex = -1;
    omlPath.className = 'native-key-bindings';
    omlPath.setAttribute('type', 'text');
    omlPath.value = globals.omlExePath;
    this.element.appendChild(omlPath);
    omlPath.addEventListener('keydown', ignoreEnter);
    omlPath.rows = 1;

    h2   = document.createElement('H1');
    txt  = document.createTextNode('OML toolboxes file path');
    h2.style = 'margin-left: 20px;resize: none;float: left;display: inline-block;' +
                'vertical-align: top;overflow: visible;width: 95%'
    h2.appendChild(txt);
    this.element.appendChild(h2);

    omlToolboxPath = document.createElement('TEXTAREA');
    omlToolboxPath.placeholder = 'Example: C:/OpenMatrix/scripts/loadtoolboxes.oml'
    omlToolboxPath.style = 'margin-left: 20px;letter-spacing:2px;font-size:14px;' +
                    'float: left;display: inline-block;vertical-align: top;' +
                    'resize: none;overflow: visible;width: 95%;border-radius: 1px;' +
                    'border-width: 1px 1px 1px;min-height: 14px;line-height: 14px;' +
                    'background-color: transparent;';
    omlToolboxPath.id  = 'omlsettings-toolboxpath';
    omlToolboxPath.tabIndex = -1;
    omlToolboxPath.className = 'native-key-bindings';
    omlToolboxPath.setAttribute('type', 'text');
    omlToolboxPath.value = globals.omlToolboxPath;
    this.element.appendChild(omlToolboxPath);
    omlToolboxPath.addEventListener('keydown', ignoreEnter);
    omlToolboxPath.rows = 1;
    function ignoreEnter(evt){
      var el = this;
        if (evt.keyCode == 13) {
          evt.preventDefault();
        }
    }
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
    return 'OML settings';
  }
  getURI() {
    return globals.omlSettingsUri;
  }
  getIconName() {
    return 'tools';
  }
  getDefaultLocation() {
    return 'center';
  }
  getAllowedLocations() {
    return ['left', 'right', 'bottom','center'];
  }
  getElement () {
    return this.element;
  }
}
