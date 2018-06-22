'use babel';

(function() {
  this.keywords         = '';
  this.keywordsArray    = [];
  this.functions        = '';
  this.functionsArray   = [];
  this.setKeywords = function (keywords) {
    this.keywords = keywords;
    this.keywordsArray = [];
    var tempArr = keywords.split(' ');
    for (var i in tempArr) {
      this.keywordsArray[i] = {text: tempArr[i], type: 'keyword'};
    }
  };
  this.setFunctions = function (functions) {
    this.functions = functions;
    this.functionsArray = [];
    var tempArr = functions.split(' ');
    for (var i in tempArr) {
      this.functionsArray[i] = {text: tempArr[i], type: 'function'};
    }
  };
  this.getAutoCompleteArr = function (prefix) {
    if(!(/^([a-z0-9_])*$/.test(prefix))){
      return [];
    }
    var completionList = this.keywordsArray.concat(this.functionsArray);
    if(!(completionList.length)){
      return [];
    }
    var reg = new RegExp('^' + prefix + '.*');
    completionList = completionList.filter(function(obj) { return reg.test(obj.text); });
    completionList.forEach(function(element){
      if (typeof element === 'object' ){
        element['replacementPrefix'] = prefix;
      }
    });
    completionList.sort(function(obj1,obj2) {
      var x = obj1['text']; var y = obj2['text'];
      return ((x < y) ? -1 : ((x > y) ? 1 : 0));
    });
    return completionList;
  };
}).call(this);
