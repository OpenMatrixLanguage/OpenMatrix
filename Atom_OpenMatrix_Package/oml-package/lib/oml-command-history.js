'use babel';

(function() {
  this.cmdHistory = [];
  this.cmdIndex = 0;
  this.addCmdHistory = function (cmd) {this.cmdHistory.unshift(cmd);this.cmdIndex = 0};
  this.getPrevious = function () {
                        this.cmdIndex = this.cmdIndex+1;
                        if(this.cmdIndex <= this.cmdHistory.length) {
                          return this.cmdHistory[this.cmdIndex-1];
                        } else {
                          this.cmdIndex = this.cmdHistory.length;
                          return '';
                        }                        
                     };
  this.getNext = function () {
                       this.cmdIndex = this.cmdIndex-1;
                       if(this.cmdIndex > 0 && this.cmdIndex <= this.cmdHistory.length) {
                         return this.cmdHistory[this.cmdIndex-1];
                       } else {
                         this.cmdIndex = this.cmdIndex+1;
                         return '';
                       }
                     };
  this.clear = function () {this.cmdHistory = [];this.cmdIndex = 0};
}).call(this);
