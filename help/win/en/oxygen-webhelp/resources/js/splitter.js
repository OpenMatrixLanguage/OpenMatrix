/**
 * jQuery.splitter.js - animated splitter plugin
 *
 * version 1.0 (2010/01/02) 
 *
 *Licensed under the MIT license:
*   http://www.opensource.org/licenses/mit-license.php
*
* Code cleaned up and made to pass jslint <scudette@google.com>. Also
* changed default behaviour to not animate when dragged.
 */
 
/**
* jQuery.splitter() plugin implements a two-pane resizable animated window, using existing DIV elements for layout.
* For more details and demo visit: http://krikus.com/js/splitter
*
* @example $("#splitterContainer").splitter({splitVertical:true,A:$('#leftPane'),B:$('#rightPane'),closeableto:0});
* @desc Create a vertical splitter with toggle button
*
* @example $("#splitterContainer").splitter({minAsize:100,maxAsize:300,splitVertical:true,A:$('#leftPane'),B:$('#rightPane'),slave:$("#rightSplitterContainer"),closeableto:0});
* @desc Create a vertical splitter with toggle button, with minimum and maximum width for plane A and bind resize event to the slave element
*
* @name splitter
* @type jQuery
* @param Object options Options for the splitter ( required)
* @cat Plugins/Splitter
* @return jQuery
* @author Kristaps Kukurs (contact@krikus.com)
*/

/*
 * Patches added by Oxygen Webhelp plugin are marked with the comments 
 * "oXygen patch start" and "oXygen patch end". The patches are the following:
 * 
 * 1. Added attribute class="splitterMask" to right side iframe
 * 2. Applied a correction to the computed size (width) of the right side iframe.
 * 3. Modified values of overflow and overflow-y CSS properties of left side iframe.
 * 4. Add URL parameter for initial state of TOC - toc.visible: true - will be visible (default) / false - TOC will not be visible.
 * 
 */

;
(function($){
    $.fn.splitter = function(args){
        args = args || {};
        return this.each(function() {
            var _ghost;		//splitbar  ghosted element 
            var splitPos;	 // current splitting position
            var _splitPos;	 // saved splitting position
            var _initPos;	//initial mouse position
            var _ismovingNow=false;	// animaton state flag

            // Default opts
            var direction = (args.splitHorizontal? 'h':'v');
            var opts = $.extend({
                minAsize:0, //minimum width/height in PX of the first (A) div.
                maxAsize:0, //maximum width/height  in PX of the first (A) div.
                minBsize:0, //minimum width/height in PX of the second (B) div.
                maxBsize:0, //maximum width/height  in PX of the second (B) div.
                ghostClass: 'working',// class name for _ghosted splitter and hovered button
                invertClass: 'invert',//class name for invert splitter button                
                animSpeed: 250 //animation speed in ms
            },{
                v:{ // Vertical
                    moving:"left",
                    sizing: "width", 
                    eventPos: "pageX",
                    splitbarClass:"splitbarV",
                    buttonClass: "splitbuttonV", 
                    cursor: "e-resize"
                },
                h: { // Horizontal 
                    moving:"top",
                    sizing: "height", 
                    eventPos: "pageY",
                    splitbarClass:"splitbarH",
                    buttonClass: "splitbuttonH",  
                    cursor: "n-resize"
                }
            }
            [direction], args);

            //setup elements
            var splitter = $(this);
            var mychilds =splitter.children(); //$(">*", splitter[0]);
            var A = args.A;	// left/top frame
            var B = args.B;// right/bottom frame
            var slave=args.slave;//optional, elemt forced to receive resize event
            //Create splitbar 
            var C=$('<div><span></span></div>');
            A.after(C);
            C.attr({
                "class": opts.splitbarClass,
                unselectable:"on"
            }).css({
                "cursor":opts.cursor,
                "user-select": "none", 
                "-webkit-user-select": "none",
                "-khtml-user-select": "none", 
                "-moz-user-select": "none"
            })
            .bind("mousedown", startDrag);
 		
            if(opts.closeableto!=undefined){
                var Bt=$('<div></div>').css("cursor",'pointer');
                C.append(Bt);
                Bt.attr({
                    "class": opts.buttonClass, 
                    unselectable: "on"
                });
                Bt.hover(function(){
                    $(this).addClass(opts.ghostClass);
                },function(){
                    $(this).removeClass(opts.ghostClass);
                });
                Bt.mousedown(function(e){
                    if(e.target != this)return;
                    Bt.toggleClass(opts.invertClass).hide();
                    splitTo((splitPos==opts.closeableto)?_splitPos:opts.closeableto,true);
                    return false;
                });
            }		
            //reset size to default.			
            var perc=(((C.position()[opts.moving]-splitter.offset()[opts.moving])/splitter[opts.sizing]())*100).toFixed(1);

            // oXygen patch start
            visibleTOC = location.search.indexOf('toc.visible=false') != -1 ? false : true;
            splitPos = visibleTOC ? perc : opts.minAsize;
            _splitPos = visibleTOC ? opts.minAsize : perc;
            var initialSplitTo = visibleTOC ? perc : opts.minAsize;
            splitTo(initialSplitTo,false,true);

            if( initialSplitTo == opts.minAsize ) {
                Bt.toggleClass(opts.invertClass);
            }
            // oXygen patch end
            
            // resize  event handlers;
            splitter.bind("resize",function(e, size){
                if(e.target!=this)return;
                splitTo(splitPos,false,true);
            });
            $(window).bind("resize",function(){
                splitTo(splitPos,false,true);
            });

            //C.onmousedown=startDrag
            function startDrag(e) {
                if(e.target != this)return;
                _ghost = _ghost || C.clone(false).insertAfter(A);
                splitter._initPos=C.position();
                splitter._initPos[opts.moving]-=C[opts.sizing]();
                _ghost.addClass(opts.ghostClass).css('position','absolute').css('top', C.offset().top+"px").css('z-index','250').css("-webkit-user-select", "none").width(C.width()).height(C.height()).css(opts.moving,splitter._initPos[opts.moving]);
                mychilds.css("-webkit-user-select", "none");	// Safari selects A/B text on a move
                // oXygen patch start
                // in order to work with right iframe
                /* need this CSS also      
                 * .splitterMask{
                 * position:absolute;
                 * top: 0;
                 * left: 0;
                 * width: 100%;
                 * height: 100%;
                 * overflow: hidden;
                 * //background-image: url(/images/spacer.gif);
                 * z-index: 20000;
                 * }
                 */
                $('<div class="splitterMask"></div>').insertAfter(C);
                // oXygen patch end
                A._posSplit = e[opts.eventPos];

                $(document).bind("mousemove", performDrag).bind("mouseup", endDrag);
            }
            //document.onmousemove=performDrag
            function performDrag(e) {
                if (!_ghost||!A) return;
                var incr = e[opts.eventPos]-A._posSplit;
                _ghost.css(opts.moving, splitter._initPos[opts.moving]+incr);
            }
            //C.onmouseup=endDrag			
            function endDrag(e){
              // oXygen patch start
              // in order to work with right iframe
              $('div.splitterMask').remove();
              // oXygen patch end
                
                var p=_ghost.position();
                _ghost.remove();
                _ghost = null;	
                mychilds.css("-webkit-user-select", "text");// let Safari select text again
                $(document).unbind("mousemove", performDrag).unbind("mouseup", endDrag);
                var perc=(((p[opts.moving]-splitter.offset()[opts.moving])/splitter[opts.sizing]())*100).toFixed(1);		 
                splitTo(perc,(splitter._initPos[opts.moving]>p[opts.moving]),false); 
                splitter._initPos=0;
            }
            //Perform actual splitting and animate it;					
            function splitTo(perc,reversedorder,fast) {
                if(_ismovingNow||perc==undefined)return;//generally MSIE problem
                _ismovingNow=true;
                if(splitPos&&splitPos>10&&splitPos<90)//do not save accidental events
                    _splitPos=splitPos;
                splitPos=perc;

                var barsize=C[opts.sizing]()+(2*parseInt(C.css('border-'+opts.moving+'-width')));//+ border. cehap&dirty
                var splitsize=splitter[opts.sizing]()-barsize;
                
                if(opts.closeableto!=perc){
                    var percpx=Math.max(parseInt((splitsize/100)*perc),opts.minAsize);
                    if(opts.maxAsize)percpx=Math.min(percpx,opts.maxAsize);
                }else{
                    var percpx=parseInt((splitsize/100)*perc,0);
                }
                if(opts.maxBsize){
                    if((splitsize-percpx)>opts.maxBsize)
                        percpx=splitsize-opts.maxBsize;
                }
                if(opts.minBsize){
                    if((splitsize-percpx)<opts.minBsize)
                        percpx=splitsize-opts.minBsize;
                }
                var sizeA=Math.max(0,(percpx-barsize));                
                // oXygen patch start
                var percpxx=percpx;
                if (sizeA==0){
                    percpxx=percpxx+barsize;
                }
                // oXygen patch end
                var sizeB=Math.max(0,(splitsize-percpxx));                        
                splitsize=(splitsize-barsize);                    
                
                //A.attr('title','- '+sizeA);  B.attr('title','- '+sizeB);
                // oXygen patch start
                if (typeof eval(debug) == 'function') {
                  //debug('Resize splitter to; A:'+sizeA+' B:'+sizeB);
                }
//                if (sizeA==0){
//                  A.css('overflow','hidden');
//                  A.css('overflow-y','hidden');
//                }else{
//                    if ($("#contentBlock").is(":visible")) {
//                        A.css('overflow', 'auto');
//                        A.css('overflow-y', 'auto');
//                    } else {
//                        A.css('overflow','hidden');
//                        A.css('overflow-y','hidden');
//                    }
//                }
                // oXygen patch end
                if(fast){
                    A.show().css(opts.sizing,sizeA+'px');
                    B.show().css(opts.sizing,sizeB+'px');
                    Bt.show();
                    if (BrowserDetect.browser!='Explorer' ){
                        mychilds.trigger("resize");
                        if(slave)slave.trigger("resize");
                    }
                    _ismovingNow=false;	
                    return true;
                }				
                if(reversedorder){//reduces flickering if total percentage becomes more  than 100 (possible while animating)
                    var anob={};
                    anob[opts.sizing]=sizeA+'px';
                    A.show().animate(anob,opts.animSpeed,function(){
                        Bt.fadeIn('fast');
                        if($(this)[opts.sizing]()<2){
                            this.style.display='none';
                            B.stop(true,true);
                            B[opts.sizing](splitsize+'px');
                        }
                    });
                var anob2={};
                anob2[opts.sizing]=sizeB+'px';
                B.show().animate(anob2,opts.animSpeed,function(){
                    Bt.fadeIn('fast');
                    if($(this)[opts.sizing]()<2){
                        this.style.display='none';
                        A.stop(true,true);
                        A[opts.sizing](splitsize+'px')
                        }
                    });
        }else{
            var anob={};
            anob[opts.sizing]=sizeB+'px';
            B.show().animate(anob,opts.animSpeed,function(){
                Bt.fadeIn('fast');
                if($(this)[opts.sizing]()<2){
                    this.style.display='none';
                    A.stop(true,true);
                    A[opts.sizing](splitsize+'px')
                    }
                });
        var anob={};
        anob[opts.sizing]=sizeA+'px';
        A.show().animate(anob,opts.animSpeed,function(){
            Bt.fadeIn('fast');
            if($(this)[opts.sizing]()<2){
                this.style.display='none';
                B.stop(true,true);
                B[opts.sizing](splitsize+'px');
            }
        });					
    }	
    //trigger resize evt
    splitter.queue(function(){  
        setTimeout(function(){  
            splitter.dequeue();
            _ismovingNow=false;
            mychilds.trigger("resize");
            if(slave)slave.trigger("resize");		
        }, opts.animSpeed+5);  
    });
 
}//end splitTo()
			
			
	


});//end each
};//end splitter

})(jQuery);