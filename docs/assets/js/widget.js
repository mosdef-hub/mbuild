function loadWidget(dataUrl, widgetSelector, width, height) {
  $.getJSON(dataUrl, function(data) {

    if(!width) { width = 300; }
    if(!height) { height = 300; }
       //console.log(data);
       var HEIGHT = height,
           WIDTH = width,
           HEIGHT_PX = height+'px',
           WIDTH_PX = width+'px';
       var canvas = $("<canvas/>").height(HEIGHT).width(WIDTH);
       var iv = new iview(canvas);
       var container = $('<div/>').css({width: HEIGHT_PX, height: WIDTH_PX});

       container.append(canvas);
       $(widgetSelector).append(container).css({
           'margin-left': 'auto',
           'margin-right': 'auto',
           'width': WIDTH
       });
       var options = {
           'camera': 'perspective',
           'background': 'white',
           'colorBy': 'spectrum',
           'primaryStructure': 'nothing',
           'secondaryStructure': 'cylinder & plate',
           'surface': 'nothing'
       };
       iv.loadTopology(data.topology);
       iv.loadCoordinates(data.frameData.coordinates);
       iv.loadAtomAttributes(data.frameData.secondaryStructure);
       iv.rebuildScene(options);
       iv.zoomInto(options);
       iv.render();
   });
}