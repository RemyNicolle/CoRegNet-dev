
var cyOutputBinding = new Shiny.OutputBinding();

  $.extend(cyOutputBinding, {
    find: function(scope) {
      return $(scope).find('.cynetwork');
    },
    renderValue: function(el, data) {
      console.log(data.type);
    if(data.type=="network"){
      createNormalNet(data);
    }else if(data.type=="pienetwork"){                
      createPieNet(data);        
    }else if(data.type=="layout"){
      $('#cy').cytoscape(function(){
      var cy = this; // now we have a reference to the core
      cy.layout(  {
                       name: 'cose',
                       refresh:10,
                       numIter:100,
                       nodeRepulsion       : 500000,
                       padding:0.1
                       });
      });
    }else if(data.type=="update"){      
      $('#cy').cytoscape(function(){
        var cy = this; // now we have a reference to the core
        cy.batchData( data.names);
      });
    }else if(data.type=="select"){
      console.log(data.tosel);
      $('#cy').cytoscape(function(){
        var cy = this; 
        console.log("zoom to");
        console.log(data.tosel)
        var ele = cy.getElementById(data.tosel);
        cy.zoom({level: 1.5});
        cy.center(ele);
      });
    }else if(data.type=="test"){
      console.log(data.names);
      $('#cy').cytoscape(function(){
      var cy = this; 
      //  cy.batchData( data.names);
      });
    }else if(data.type =="url"){
     var mydiv = document.getElementById("myDiv");
    var aTag = document.createElement('a');
    aTag.setAttribute('href',"yourlink.htm");
    aTag.innerHTML = "link text";
    mydiv.appendChild(aTag);
    }else if(data.type =="png"){
      $('#cy').cytoscape(function(){
        var cy = this;                 
        
        document.getElementById("snap").setAttribute('src',cy.png({full:true}));
        
        document.getElementById("snap").setAttribute('src',cy.png({full:true}));
      });    
    }
  }
});


jQuery(document).ready(function($) {
    $('#legendNetwork_slick').dcSlick({
            location: 'top',
            align: 'right',
            offset: '0px',
            speed: 'fast',
            tabText: 'Network Legend',
            autoClose: true
    });
});


jQuery(document).ready(function($) {
    $('#legend_slick').dcSlick({
            location: 'top',
            align: 'right',
            offset: '130px',
            speed: 'fast',
            tabText: 'Legend',
            autoClose: true
    });
});


$(window).resize(function() {
    $('#cy').cytoscape(function(){
        var cy = this;
//        cy.fit();
        cy.center();
//        cy.resize();
    });
});


Shiny.outputBindings.register(cyOutputBinding);







var createNormalNet = function(data){


    var myNodes = [];
    var myEdges = [];

     for (var i = 0; i < data.nodes.length; i++) {
      myNodes.push({ data: data.nodes[i] });}      
 
    for (var i=0;  i < data.links.id.length; i++) {

        myEdges.push({
                       data: {
                       id:  data.links.id[i],
                       source: data.links.Reg1[i],
                       target: data.links.Reg2[i],
                       weight: data.links.Support[i],
                       color: data.links.color[i],
                       arrow :data.links.arrow[i]
                       }
                       });
    }

    
    $('#cy').cytoscape({
                       minZoom: 0.05,
                       maxZoom: 1.5,
                       layout: {
                       name: 'cose',
                       numIter:50,
                       nodeRepulsion       : 500000,
                       padding:0.1
                       },                       
                       style: cytoscape.stylesheet()
                       .selector('node')
                       .css({
                            'shape': 'ellipse',
                            'height': 'mapData(size, 0, 1, 20, 50)',
                            'width': 'mapData(size, 0, 1, 20, 50)',
                            'color':'white',
                            'content': 'data(name)',
                            'text-valign': 'center',
                            'text-outline-width': 2,
                            'text-outline-color': 'mapData(color, -10, 10, blue, red)',
                            'background-color': 'mapData(color, -10, 10, blue, red)'
                            })
                       .selector(':selected')
                       .css({
                            'border-width': 4,
                            'border-color': '#ccc'
                            })
                       .selector('edge')
                       .css({
                           'width': 'mapData(weight, 0.001, 1, 1, 5)',
//                            'width': '2',
                            'target-arrow-shape': 'data(arrow)',
                            'line-color': 'data(color)',
                            'source-arrow-color': 'data(color)',
                            'target-arrow-color': 'data(color)'
                            })
                       .selector('.faded')
                       .css({
                            'opacity': 0.5,
                            'text-opacity': 0
                            }),
                       elements: {
                       nodes: myNodes,
                       edges: myEdges,
                       },
                       ready: function(){
                       window.cy = this;
                       cy.on('tap', 'node', function(e){
                             var node = e.cyTarget;
                             var neighborhood = node.neighborhood().add(node);

                            Shiny.onInputChange("cy", node.data("name"));

                             cy.elements().addClass('faded');
                             neighborhood.removeClass('faded');
                             });
                cy.on('select', 'node', function(event){
                     var selnodes = [];
                     cy.$(':selected').each(function(i,ele){
                                        selnodes.push({
                                                name: ele.data("name") });
                                            });                                            
                             console.log( cy.$(':selected').length);
                     Shiny.onInputChange("cy",selnodes);
                });
                       cy.on('tap', function(e){
                             if( e.cyTarget === cy ){
                             cy.elements().removeClass('faded');
                              Shiny.onInputChange("cy", "NULL");
                             }
                          
                             });
                       }
                       });
   }





var createPieNet = function(data){
    var myNodes = [];
    var myEdges = [];
 
      for (var i = 0; i < data.nodes.length; i++) {
      myNodes.push({ data: data.nodes[i] });}
    for (var i=0;  i < data.links.id.length; i++) {
        myEdges.push({
                       data: {
                       id:  data.links.id[i],
                       source: data.links.Reg1[i],
                       target: data.links.Reg2[i],
                       weight: data.links.Support[i],
                       color: data.links.color[i],
                       arrow :data.links.arrow[i]
                       }
                       });
    }
    
    $('#cy').cytoscape({
                       minZoom: 0.05,
                       maxZoom: 1.5,
                       hideEdgesOnViewport: true,
                       hideLabelsOnViewport: true,
                       layout: {
                       name: 'cose',
                       numIter:50,
                       nodeRepulsion       : 500000,
                       padding:0.1
                       },                       
                       style: cytoscape.stylesheet()
                       .selector('node')
                       .css({
                            'shape': 'ellipse',
                            'height': 'mapData(size, 0, 1, 20, 50)',
                            'width': 'mapData(size, 0, 1, 20, 50)',
                            'color':'white',
                            'content': 'data(name)',
                            'text-outline-width': 2,
                            'text-outline-color': 'mapData(color, -10, 10, blue, red)',
                            'background-color': 'mapData(color, -10, 10, blue, red)',
                            'pie-size': '80%',
                            'pie-1-background-color': '#2b83ba',
                            'pie-1-background-size': 'mapData(delet, 0, 100, 0, 100)',
                            'pie-2-background-color': '#abdda4',
                            'pie-2-background-size': 'mapData(loss, 0, 100, 0, 100)',
                            'pie-3-background-color': '#ffffbf',
                            'pie-3-background-size': 'mapData(diplo, 0, 100, 0, 100)',
                            'pie-4-background-color': '#fdae61',
                            'pie-4-background-size': 'mapData(gain, 0, 100, 0, 100)',
                            'pie-5-background-color': '#d7191c',
                            'pie-5-background-size': 'mapData(ampli, 0, 100, 0, 100)',
                            'pie-6-background-color': '#888888',
                            'pie-6-background-size': 'mapData(unknown, 0, 100, 0, 100)'
                            })
                       .selector(':selected')
                       .css({
                            'border-width': 4,
                            'border-color': '#ccc'
                            })
                       .selector('edge')
                       .css({
                           'width': 'mapData(weight, 0.001, 1, 1, 5)',
                            'target-arrow-shape': 'data(arrow)',
                            'line-color': 'data(color)',
                            'source-arrow-color': 'data(color)',
                            'target-arrow-color': 'data(color)'
                            })
                       .selector('.faded')
                       .css({
                            'opacity': 0.5,
                            'text-opacity': 0
                            }),
                       
                       elements: {
                       nodes: myNodes,
                       edges: myEdges,
                       },
                       ready: function(){
                       window.cy = this;
            
                       
                       
                       cy.on('tap', 'node', function(e){
                             var node = e.cyTarget;
                             var neighborhood = node.neighborhood().add(node);

                            Shiny.onInputChange("cy", node.data("name"));

                             cy.elements().addClass('faded');
                             neighborhood.removeClass('faded');
                             });
                
                
                cy.on('select', 'node', function(event){
                     var selnodes = [];
                     cy.$(':selected').each(function(i,ele){
                                        selnodes.push({
                                                name: ele.data("name") });
                                            });                                            
                             console.log( cy.$(':selected').length);
                     Shiny.onInputChange("cy",selnodes);
                });

                       cy.on('tap', function(e){
                             if( e.cyTarget === cy ){
                             cy.elements().removeClass('faded');
                              Shiny.onInputChange("cy", "NULL");
                             }
                          
                             });
                       }
                       });
   }






