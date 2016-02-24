var kriging = require('./kriging');

describe('kriging', function(){
    it('makes a basic prediction', function(done){

        var x = [];
        var y = [];
        var t = [];

        // create a fake dataset where t is 100 when lat is east
        for (var i = 0; i < 100; i++){
            x[i] = (-180)+Math.random()*360;
            y[i] = (-90)+Math.random()*180;
            t[i] = (x[i] > 0) ? 100 : 0;
        }

        var variogram = kriging.train(t,x,y,"exponential", 0, 10);

        if (!variogram) return done("variogram is null");

        if (kriging.predict(180, 0, variogram) < 50) done("unexpected result (<50)");
        if (kriging.predict(-180, 0, variogram) > 50) done("unexpected result (>50)");

        done();
    });
});