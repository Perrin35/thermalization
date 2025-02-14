OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0316226) q[0];
sx q[0];
rz(-0.87183824) q[0];
sx q[0];
rz(-0.81575552) q[0];
rz(-1.8176796) q[1];
sx q[1];
rz(-1.7426999) q[1];
sx q[1];
rz(0.89371347) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92233673) q[0];
sx q[0];
rz(-1.3748264) q[0];
sx q[0];
rz(0.45536572) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42025153) q[2];
sx q[2];
rz(-0.21152341) q[2];
sx q[2];
rz(-0.45892059) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0055157) q[1];
sx q[1];
rz(-1.3817594) q[1];
sx q[1];
rz(2.7017016) q[1];
x q[2];
rz(-2.9245236) q[3];
sx q[3];
rz(-2.8962629) q[3];
sx q[3];
rz(-2.7360578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2414503) q[2];
sx q[2];
rz(-2.9505079) q[2];
sx q[2];
rz(0.07766635) q[2];
rz(-3.0302454) q[3];
sx q[3];
rz(-2.1249378) q[3];
sx q[3];
rz(-0.98285037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0351008) q[0];
sx q[0];
rz(-2.3973871) q[0];
sx q[0];
rz(-0.27150723) q[0];
rz(-2.5356281) q[1];
sx q[1];
rz(-2.7022336) q[1];
sx q[1];
rz(2.4966168) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5640372) q[0];
sx q[0];
rz(-0.63038578) q[0];
sx q[0];
rz(-2.1255998) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66627247) q[2];
sx q[2];
rz(-1.4372908) q[2];
sx q[2];
rz(-1.8174416) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6771705) q[1];
sx q[1];
rz(-1.7126516) q[1];
sx q[1];
rz(2.9173098) q[1];
rz(-2.2521001) q[3];
sx q[3];
rz(-0.23372641) q[3];
sx q[3];
rz(-0.97433486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6509387) q[2];
sx q[2];
rz(-1.1488612) q[2];
sx q[2];
rz(0.33565721) q[2];
rz(-0.13441864) q[3];
sx q[3];
rz(-0.69177827) q[3];
sx q[3];
rz(0.604983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26948872) q[0];
sx q[0];
rz(-1.1815) q[0];
sx q[0];
rz(0.14542018) q[0];
rz(-2.2272286) q[1];
sx q[1];
rz(-1.7040355) q[1];
sx q[1];
rz(2.5252555) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4805883) q[0];
sx q[0];
rz(-1.0618774) q[0];
sx q[0];
rz(0.53553217) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0923244) q[2];
sx q[2];
rz(-1.0596837) q[2];
sx q[2];
rz(0.12898239) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7246252) q[1];
sx q[1];
rz(-2.6706302) q[1];
sx q[1];
rz(-0.70711406) q[1];
rz(-pi) q[2];
rz(1.3389204) q[3];
sx q[3];
rz(-2.1926741) q[3];
sx q[3];
rz(0.14280126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1472037) q[2];
sx q[2];
rz(-0.92999593) q[2];
sx q[2];
rz(0.27929107) q[2];
rz(-2.764737) q[3];
sx q[3];
rz(-2.2241204) q[3];
sx q[3];
rz(-2.3963624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84173161) q[0];
sx q[0];
rz(-2.0722516) q[0];
sx q[0];
rz(1.3413062) q[0];
rz(2.3934441) q[1];
sx q[1];
rz(-0.52296269) q[1];
sx q[1];
rz(-2.0789007) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5077909) q[0];
sx q[0];
rz(-1.6865382) q[0];
sx q[0];
rz(-2.1383889) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31514374) q[2];
sx q[2];
rz(-0.82568554) q[2];
sx q[2];
rz(0.43154237) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6036883) q[1];
sx q[1];
rz(-0.66715566) q[1];
sx q[1];
rz(-1.0546507) q[1];
rz(0.76983863) q[3];
sx q[3];
rz(-2.0645118) q[3];
sx q[3];
rz(-2.5635886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5144389) q[2];
sx q[2];
rz(-1.9095162) q[2];
sx q[2];
rz(2.3717234) q[2];
rz(1.7047966) q[3];
sx q[3];
rz(-1.0154279) q[3];
sx q[3];
rz(2.3427486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.964426) q[0];
sx q[0];
rz(-1.3097958) q[0];
sx q[0];
rz(0.31846309) q[0];
rz(0.98694363) q[1];
sx q[1];
rz(-1.3401597) q[1];
sx q[1];
rz(2.9950704) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28956911) q[0];
sx q[0];
rz(-0.98678723) q[0];
sx q[0];
rz(-2.7147074) q[0];
rz(-pi) q[1];
rz(-0.96168965) q[2];
sx q[2];
rz(-0.95476788) q[2];
sx q[2];
rz(2.2005657) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1825936) q[1];
sx q[1];
rz(-1.800391) q[1];
sx q[1];
rz(2.793538) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9349443) q[3];
sx q[3];
rz(-0.86508646) q[3];
sx q[3];
rz(-3.074086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5557308) q[2];
sx q[2];
rz(-2.9572058) q[2];
sx q[2];
rz(1.2510074) q[2];
rz(0.74293724) q[3];
sx q[3];
rz(-1.6299959) q[3];
sx q[3];
rz(-1.9353297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10155216) q[0];
sx q[0];
rz(-2.0843625) q[0];
sx q[0];
rz(1.936116) q[0];
rz(0.34791738) q[1];
sx q[1];
rz(-1.1707062) q[1];
sx q[1];
rz(-1.233137) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3065227) q[0];
sx q[0];
rz(-0.26598922) q[0];
sx q[0];
rz(2.8138922) q[0];
x q[1];
rz(-0.50384028) q[2];
sx q[2];
rz(-2.2723393) q[2];
sx q[2];
rz(-1.9761249) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.55716276) q[1];
sx q[1];
rz(-1.129377) q[1];
sx q[1];
rz(-0.95870383) q[1];
rz(-pi) q[2];
rz(-1.2258271) q[3];
sx q[3];
rz(-2.0668092) q[3];
sx q[3];
rz(1.8535875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.39522383) q[2];
sx q[2];
rz(-0.93337983) q[2];
sx q[2];
rz(-2.7044435) q[2];
rz(2.7002667) q[3];
sx q[3];
rz(-0.91259846) q[3];
sx q[3];
rz(0.92379409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8298892) q[0];
sx q[0];
rz(-1.403911) q[0];
sx q[0];
rz(-0.29790685) q[0];
rz(-2.9333313) q[1];
sx q[1];
rz(-0.88349897) q[1];
sx q[1];
rz(2.7347402) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9637311) q[0];
sx q[0];
rz(-3.047982) q[0];
sx q[0];
rz(-1.2340182) q[0];
x q[1];
rz(-0.88201875) q[2];
sx q[2];
rz(-1.758496) q[2];
sx q[2];
rz(-1.3903914) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6089203) q[1];
sx q[1];
rz(-2.3708276) q[1];
sx q[1];
rz(1.9449598) q[1];
rz(-0.49439821) q[3];
sx q[3];
rz(-2.0819849) q[3];
sx q[3];
rz(-0.4292115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3901861) q[2];
sx q[2];
rz(-2.4441256) q[2];
sx q[2];
rz(-2.3263399) q[2];
rz(0.32663545) q[3];
sx q[3];
rz(-1.1465466) q[3];
sx q[3];
rz(0.013464125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1470452) q[0];
sx q[0];
rz(-0.70962405) q[0];
sx q[0];
rz(-0.56353322) q[0];
rz(-1.7067478) q[1];
sx q[1];
rz(-0.98054612) q[1];
sx q[1];
rz(0.30418667) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2741843) q[0];
sx q[0];
rz(-1.6049258) q[0];
sx q[0];
rz(0.64737583) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5553914) q[2];
sx q[2];
rz(-2.1888141) q[2];
sx q[2];
rz(1.5635334) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1796602) q[1];
sx q[1];
rz(-1.0294401) q[1];
sx q[1];
rz(-0.24321817) q[1];
x q[2];
rz(-1.1564141) q[3];
sx q[3];
rz(-1.3638921) q[3];
sx q[3];
rz(0.072553886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7131416) q[2];
sx q[2];
rz(-2.9623172) q[2];
sx q[2];
rz(0.98503867) q[2];
rz(-1.8200412) q[3];
sx q[3];
rz(-1.9165087) q[3];
sx q[3];
rz(0.2329181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.709885) q[0];
sx q[0];
rz(-0.49254492) q[0];
sx q[0];
rz(2.7370969) q[0];
rz(2.3996023) q[1];
sx q[1];
rz(-1.9335577) q[1];
sx q[1];
rz(-0.56125364) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7405734) q[0];
sx q[0];
rz(-0.74441806) q[0];
sx q[0];
rz(-0.94835001) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91902952) q[2];
sx q[2];
rz(-0.56927753) q[2];
sx q[2];
rz(2.809466) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8617647) q[1];
sx q[1];
rz(-2.7142254) q[1];
sx q[1];
rz(-0.98160569) q[1];
rz(-pi) q[2];
rz(0.75144479) q[3];
sx q[3];
rz(-1.8490233) q[3];
sx q[3];
rz(-3.1415526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3227417) q[2];
sx q[2];
rz(-0.35917869) q[2];
sx q[2];
rz(0.90510577) q[2];
rz(-0.94308606) q[3];
sx q[3];
rz(-1.9102996) q[3];
sx q[3];
rz(2.4020307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1484225) q[0];
sx q[0];
rz(-0.44023308) q[0];
sx q[0];
rz(-2.7711208) q[0];
rz(2.2263893) q[1];
sx q[1];
rz(-1.7135432) q[1];
sx q[1];
rz(0.76348335) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8158096) q[0];
sx q[0];
rz(-1.9535583) q[0];
sx q[0];
rz(2.4494657) q[0];
rz(-0.47637003) q[2];
sx q[2];
rz(-0.61846501) q[2];
sx q[2];
rz(-2.9805984) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1814628) q[1];
sx q[1];
rz(-2.4467111) q[1];
sx q[1];
rz(-1.927666) q[1];
rz(-2.5666272) q[3];
sx q[3];
rz(-0.063474726) q[3];
sx q[3];
rz(-0.50395645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3308328) q[2];
sx q[2];
rz(-1.5629038) q[2];
sx q[2];
rz(1.8806489) q[2];
rz(1.8597982) q[3];
sx q[3];
rz(-2.1112879) q[3];
sx q[3];
rz(-1.9542046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24487615) q[0];
sx q[0];
rz(-1.5302932) q[0];
sx q[0];
rz(1.2320919) q[0];
rz(2.2121519) q[1];
sx q[1];
rz(-2.1771931) q[1];
sx q[1];
rz(-2.8246115) q[1];
rz(-1.4030761) q[2];
sx q[2];
rz(-2.8100596) q[2];
sx q[2];
rz(2.6505627) q[2];
rz(-1.2019602) q[3];
sx q[3];
rz(-1.0452253) q[3];
sx q[3];
rz(-0.54946231) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
