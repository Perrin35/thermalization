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
rz(1.323913) q[1];
sx q[1];
rz(-1.3988928) q[1];
sx q[1];
rz(-0.89371347) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0269256) q[0];
sx q[0];
rz(-2.6485778) q[0];
sx q[0];
rz(0.42401029) q[0];
rz(-pi) q[1];
rz(1.6581832) q[2];
sx q[2];
rz(-1.7636711) q[2];
sx q[2];
rz(-0.88763104) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1360769) q[1];
sx q[1];
rz(-1.7598333) q[1];
sx q[1];
rz(0.4398911) q[1];
rz(2.9017994) q[3];
sx q[3];
rz(-1.6231281) q[3];
sx q[3];
rz(-2.1870942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.90014234) q[2];
sx q[2];
rz(-0.19108471) q[2];
sx q[2];
rz(3.0639263) q[2];
rz(3.0302454) q[3];
sx q[3];
rz(-1.0166549) q[3];
sx q[3];
rz(-0.98285037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10649189) q[0];
sx q[0];
rz(-0.74420559) q[0];
sx q[0];
rz(2.8700854) q[0];
rz(2.5356281) q[1];
sx q[1];
rz(-2.7022336) q[1];
sx q[1];
rz(0.6449759) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5640372) q[0];
sx q[0];
rz(-2.5112069) q[0];
sx q[0];
rz(-2.1255998) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9276191) q[2];
sx q[2];
rz(-0.67751086) q[2];
sx q[2];
rz(-3.0626631) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6771705) q[1];
sx q[1];
rz(-1.7126516) q[1];
sx q[1];
rz(-0.22428288) q[1];
rz(-pi) q[2];
rz(-2.2521001) q[3];
sx q[3];
rz(-2.9078662) q[3];
sx q[3];
rz(-2.1672578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.490654) q[2];
sx q[2];
rz(-1.1488612) q[2];
sx q[2];
rz(0.33565721) q[2];
rz(0.13441864) q[3];
sx q[3];
rz(-2.4498144) q[3];
sx q[3];
rz(-2.5366096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8721039) q[0];
sx q[0];
rz(-1.1815) q[0];
sx q[0];
rz(2.9961725) q[0];
rz(-2.2272286) q[1];
sx q[1];
rz(-1.4375571) q[1];
sx q[1];
rz(-2.5252555) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59728226) q[0];
sx q[0];
rz(-2.4203886) q[0];
sx q[0];
rz(2.3116259) q[0];
x q[1];
rz(-0.049268289) q[2];
sx q[2];
rz(-2.081909) q[2];
sx q[2];
rz(-0.12898239) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3369477) q[1];
sx q[1];
rz(-1.27158) q[1];
sx q[1];
rz(-2.772259) q[1];
rz(-0.31026526) q[3];
sx q[3];
rz(-0.65830001) q[3];
sx q[3];
rz(-0.24228111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1472037) q[2];
sx q[2];
rz(-2.2115967) q[2];
sx q[2];
rz(-2.8623016) q[2];
rz(-2.764737) q[3];
sx q[3];
rz(-2.2241204) q[3];
sx q[3];
rz(0.74523029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84173161) q[0];
sx q[0];
rz(-2.0722516) q[0];
sx q[0];
rz(1.3413062) q[0];
rz(0.74814859) q[1];
sx q[1];
rz(-0.52296269) q[1];
sx q[1];
rz(2.0789007) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13650249) q[0];
sx q[0];
rz(-1.0074638) q[0];
sx q[0];
rz(0.13701771) q[0];
rz(-0.31514374) q[2];
sx q[2];
rz(-0.82568554) q[2];
sx q[2];
rz(-2.7100503) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6036883) q[1];
sx q[1];
rz(-0.66715566) q[1];
sx q[1];
rz(2.0869419) q[1];
rz(-2.483401) q[3];
sx q[3];
rz(-2.255283) q[3];
sx q[3];
rz(1.6940862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.62715379) q[2];
sx q[2];
rz(-1.9095162) q[2];
sx q[2];
rz(2.3717234) q[2];
rz(1.7047966) q[3];
sx q[3];
rz(-2.1261647) q[3];
sx q[3];
rz(-2.3427486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1771667) q[0];
sx q[0];
rz(-1.3097958) q[0];
sx q[0];
rz(2.8231296) q[0];
rz(0.98694363) q[1];
sx q[1];
rz(-1.3401597) q[1];
sx q[1];
rz(-0.14652227) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7413899) q[0];
sx q[0];
rz(-2.4331749) q[0];
sx q[0];
rz(1.0111459) q[0];
rz(0.96168965) q[2];
sx q[2];
rz(-0.95476788) q[2];
sx q[2];
rz(-2.2005657) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6121813) q[1];
sx q[1];
rz(-1.9093462) q[1];
sx q[1];
rz(1.3271133) q[1];
x q[2];
rz(0.39590368) q[3];
sx q[3];
rz(-0.77953458) q[3];
sx q[3];
rz(-2.5427839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5557308) q[2];
sx q[2];
rz(-0.18438688) q[2];
sx q[2];
rz(1.2510074) q[2];
rz(-2.3986554) q[3];
sx q[3];
rz(-1.6299959) q[3];
sx q[3];
rz(-1.9353297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0400405) q[0];
sx q[0];
rz(-2.0843625) q[0];
sx q[0];
rz(-1.936116) q[0];
rz(-2.7936753) q[1];
sx q[1];
rz(-1.1707062) q[1];
sx q[1];
rz(-1.233137) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41878149) q[0];
sx q[0];
rz(-1.486088) q[0];
sx q[0];
rz(-2.8891488) q[0];
x q[1];
rz(-2.6377524) q[2];
sx q[2];
rz(-0.8692534) q[2];
sx q[2];
rz(-1.9761249) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.72224381) q[1];
sx q[1];
rz(-1.0245242) q[1];
sx q[1];
rz(-0.52358475) q[1];
rz(-pi) q[2];
rz(0.55854256) q[3];
sx q[3];
rz(-0.59584801) q[3];
sx q[3];
rz(2.5003025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.39522383) q[2];
sx q[2];
rz(-0.93337983) q[2];
sx q[2];
rz(0.4371492) q[2];
rz(-0.44132597) q[3];
sx q[3];
rz(-2.2289942) q[3];
sx q[3];
rz(-0.92379409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8298892) q[0];
sx q[0];
rz(-1.7376816) q[0];
sx q[0];
rz(-2.8436858) q[0];
rz(-2.9333313) q[1];
sx q[1];
rz(-0.88349897) q[1];
sx q[1];
rz(-0.40685245) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7283467) q[0];
sx q[0];
rz(-1.5399031) q[0];
sx q[0];
rz(-1.4824162) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88201875) q[2];
sx q[2];
rz(-1.3830966) q[2];
sx q[2];
rz(1.3903914) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1081097) q[1];
sx q[1];
rz(-0.86520608) q[1];
sx q[1];
rz(-0.34107855) q[1];
x q[2];
rz(0.49439821) q[3];
sx q[3];
rz(-1.0596078) q[3];
sx q[3];
rz(2.7123812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.75140658) q[2];
sx q[2];
rz(-2.4441256) q[2];
sx q[2];
rz(0.81525272) q[2];
rz(-2.8149572) q[3];
sx q[3];
rz(-1.1465466) q[3];
sx q[3];
rz(-3.1281285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99454749) q[0];
sx q[0];
rz(-0.70962405) q[0];
sx q[0];
rz(0.56353322) q[0];
rz(-1.7067478) q[1];
sx q[1];
rz(-2.1610465) q[1];
sx q[1];
rz(-0.30418667) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4124111) q[0];
sx q[0];
rz(-0.92386073) q[0];
sx q[0];
rz(-1.5280194) q[0];
rz(-2.2772592) q[2];
sx q[2];
rz(-2.0385409) q[2];
sx q[2];
rz(-2.7815429) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1796602) q[1];
sx q[1];
rz(-1.0294401) q[1];
sx q[1];
rz(2.8983745) q[1];
x q[2];
rz(-1.0902152) q[3];
sx q[3];
rz(-2.6811185) q[3];
sx q[3];
rz(-2.0802896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4284511) q[2];
sx q[2];
rz(-2.9623172) q[2];
sx q[2];
rz(-0.98503867) q[2];
rz(1.8200412) q[3];
sx q[3];
rz(-1.225084) q[3];
sx q[3];
rz(-2.9086746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4317076) q[0];
sx q[0];
rz(-0.49254492) q[0];
sx q[0];
rz(-0.4044958) q[0];
rz(0.74199039) q[1];
sx q[1];
rz(-1.208035) q[1];
sx q[1];
rz(2.580339) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68414664) q[0];
sx q[0];
rz(-1.976891) q[0];
sx q[0];
rz(-2.2133064) q[0];
x q[1];
rz(2.2225631) q[2];
sx q[2];
rz(-0.56927753) q[2];
sx q[2];
rz(0.33212663) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.83744637) q[1];
sx q[1];
rz(-1.8032025) q[1];
sx q[1];
rz(1.2088319) q[1];
x q[2];
rz(1.1981549) q[3];
sx q[3];
rz(-0.85470457) q[3];
sx q[3];
rz(1.3195697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3227417) q[2];
sx q[2];
rz(-0.35917869) q[2];
sx q[2];
rz(-0.90510577) q[2];
rz(0.94308606) q[3];
sx q[3];
rz(-1.9102996) q[3];
sx q[3];
rz(0.73956195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1484225) q[0];
sx q[0];
rz(-0.44023308) q[0];
sx q[0];
rz(-2.7711208) q[0];
rz(-0.91520339) q[1];
sx q[1];
rz(-1.4280495) q[1];
sx q[1];
rz(-0.76348335) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.59635) q[0];
sx q[0];
rz(-0.93726369) q[0];
sx q[0];
rz(-1.088953) q[0];
x q[1];
rz(-2.6652226) q[2];
sx q[2];
rz(-0.61846501) q[2];
sx q[2];
rz(-0.16099421) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50824805) q[1];
sx q[1];
rz(-0.92735433) q[1];
sx q[1];
rz(0.28336759) q[1];
rz(-pi) q[2];
rz(-0.053289847) q[3];
sx q[3];
rz(-1.6052979) q[3];
sx q[3];
rz(0.49279398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8107599) q[2];
sx q[2];
rz(-1.5786889) q[2];
sx q[2];
rz(1.8806489) q[2];
rz(1.8597982) q[3];
sx q[3];
rz(-1.0303048) q[3];
sx q[3];
rz(1.9542046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8967165) q[0];
sx q[0];
rz(-1.5302932) q[0];
sx q[0];
rz(1.2320919) q[0];
rz(2.2121519) q[1];
sx q[1];
rz(-2.1771931) q[1];
sx q[1];
rz(-2.8246115) q[1];
rz(-1.2435883) q[2];
sx q[2];
rz(-1.6251593) q[2];
sx q[2];
rz(-1.9030824) q[2];
rz(-0.55616641) q[3];
sx q[3];
rz(-0.63197613) q[3];
sx q[3];
rz(1.9356884) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
