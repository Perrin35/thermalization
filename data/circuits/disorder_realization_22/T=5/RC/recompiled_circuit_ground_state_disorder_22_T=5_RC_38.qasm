OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8393505) q[0];
sx q[0];
rz(-1.8678764) q[0];
sx q[0];
rz(-0.94925517) q[0];
rz(-1.0957837) q[1];
sx q[1];
rz(2.1646808) q[1];
sx q[1];
rz(10.77471) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95800864) q[0];
sx q[0];
rz(-1.6537979) q[0];
sx q[0];
rz(-0.50358332) q[0];
rz(0.33265661) q[2];
sx q[2];
rz(-1.7094618) q[2];
sx q[2];
rz(2.0227416) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9266859) q[1];
sx q[1];
rz(-1.5027355) q[1];
sx q[1];
rz(2.7665274) q[1];
rz(2.1175577) q[3];
sx q[3];
rz(-2.6461678) q[3];
sx q[3];
rz(-1.6273496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.533941) q[2];
sx q[2];
rz(-2.0484296) q[2];
sx q[2];
rz(2.9489813) q[2];
rz(1.8588148) q[3];
sx q[3];
rz(-1.1056113) q[3];
sx q[3];
rz(2.5530596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.823536) q[0];
sx q[0];
rz(-0.41724351) q[0];
sx q[0];
rz(0.70660025) q[0];
rz(-2.508714) q[1];
sx q[1];
rz(-2.0426079) q[1];
sx q[1];
rz(0.28876567) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5486886) q[0];
sx q[0];
rz(-0.0026772896) q[0];
sx q[0];
rz(-0.38627465) q[0];
rz(-pi) q[1];
rz(-1.7328506) q[2];
sx q[2];
rz(-1.7758435) q[2];
sx q[2];
rz(2.1265202) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.56882492) q[1];
sx q[1];
rz(-1.8563885) q[1];
sx q[1];
rz(1.7920124) q[1];
rz(1.1488647) q[3];
sx q[3];
rz(-0.85638085) q[3];
sx q[3];
rz(-1.8569291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4660945) q[2];
sx q[2];
rz(-1.0885295) q[2];
sx q[2];
rz(1.6509854) q[2];
rz(2.364482) q[3];
sx q[3];
rz(-1.4584352) q[3];
sx q[3];
rz(1.3597663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31558388) q[0];
sx q[0];
rz(-1.6827787) q[0];
sx q[0];
rz(-2.7050731) q[0];
rz(0.79477683) q[1];
sx q[1];
rz(-1.7710641) q[1];
sx q[1];
rz(2.2025542) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7815112) q[0];
sx q[0];
rz(-1.0699125) q[0];
sx q[0];
rz(-1.341218) q[0];
rz(-pi) q[1];
rz(-0.58231488) q[2];
sx q[2];
rz(-2.2852201) q[2];
sx q[2];
rz(2.0272875) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8104097) q[1];
sx q[1];
rz(-2.4723834) q[1];
sx q[1];
rz(-2.616908) q[1];
x q[2];
rz(1.9179929) q[3];
sx q[3];
rz(-2.5592862) q[3];
sx q[3];
rz(0.83707419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7006435) q[2];
sx q[2];
rz(-1.0627397) q[2];
sx q[2];
rz(-0.6353333) q[2];
rz(-2.3144531) q[3];
sx q[3];
rz(-0.45378903) q[3];
sx q[3];
rz(1.2686096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24469911) q[0];
sx q[0];
rz(-0.06299717) q[0];
sx q[0];
rz(-2.6196106) q[0];
rz(-0.43102795) q[1];
sx q[1];
rz(-1.6072175) q[1];
sx q[1];
rz(-0.65188754) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7303607) q[0];
sx q[0];
rz(-1.8069977) q[0];
sx q[0];
rz(-0.51933164) q[0];
rz(-pi) q[1];
rz(-0.84649936) q[2];
sx q[2];
rz(-1.9545467) q[2];
sx q[2];
rz(-2.4946314) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.83623278) q[1];
sx q[1];
rz(-1.9499602) q[1];
sx q[1];
rz(-0.01474488) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96792696) q[3];
sx q[3];
rz(-0.9612007) q[3];
sx q[3];
rz(2.2948007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.42644694) q[2];
sx q[2];
rz(-1.8603674) q[2];
sx q[2];
rz(2.35671) q[2];
rz(-2.6134885) q[3];
sx q[3];
rz(-1.2930861) q[3];
sx q[3];
rz(-1.8538792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89609471) q[0];
sx q[0];
rz(-0.69480768) q[0];
sx q[0];
rz(-2.3520663) q[0];
rz(-0.49742571) q[1];
sx q[1];
rz(-2.1010294) q[1];
sx q[1];
rz(-1.3015889) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4052491) q[0];
sx q[0];
rz(-2.1642125) q[0];
sx q[0];
rz(1.7878535) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.80255) q[2];
sx q[2];
rz(-1.6280988) q[2];
sx q[2];
rz(-0.49446854) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9250946) q[1];
sx q[1];
rz(-1.3984233) q[1];
sx q[1];
rz(-2.4045776) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9381125) q[3];
sx q[3];
rz(-1.5496033) q[3];
sx q[3];
rz(-0.68875865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.54399458) q[2];
sx q[2];
rz(-1.5388637) q[2];
sx q[2];
rz(-0.27628118) q[2];
rz(-0.9497408) q[3];
sx q[3];
rz(-0.26849982) q[3];
sx q[3];
rz(-0.55324078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62040579) q[0];
sx q[0];
rz(-0.036245417) q[0];
sx q[0];
rz(-0.94394839) q[0];
rz(-1.8432519) q[1];
sx q[1];
rz(-1.0669758) q[1];
sx q[1];
rz(-2.3640769) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80332887) q[0];
sx q[0];
rz(-1.216812) q[0];
sx q[0];
rz(0.5838809) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73410122) q[2];
sx q[2];
rz(-1.7673552) q[2];
sx q[2];
rz(0.71069709) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0447415) q[1];
sx q[1];
rz(-1.4440084) q[1];
sx q[1];
rz(-0.37146588) q[1];
rz(0.10482444) q[3];
sx q[3];
rz(-0.30571211) q[3];
sx q[3];
rz(1.0979528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5139318) q[2];
sx q[2];
rz(-1.7646503) q[2];
sx q[2];
rz(-2.7505006) q[2];
rz(-0.62134653) q[3];
sx q[3];
rz(-2.4070599) q[3];
sx q[3];
rz(-0.77596107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4710627) q[0];
sx q[0];
rz(-0.10469086) q[0];
sx q[0];
rz(-1.712557) q[0];
rz(0.96587005) q[1];
sx q[1];
rz(-1.8873676) q[1];
sx q[1];
rz(2.3557854) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87596506) q[0];
sx q[0];
rz(-0.64564359) q[0];
sx q[0];
rz(3.0558048) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9740392) q[2];
sx q[2];
rz(-2.2637562) q[2];
sx q[2];
rz(-0.6374109) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.475226) q[1];
sx q[1];
rz(-2.6738648) q[1];
sx q[1];
rz(-0.18233129) q[1];
rz(-pi) q[2];
rz(1.2811411) q[3];
sx q[3];
rz(-1.148734) q[3];
sx q[3];
rz(0.45988032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.69661951) q[2];
sx q[2];
rz(-2.5225621) q[2];
sx q[2];
rz(-2.6444198) q[2];
rz(-0.87388006) q[3];
sx q[3];
rz(-2.0495575) q[3];
sx q[3];
rz(1.2018275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9484321) q[0];
sx q[0];
rz(-0.30696294) q[0];
sx q[0];
rz(-0.69751414) q[0];
rz(2.7613617) q[1];
sx q[1];
rz(-2.0262599) q[1];
sx q[1];
rz(3.0013705) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2220737) q[0];
sx q[0];
rz(-1.4937048) q[0];
sx q[0];
rz(-0.18203966) q[0];
x q[1];
rz(-1.5319139) q[2];
sx q[2];
rz(-1.1302396) q[2];
sx q[2];
rz(-1.3300542) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6483874) q[1];
sx q[1];
rz(-0.81118656) q[1];
sx q[1];
rz(0.56583515) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7305364) q[3];
sx q[3];
rz(-1.7470164) q[3];
sx q[3];
rz(-0.4057708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2908638) q[2];
sx q[2];
rz(-1.2173434) q[2];
sx q[2];
rz(2.6503837) q[2];
rz(0.20719191) q[3];
sx q[3];
rz(-2.5184293) q[3];
sx q[3];
rz(-1.4108968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9687013) q[0];
sx q[0];
rz(-1.7270813) q[0];
sx q[0];
rz(-0.15039314) q[0];
rz(2.4405759) q[1];
sx q[1];
rz(-2.0471408) q[1];
sx q[1];
rz(1.6640123) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22035881) q[0];
sx q[0];
rz(-1.4495175) q[0];
sx q[0];
rz(-0.63684271) q[0];
rz(-pi) q[1];
rz(-2.5186164) q[2];
sx q[2];
rz(-2.9439363) q[2];
sx q[2];
rz(2.4923639) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8898082) q[1];
sx q[1];
rz(-1.4191333) q[1];
sx q[1];
rz(1.4633578) q[1];
x q[2];
rz(-1.2712237) q[3];
sx q[3];
rz(-0.75787395) q[3];
sx q[3];
rz(-1.8488499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.56665862) q[2];
sx q[2];
rz(-2.7329972) q[2];
sx q[2];
rz(-1.5236141) q[2];
rz(-0.59897113) q[3];
sx q[3];
rz(-2.6409769) q[3];
sx q[3];
rz(-0.18794255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69342518) q[0];
sx q[0];
rz(-2.7126815) q[0];
sx q[0];
rz(1.5018139) q[0];
rz(-2.2046454) q[1];
sx q[1];
rz(-2.1138771) q[1];
sx q[1];
rz(-2.1243336) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27148025) q[0];
sx q[0];
rz(-1.1944298) q[0];
sx q[0];
rz(0.41618213) q[0];
rz(-pi) q[1];
rz(1.4849365) q[2];
sx q[2];
rz(-1.6756264) q[2];
sx q[2];
rz(-1.3121999) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3268633) q[1];
sx q[1];
rz(-1.5376602) q[1];
sx q[1];
rz(-1.3087981) q[1];
rz(-pi) q[2];
x q[2];
rz(0.42933257) q[3];
sx q[3];
rz(-1.7547914) q[3];
sx q[3];
rz(1.3972549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.1169869) q[2];
sx q[2];
rz(-0.6822497) q[2];
sx q[2];
rz(1.5218706) q[2];
rz(1.4874602) q[3];
sx q[3];
rz(-1.4922851) q[3];
sx q[3];
rz(-1.3967995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7617154) q[0];
sx q[0];
rz(-1.3991671) q[0];
sx q[0];
rz(0.73200926) q[0];
rz(-1.6666182) q[1];
sx q[1];
rz(-1.9261618) q[1];
sx q[1];
rz(0.7934657) q[1];
rz(2.7421065) q[2];
sx q[2];
rz(-2.1151092) q[2];
sx q[2];
rz(-1.5110037) q[2];
rz(1.8745244) q[3];
sx q[3];
rz(-1.2645742) q[3];
sx q[3];
rz(-3.054061) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
