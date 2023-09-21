OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7005641) q[0];
sx q[0];
rz(1.1428042) q[0];
sx q[0];
rz(11.354843) q[0];
rz(2.9149574) q[1];
sx q[1];
rz(-1.5645138) q[1];
sx q[1];
rz(-0.29830631) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0878108) q[0];
sx q[0];
rz(-2.2349173) q[0];
sx q[0];
rz(1.640663) q[0];
rz(1.9036129) q[2];
sx q[2];
rz(-1.6593854) q[2];
sx q[2];
rz(0.36117902) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5391985) q[1];
sx q[1];
rz(-2.1089923) q[1];
sx q[1];
rz(2.5799275) q[1];
rz(-pi) q[2];
rz(-0.88793036) q[3];
sx q[3];
rz(-0.13266064) q[3];
sx q[3];
rz(2.5761029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4937218) q[2];
sx q[2];
rz(-1.2080668) q[2];
sx q[2];
rz(0.99386627) q[2];
rz(0.99938756) q[3];
sx q[3];
rz(-1.9013654) q[3];
sx q[3];
rz(-0.58888155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4988929) q[0];
sx q[0];
rz(-1.9357341) q[0];
sx q[0];
rz(0.018218536) q[0];
rz(2.3253564) q[1];
sx q[1];
rz(-1.0304334) q[1];
sx q[1];
rz(0.47168628) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52662151) q[0];
sx q[0];
rz(-1.3437628) q[0];
sx q[0];
rz(1.4814266) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6623146) q[2];
sx q[2];
rz(-2.5844378) q[2];
sx q[2];
rz(-0.011475871) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8397123) q[1];
sx q[1];
rz(-2.1228959) q[1];
sx q[1];
rz(1.8865955) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3470115) q[3];
sx q[3];
rz(-0.61962485) q[3];
sx q[3];
rz(-2.1686045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.54962426) q[2];
sx q[2];
rz(-1.1529808) q[2];
sx q[2];
rz(-0.96898752) q[2];
rz(0.5747059) q[3];
sx q[3];
rz(-2.59022) q[3];
sx q[3];
rz(2.1000752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7085003) q[0];
sx q[0];
rz(-1.0555462) q[0];
sx q[0];
rz(2.95978) q[0];
rz(-2.0388942) q[1];
sx q[1];
rz(-1.6405374) q[1];
sx q[1];
rz(1.4556494) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11765471) q[0];
sx q[0];
rz(-1.0929937) q[0];
sx q[0];
rz(2.2295203) q[0];
rz(-1.1481029) q[2];
sx q[2];
rz(-0.87562497) q[2];
sx q[2];
rz(-2.9425651) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.94565832) q[1];
sx q[1];
rz(-3.0155026) q[1];
sx q[1];
rz(0.083421589) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5163243) q[3];
sx q[3];
rz(-1.4843974) q[3];
sx q[3];
rz(0.82739917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87749798) q[2];
sx q[2];
rz(-1.654518) q[2];
sx q[2];
rz(0.96763119) q[2];
rz(-0.72757059) q[3];
sx q[3];
rz(-1.2604159) q[3];
sx q[3];
rz(-0.23770604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2847292) q[0];
sx q[0];
rz(-0.52607042) q[0];
sx q[0];
rz(-2.5033584) q[0];
rz(-1.1278641) q[1];
sx q[1];
rz(-2.3141839) q[1];
sx q[1];
rz(-1.2329873) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9905332) q[0];
sx q[0];
rz(-0.49976832) q[0];
sx q[0];
rz(2.997056) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5523283) q[2];
sx q[2];
rz(-2.0406084) q[2];
sx q[2];
rz(1.1772917) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7113263) q[1];
sx q[1];
rz(-1.3923936) q[1];
sx q[1];
rz(3.1049411) q[1];
x q[2];
rz(2.9144985) q[3];
sx q[3];
rz(-0.89082754) q[3];
sx q[3];
rz(1.029315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2234852) q[2];
sx q[2];
rz(-1.2780259) q[2];
sx q[2];
rz(-2.3275862) q[2];
rz(2.0984086) q[3];
sx q[3];
rz(-0.63101763) q[3];
sx q[3];
rz(-1.957318) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84213132) q[0];
sx q[0];
rz(-1.3548387) q[0];
sx q[0];
rz(2.2498851) q[0];
rz(1.2437598) q[1];
sx q[1];
rz(-1.3777106) q[1];
sx q[1];
rz(-0.2125425) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4310303) q[0];
sx q[0];
rz(-3.115603) q[0];
sx q[0];
rz(0.89528577) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6264621) q[2];
sx q[2];
rz(-2.5362483) q[2];
sx q[2];
rz(2.2948613) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5383496) q[1];
sx q[1];
rz(-2.1117359) q[1];
sx q[1];
rz(2.4451838) q[1];
x q[2];
rz(-0.84380031) q[3];
sx q[3];
rz(-1.2146597) q[3];
sx q[3];
rz(-0.86405495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1084958) q[2];
sx q[2];
rz(-0.97110811) q[2];
sx q[2];
rz(0.5212211) q[2];
rz(-1.3850348) q[3];
sx q[3];
rz(-1.9135467) q[3];
sx q[3];
rz(1.2683755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824771) q[0];
sx q[0];
rz(-0.92882597) q[0];
sx q[0];
rz(-0.22512063) q[0];
rz(1.7865932) q[1];
sx q[1];
rz(-1.0083818) q[1];
sx q[1];
rz(2.7640142) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7901944) q[0];
sx q[0];
rz(-1.5488008) q[0];
sx q[0];
rz(-1.8368594) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2874561) q[2];
sx q[2];
rz(-1.6171347) q[2];
sx q[2];
rz(0.031678274) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7403455) q[1];
sx q[1];
rz(-2.0013072) q[1];
sx q[1];
rz(0.11534782) q[1];
x q[2];
rz(-2.3466831) q[3];
sx q[3];
rz(-1.5494293) q[3];
sx q[3];
rz(1.4710609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.98465115) q[2];
sx q[2];
rz(-0.7876544) q[2];
sx q[2];
rz(2.6605576) q[2];
rz(-0.40361079) q[3];
sx q[3];
rz(-2.1026881) q[3];
sx q[3];
rz(0.31479442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0890546) q[0];
sx q[0];
rz(-0.60482329) q[0];
sx q[0];
rz(2.9470434) q[0];
rz(2.9220707) q[1];
sx q[1];
rz(-1.4621282) q[1];
sx q[1];
rz(-0.25442466) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62691488) q[0];
sx q[0];
rz(-1.8906381) q[0];
sx q[0];
rz(0.11492782) q[0];
rz(-pi) q[1];
rz(1.3527855) q[2];
sx q[2];
rz(-2.5180452) q[2];
sx q[2];
rz(-0.64507285) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.88528819) q[1];
sx q[1];
rz(-1.8418152) q[1];
sx q[1];
rz(-2.6726252) q[1];
rz(-pi) q[2];
rz(0.80612225) q[3];
sx q[3];
rz(-2.5297909) q[3];
sx q[3];
rz(0.34477371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1640132) q[2];
sx q[2];
rz(-1.3914725) q[2];
sx q[2];
rz(1.3158201) q[2];
rz(-0.87604648) q[3];
sx q[3];
rz(-0.13893572) q[3];
sx q[3];
rz(-2.1379437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9309689) q[0];
sx q[0];
rz(-2.7656778) q[0];
sx q[0];
rz(1.4550495) q[0];
rz(0.82398206) q[1];
sx q[1];
rz(-0.95183698) q[1];
sx q[1];
rz(-1.5751858) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63021916) q[0];
sx q[0];
rz(-0.89238088) q[0];
sx q[0];
rz(1.1648965) q[0];
x q[1];
rz(-2.1112061) q[2];
sx q[2];
rz(-1.4166797) q[2];
sx q[2];
rz(2.266778) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.52289256) q[1];
sx q[1];
rz(-2.61781) q[1];
sx q[1];
rz(0.79843847) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.992222) q[3];
sx q[3];
rz(-2.0391658) q[3];
sx q[3];
rz(0.8120265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2660797) q[2];
sx q[2];
rz(-1.9078887) q[2];
sx q[2];
rz(2.8640462) q[2];
rz(-1.6905486) q[3];
sx q[3];
rz(-0.45193672) q[3];
sx q[3];
rz(1.014876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16185109) q[0];
sx q[0];
rz(-0.45270544) q[0];
sx q[0];
rz(1.45654) q[0];
rz(-2.5121571) q[1];
sx q[1];
rz(-1.1673735) q[1];
sx q[1];
rz(1.1368407) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6900078) q[0];
sx q[0];
rz(-1.0613872) q[0];
sx q[0];
rz(0.33126979) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9929664) q[2];
sx q[2];
rz(-1.9151701) q[2];
sx q[2];
rz(-0.070377199) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6362308) q[1];
sx q[1];
rz(-2.8103235) q[1];
sx q[1];
rz(-1.2760217) q[1];
x q[2];
rz(2.4056899) q[3];
sx q[3];
rz(-1.8527485) q[3];
sx q[3];
rz(1.499093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4920766) q[2];
sx q[2];
rz(-1.4294383) q[2];
sx q[2];
rz(2.1006404) q[2];
rz(-3.1395636) q[3];
sx q[3];
rz(-0.40015951) q[3];
sx q[3];
rz(2.8295529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8390389) q[0];
sx q[0];
rz(-0.22452393) q[0];
sx q[0];
rz(-2.1955406) q[0];
rz(0.91167766) q[1];
sx q[1];
rz(-1.2152351) q[1];
sx q[1];
rz(2.5295703) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34110585) q[0];
sx q[0];
rz(-0.90011156) q[0];
sx q[0];
rz(0.97408803) q[0];
rz(1.3304747) q[2];
sx q[2];
rz(-1.9069306) q[2];
sx q[2];
rz(-1.0898332) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7358688) q[1];
sx q[1];
rz(-1.6940261) q[1];
sx q[1];
rz(0.65995364) q[1];
rz(1.9480115) q[3];
sx q[3];
rz(-1.0672788) q[3];
sx q[3];
rz(-2.6138888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0845906) q[2];
sx q[2];
rz(-0.64208639) q[2];
sx q[2];
rz(-0.65336147) q[2];
rz(0.35081321) q[3];
sx q[3];
rz(-1.5272798) q[3];
sx q[3];
rz(-2.4408834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54031298) q[0];
sx q[0];
rz(-1.5355587) q[0];
sx q[0];
rz(-3.0328947) q[0];
rz(0.75469771) q[1];
sx q[1];
rz(-1.8221868) q[1];
sx q[1];
rz(1.6356161) q[1];
rz(2.4621261) q[2];
sx q[2];
rz(-2.0308528) q[2];
sx q[2];
rz(0.46926342) q[2];
rz(0.084372088) q[3];
sx q[3];
rz(-1.2481239) q[3];
sx q[3];
rz(0.39831755) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
