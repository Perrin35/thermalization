OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1144855) q[0];
sx q[0];
rz(-0.067582421) q[0];
sx q[0];
rz(-0.58900589) q[0];
rz(2.0677805) q[1];
sx q[1];
rz(-1.685073) q[1];
sx q[1];
rz(-2.9340802) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5440581) q[0];
sx q[0];
rz(-1.7701214) q[0];
sx q[0];
rz(1.8117649) q[0];
x q[1];
rz(-2.7314721) q[2];
sx q[2];
rz(-0.1995351) q[2];
sx q[2];
rz(-1.4360957) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.21440014) q[1];
sx q[1];
rz(-1.542558) q[1];
sx q[1];
rz(-0.0035889665) q[1];
x q[2];
rz(0.21353586) q[3];
sx q[3];
rz(-0.47225228) q[3];
sx q[3];
rz(2.6708598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8091858) q[2];
sx q[2];
rz(-3.1250592) q[2];
sx q[2];
rz(-0.22745505) q[2];
rz(-2.4849232) q[3];
sx q[3];
rz(-2.4469817) q[3];
sx q[3];
rz(-2.606126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.6856132) q[0];
sx q[0];
rz(-0.032129012) q[0];
sx q[0];
rz(0.44422126) q[0];
rz(1.1086858) q[1];
sx q[1];
rz(-1.8007092) q[1];
sx q[1];
rz(-1.1681555) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5801808) q[0];
sx q[0];
rz(-1.7636443) q[0];
sx q[0];
rz(2.4417402) q[0];
rz(-pi) q[1];
rz(1.6158478) q[2];
sx q[2];
rz(-1.117393) q[2];
sx q[2];
rz(3.002883) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8853569) q[1];
sx q[1];
rz(-2.0200122) q[1];
sx q[1];
rz(-3.0577809) q[1];
rz(-0.95529629) q[3];
sx q[3];
rz(-1.9855193) q[3];
sx q[3];
rz(2.6956468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.68219677) q[2];
sx q[2];
rz(-1.0707868) q[2];
sx q[2];
rz(-3.0219141) q[2];
rz(2.1648572) q[3];
sx q[3];
rz(-0.027438199) q[3];
sx q[3];
rz(-1.9500835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.479849) q[0];
sx q[0];
rz(-2.9763344) q[0];
sx q[0];
rz(1.4953493) q[0];
rz(-2.7961075) q[1];
sx q[1];
rz(-2.5472239) q[1];
sx q[1];
rz(-2.3631309) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0986106) q[0];
sx q[0];
rz(-1.6263282) q[0];
sx q[0];
rz(1.5776724) q[0];
x q[1];
rz(-2.3216281) q[2];
sx q[2];
rz(-1.1021309) q[2];
sx q[2];
rz(-2.6033273) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.567019) q[1];
sx q[1];
rz(-1.6751833) q[1];
sx q[1];
rz(2.0135572) q[1];
rz(-pi) q[2];
rz(2.4069294) q[3];
sx q[3];
rz(-2.0941995) q[3];
sx q[3];
rz(-1.911834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8567132) q[2];
sx q[2];
rz(-0.036980193) q[2];
sx q[2];
rz(-1.3899089) q[2];
rz(-0.022627929) q[3];
sx q[3];
rz(-0.32630625) q[3];
sx q[3];
rz(-2.7271395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49700272) q[0];
sx q[0];
rz(-3.1019326) q[0];
sx q[0];
rz(0.52763754) q[0];
rz(0.35609326) q[1];
sx q[1];
rz(-1.6368607) q[1];
sx q[1];
rz(-1.5632695) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1183076) q[0];
sx q[0];
rz(-1.0480799) q[0];
sx q[0];
rz(2.0705332) q[0];
rz(-pi) q[1];
rz(-3.061297) q[2];
sx q[2];
rz(-0.67932898) q[2];
sx q[2];
rz(1.7905362) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7303347) q[1];
sx q[1];
rz(-2.6739542) q[1];
sx q[1];
rz(1.1958666) q[1];
rz(-2.0464055) q[3];
sx q[3];
rz(-1.4970253) q[3];
sx q[3];
rz(-2.1493585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.57033527) q[2];
sx q[2];
rz(-0.0084849914) q[2];
sx q[2];
rz(-3.0961032) q[2];
rz(-0.82413834) q[3];
sx q[3];
rz(-0.51826158) q[3];
sx q[3];
rz(-0.99388188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6028041) q[0];
sx q[0];
rz(-2.0863918) q[0];
sx q[0];
rz(1.7428727) q[0];
rz(0.16821965) q[1];
sx q[1];
rz(-1.3361479) q[1];
sx q[1];
rz(2.9199563) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5835041) q[0];
sx q[0];
rz(-0.58204356) q[0];
sx q[0];
rz(2.01016) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8759769) q[2];
sx q[2];
rz(-0.22394315) q[2];
sx q[2];
rz(2.4502886) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0757743) q[1];
sx q[1];
rz(-1.5950959) q[1];
sx q[1];
rz(-1.1323117) q[1];
x q[2];
rz(-0.92042376) q[3];
sx q[3];
rz(-1.2022965) q[3];
sx q[3];
rz(-0.96741048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.98442709) q[2];
sx q[2];
rz(-0.027313622) q[2];
sx q[2];
rz(1.0013162) q[2];
rz(-2.2973513) q[3];
sx q[3];
rz(-0.068035754) q[3];
sx q[3];
rz(0.74725738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35185128) q[0];
sx q[0];
rz(-0.25775596) q[0];
sx q[0];
rz(-1.7716273) q[0];
rz(0.17579707) q[1];
sx q[1];
rz(-1.628123) q[1];
sx q[1];
rz(2.083185) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.269286) q[0];
sx q[0];
rz(-0.71505419) q[0];
sx q[0];
rz(2.0096409) q[0];
rz(-pi) q[1];
rz(-0.01137017) q[2];
sx q[2];
rz(-2.1595213) q[2];
sx q[2];
rz(-2.7259072) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5652983) q[1];
sx q[1];
rz(-0.75913069) q[1];
sx q[1];
rz(-2.7783093) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7806272) q[3];
sx q[3];
rz(-0.95300337) q[3];
sx q[3];
rz(1.3124026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.83219641) q[2];
sx q[2];
rz(-0.0038853566) q[2];
sx q[2];
rz(2.2957809) q[2];
rz(-2.5064365) q[3];
sx q[3];
rz(-0.35717765) q[3];
sx q[3];
rz(-0.14491189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7307067) q[0];
sx q[0];
rz(-2.9697953) q[0];
sx q[0];
rz(-0.26096499) q[0];
rz(-1.4313401) q[1];
sx q[1];
rz(-0.15161082) q[1];
sx q[1];
rz(1.8401015) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9985913) q[0];
sx q[0];
rz(-2.3218985) q[0];
sx q[0];
rz(-2.2065333) q[0];
rz(-pi) q[1];
rz(-0.13589047) q[2];
sx q[2];
rz(-0.11213485) q[2];
sx q[2];
rz(2.2035315) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3095418) q[1];
sx q[1];
rz(-2.0027035) q[1];
sx q[1];
rz(-1.5276272) q[1];
rz(-pi) q[2];
x q[2];
rz(1.306862) q[3];
sx q[3];
rz(-1.571705) q[3];
sx q[3];
rz(-3.0543882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6053551) q[2];
sx q[2];
rz(-1.321512) q[2];
sx q[2];
rz(-0.028667299) q[2];
rz(3.0981433) q[3];
sx q[3];
rz(-0.236792) q[3];
sx q[3];
rz(-1.4437599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0726149) q[0];
sx q[0];
rz(-2.7942939) q[0];
sx q[0];
rz(2.8436227) q[0];
rz(1.7295674) q[1];
sx q[1];
rz(-2.7822918) q[1];
sx q[1];
rz(-1.5493468) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9095347) q[0];
sx q[0];
rz(-1.8804872) q[0];
sx q[0];
rz(-2.880611) q[0];
x q[1];
rz(1.554073) q[2];
sx q[2];
rz(-1.5771167) q[2];
sx q[2];
rz(-2.9956452) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(5/(7*pi)) q[1];
sx q[1];
rz(-2.4691205) q[1];
sx q[1];
rz(0.48788957) q[1];
x q[2];
rz(1.8531591) q[3];
sx q[3];
rz(-1.3579277) q[3];
sx q[3];
rz(2.4119056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.90184244) q[2];
sx q[2];
rz(-0.032278927) q[2];
sx q[2];
rz(0.65995222) q[2];
rz(2.3615071) q[3];
sx q[3];
rz(-1.2939021) q[3];
sx q[3];
rz(-1.1680781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50769794) q[0];
sx q[0];
rz(-0.036490353) q[0];
sx q[0];
rz(2.3376035) q[0];
rz(-0.2313624) q[1];
sx q[1];
rz(-1.4120801) q[1];
sx q[1];
rz(-3.0661327) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20310878) q[0];
sx q[0];
rz(-1.9265947) q[0];
sx q[0];
rz(1.3514765) q[0];
rz(-1.6310591) q[2];
sx q[2];
rz(-0.00060877006) q[2];
sx q[2];
rz(-3.0920467) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4659404) q[1];
sx q[1];
rz(-0.91788252) q[1];
sx q[1];
rz(0.5799579) q[1];
x q[2];
rz(1.1832123) q[3];
sx q[3];
rz(-1.5389666) q[3];
sx q[3];
rz(-1.4393161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.017645322) q[2];
sx q[2];
rz(-2.6804774) q[2];
sx q[2];
rz(0.41575113) q[2];
rz(-1.4990643) q[3];
sx q[3];
rz(-3.1406904) q[3];
sx q[3];
rz(0.7846964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7163664) q[0];
sx q[0];
rz(-3.0854736) q[0];
sx q[0];
rz(2.8477493) q[0];
rz(1.6055239) q[1];
sx q[1];
rz(-0.87296456) q[1];
sx q[1];
rz(-1.6764838) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7820202) q[0];
sx q[0];
rz(-1.5806507) q[0];
sx q[0];
rz(0.1043679) q[0];
rz(-2.5757786) q[2];
sx q[2];
rz(-0.68556684) q[2];
sx q[2];
rz(-3.1367347) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5398354) q[1];
sx q[1];
rz(-0.64714995) q[1];
sx q[1];
rz(0.045250968) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1344101) q[3];
sx q[3];
rz(-1.5663538) q[3];
sx q[3];
rz(2.9053807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6426223) q[2];
sx q[2];
rz(-0.62752807) q[2];
sx q[2];
rz(1.9210531) q[2];
rz(2.7484861) q[3];
sx q[3];
rz(-3.1246287) q[3];
sx q[3];
rz(-0.19256798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4569693) q[0];
sx q[0];
rz(-1.4263117) q[0];
sx q[0];
rz(-0.34761467) q[0];
rz(-2.5034703) q[1];
sx q[1];
rz(-0.77824021) q[1];
sx q[1];
rz(2.4660769) q[1];
rz(-2.9478922) q[2];
sx q[2];
rz(-1.6483518) q[2];
sx q[2];
rz(-1.8087587) q[2];
rz(-1.5638696) q[3];
sx q[3];
rz(-1.6835811) q[3];
sx q[3];
rz(-1.2006105) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
