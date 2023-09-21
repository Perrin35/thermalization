OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4560661) q[0];
sx q[0];
rz(-0.38903061) q[0];
sx q[0];
rz(2.2580137) q[0];
rz(3.1318624) q[1];
sx q[1];
rz(-1.6844123) q[1];
sx q[1];
rz(-1.943346) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0966914) q[0];
sx q[0];
rz(-2.9400819) q[0];
sx q[0];
rz(2.6411396) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2486357) q[2];
sx q[2];
rz(-0.82167168) q[2];
sx q[2];
rz(1.2920213) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7361684) q[1];
sx q[1];
rz(-0.6286469) q[1];
sx q[1];
rz(-2.5581762) q[1];
x q[2];
rz(1.4000113) q[3];
sx q[3];
rz(-1.5251953) q[3];
sx q[3];
rz(0.036269773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1951695) q[2];
sx q[2];
rz(-2.158458) q[2];
sx q[2];
rz(-2.9602489) q[2];
rz(-0.26120734) q[3];
sx q[3];
rz(-1.2657335) q[3];
sx q[3];
rz(-0.75631022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8392035) q[0];
sx q[0];
rz(-0.2897245) q[0];
sx q[0];
rz(2.7547577) q[0];
rz(2.6392) q[1];
sx q[1];
rz(-0.97351176) q[1];
sx q[1];
rz(-1.5997255) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4138448) q[0];
sx q[0];
rz(-1.9808597) q[0];
sx q[0];
rz(-2.4983665) q[0];
rz(-pi) q[1];
rz(-0.16427152) q[2];
sx q[2];
rz(-2.0220244) q[2];
sx q[2];
rz(2.3300366) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4914354) q[1];
sx q[1];
rz(-0.037089247) q[1];
sx q[1];
rz(-0.27122072) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3719671) q[3];
sx q[3];
rz(-0.63109055) q[3];
sx q[3];
rz(-1.0069932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5388422) q[2];
sx q[2];
rz(-2.2637612) q[2];
sx q[2];
rz(1.7269469) q[2];
rz(0.85033068) q[3];
sx q[3];
rz(-2.7089705) q[3];
sx q[3];
rz(-1.6833646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26281115) q[0];
sx q[0];
rz(-1.6101863) q[0];
sx q[0];
rz(-0.61022726) q[0];
rz(-1.8521076) q[1];
sx q[1];
rz(-0.97924966) q[1];
sx q[1];
rz(2.1496225) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45527601) q[0];
sx q[0];
rz(-2.1222097) q[0];
sx q[0];
rz(-2.9247012) q[0];
rz(-pi) q[1];
rz(1.7730373) q[2];
sx q[2];
rz(-2.2291406) q[2];
sx q[2];
rz(-0.973268) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.62832824) q[1];
sx q[1];
rz(-0.4497512) q[1];
sx q[1];
rz(0.20723923) q[1];
rz(-1.8157585) q[3];
sx q[3];
rz(-1.1253329) q[3];
sx q[3];
rz(0.20416343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.38301864) q[2];
sx q[2];
rz(-0.16123161) q[2];
sx q[2];
rz(0.26829159) q[2];
rz(-0.39408436) q[3];
sx q[3];
rz(-1.2309309) q[3];
sx q[3];
rz(-2.9530853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85369337) q[0];
sx q[0];
rz(-2.6278966) q[0];
sx q[0];
rz(0.40507856) q[0];
rz(2.4515117) q[1];
sx q[1];
rz(-1.9837374) q[1];
sx q[1];
rz(2.4437723) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21701248) q[0];
sx q[0];
rz(-3.0487061) q[0];
sx q[0];
rz(1.3148944) q[0];
x q[1];
rz(1.1890829) q[2];
sx q[2];
rz(-2.082798) q[2];
sx q[2];
rz(-0.94101671) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7792369) q[1];
sx q[1];
rz(-0.52473611) q[1];
sx q[1];
rz(-1.1974105) q[1];
x q[2];
rz(0.29699765) q[3];
sx q[3];
rz(-1.1668491) q[3];
sx q[3];
rz(-1.3511853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4219389) q[2];
sx q[2];
rz(-1.3976588) q[2];
sx q[2];
rz(1.5092124) q[2];
rz(0.40431067) q[3];
sx q[3];
rz(-2.4590838) q[3];
sx q[3];
rz(1.6633165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25800911) q[0];
sx q[0];
rz(-1.785935) q[0];
sx q[0];
rz(-2.1160545) q[0];
rz(-0.57199663) q[1];
sx q[1];
rz(-2.0472186) q[1];
sx q[1];
rz(0.62932032) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48110163) q[0];
sx q[0];
rz(-1.8697303) q[0];
sx q[0];
rz(0.42994182) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2467381) q[2];
sx q[2];
rz(-1.7712777) q[2];
sx q[2];
rz(-1.5948053) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1197966) q[1];
sx q[1];
rz(-2.6152059) q[1];
sx q[1];
rz(0.22967931) q[1];
x q[2];
rz(2.1719645) q[3];
sx q[3];
rz(-1.3519577) q[3];
sx q[3];
rz(0.18274433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5489674) q[2];
sx q[2];
rz(-2.7719345) q[2];
sx q[2];
rz(-0.34234753) q[2];
rz(-1.7957548) q[3];
sx q[3];
rz(-1.6941518) q[3];
sx q[3];
rz(-0.19601823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.489007) q[0];
sx q[0];
rz(-1.9056029) q[0];
sx q[0];
rz(0.18519369) q[0];
rz(1.406503) q[1];
sx q[1];
rz(-1.0909189) q[1];
sx q[1];
rz(1.7746183) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0104116) q[0];
sx q[0];
rz(-2.512012) q[0];
sx q[0];
rz(0.38139831) q[0];
rz(-pi) q[1];
rz(2.30079) q[2];
sx q[2];
rz(-1.3422988) q[2];
sx q[2];
rz(2.218354) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0095014) q[1];
sx q[1];
rz(-0.17772929) q[1];
sx q[1];
rz(-1.0264261) q[1];
x q[2];
rz(-1.9983415) q[3];
sx q[3];
rz(-1.0435836) q[3];
sx q[3];
rz(-0.15641071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.24484816) q[2];
sx q[2];
rz(-0.5534133) q[2];
sx q[2];
rz(0.883376) q[2];
rz(1.4128489) q[3];
sx q[3];
rz(-0.69245517) q[3];
sx q[3];
rz(-2.3278055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.8900523) q[0];
sx q[0];
rz(-0.10245704) q[0];
sx q[0];
rz(1.2782156) q[0];
rz(-3.1037519) q[1];
sx q[1];
rz(-2.3262639) q[1];
sx q[1];
rz(1.7657123) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85526953) q[0];
sx q[0];
rz(-0.77676847) q[0];
sx q[0];
rz(1.1458678) q[0];
rz(-pi) q[1];
rz(-1.0972294) q[2];
sx q[2];
rz(-1.5175022) q[2];
sx q[2];
rz(-3.1380944) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.618305) q[1];
sx q[1];
rz(-1.6860645) q[1];
sx q[1];
rz(1.6769888) q[1];
rz(1.2405346) q[3];
sx q[3];
rz(-1.5152647) q[3];
sx q[3];
rz(-2.8020669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6841131) q[2];
sx q[2];
rz(-2.423968) q[2];
sx q[2];
rz(2.4105371) q[2];
rz(-3.030792) q[3];
sx q[3];
rz(-1.5564857) q[3];
sx q[3];
rz(2.4462162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5450127) q[0];
sx q[0];
rz(-0.87485635) q[0];
sx q[0];
rz(-2.4080283) q[0];
rz(-0.60797524) q[1];
sx q[1];
rz(-1.1939476) q[1];
sx q[1];
rz(0.2342934) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60093555) q[0];
sx q[0];
rz(-1.0957076) q[0];
sx q[0];
rz(0.29151543) q[0];
x q[1];
rz(0.73306002) q[2];
sx q[2];
rz(-0.86793938) q[2];
sx q[2];
rz(1.5066063) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.56470796) q[1];
sx q[1];
rz(-1.7600093) q[1];
sx q[1];
rz(-2.3343759) q[1];
x q[2];
rz(-0.33741823) q[3];
sx q[3];
rz(-2.1302345) q[3];
sx q[3];
rz(-2.9448201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0155448) q[2];
sx q[2];
rz(-1.8323106) q[2];
sx q[2];
rz(1.7162494) q[2];
rz(1.4632633) q[3];
sx q[3];
rz(-2.3571456) q[3];
sx q[3];
rz(-1.4956168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54206806) q[0];
sx q[0];
rz(-2.8064089) q[0];
sx q[0];
rz(1.9482127) q[0];
rz(-1.2606196) q[1];
sx q[1];
rz(-1.7766989) q[1];
sx q[1];
rz(2.1967922) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017100447) q[0];
sx q[0];
rz(-1.955535) q[0];
sx q[0];
rz(2.4186633) q[0];
rz(1.9656885) q[2];
sx q[2];
rz(-1.2972304) q[2];
sx q[2];
rz(0.045189518) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.33465696) q[1];
sx q[1];
rz(-1.2922704) q[1];
sx q[1];
rz(1.7486228) q[1];
rz(-1.9267843) q[3];
sx q[3];
rz(-1.755852) q[3];
sx q[3];
rz(1.9389648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9251359) q[2];
sx q[2];
rz(-2.0704806) q[2];
sx q[2];
rz(-2.8533868) q[2];
rz(2.6618585) q[3];
sx q[3];
rz(-2.0917442) q[3];
sx q[3];
rz(-1.6335999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0297246) q[0];
sx q[0];
rz(-2.8653963) q[0];
sx q[0];
rz(0.9129886) q[0];
rz(-2.7669725) q[1];
sx q[1];
rz(-1.4034142) q[1];
sx q[1];
rz(0.8909117) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3324146) q[0];
sx q[0];
rz(-0.42254585) q[0];
sx q[0];
rz(2.4992864) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28384039) q[2];
sx q[2];
rz(-2.1643057) q[2];
sx q[2];
rz(-0.80988353) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3799694) q[1];
sx q[1];
rz(-2.7459811) q[1];
sx q[1];
rz(-2.604894) q[1];
rz(2.9534146) q[3];
sx q[3];
rz(-1.4965701) q[3];
sx q[3];
rz(2.1256413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9051819) q[2];
sx q[2];
rz(-0.85835251) q[2];
sx q[2];
rz(-2.6399844) q[2];
rz(1.2891399) q[3];
sx q[3];
rz(-1.4533549) q[3];
sx q[3];
rz(-2.6954209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5065153) q[0];
sx q[0];
rz(-1.4415393) q[0];
sx q[0];
rz(-2.517979) q[0];
rz(2.0093909) q[1];
sx q[1];
rz(-2.3846346) q[1];
sx q[1];
rz(0.089288575) q[1];
rz(-1.7239465) q[2];
sx q[2];
rz(-1.8344804) q[2];
sx q[2];
rz(-1.9305965) q[2];
rz(1.489747) q[3];
sx q[3];
rz(-0.54809082) q[3];
sx q[3];
rz(-2.5561668) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
