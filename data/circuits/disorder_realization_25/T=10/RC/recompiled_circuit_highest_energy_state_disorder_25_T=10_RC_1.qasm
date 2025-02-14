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
rz(0.99474466) q[0];
sx q[0];
rz(-2.5543307) q[0];
sx q[0];
rz(0.62556148) q[0];
rz(-2.5211531) q[1];
sx q[1];
rz(-0.99045366) q[1];
sx q[1];
rz(2.3576696) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1358937) q[0];
sx q[0];
rz(-1.5672475) q[0];
sx q[0];
rz(-0.43412848) q[0];
rz(-1.6541846) q[2];
sx q[2];
rz(-1.4321799) q[2];
sx q[2];
rz(-2.1237789) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.82134889) q[1];
sx q[1];
rz(-2.26198) q[1];
sx q[1];
rz(-0.76637474) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0524527) q[3];
sx q[3];
rz(-2.3618792) q[3];
sx q[3];
rz(2.4453795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1938532) q[2];
sx q[2];
rz(-2.096602) q[2];
sx q[2];
rz(-0.6673153) q[2];
rz(-2.8269178) q[3];
sx q[3];
rz(-0.59942013) q[3];
sx q[3];
rz(-2.2144894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2941403) q[0];
sx q[0];
rz(-2.2951234) q[0];
sx q[0];
rz(0.2997998) q[0];
rz(-1.0046129) q[1];
sx q[1];
rz(-0.71669465) q[1];
sx q[1];
rz(0.85003781) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3687821) q[0];
sx q[0];
rz(-1.6666935) q[0];
sx q[0];
rz(0.96536388) q[0];
rz(-pi) q[1];
rz(-2.9710567) q[2];
sx q[2];
rz(-1.4922304) q[2];
sx q[2];
rz(-0.54135347) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5932754) q[1];
sx q[1];
rz(-1.9276345) q[1];
sx q[1];
rz(0.22290454) q[1];
x q[2];
rz(1.5845234) q[3];
sx q[3];
rz(-1.9753595) q[3];
sx q[3];
rz(1.3409529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5726418) q[2];
sx q[2];
rz(-0.53223842) q[2];
sx q[2];
rz(2.4786095) q[2];
rz(2.3891383) q[3];
sx q[3];
rz(-1.821725) q[3];
sx q[3];
rz(-0.19645709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3982584) q[0];
sx q[0];
rz(-2.7404116) q[0];
sx q[0];
rz(-2.4617526) q[0];
rz(-0.72194779) q[1];
sx q[1];
rz(-2.5847021) q[1];
sx q[1];
rz(-2.360875) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74405438) q[0];
sx q[0];
rz(-2.7408532) q[0];
sx q[0];
rz(1.5903872) q[0];
rz(1.5268097) q[2];
sx q[2];
rz(-3.0553748) q[2];
sx q[2];
rz(1.1607064) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.91619766) q[1];
sx q[1];
rz(-2.19529) q[1];
sx q[1];
rz(1.1033113) q[1];
x q[2];
rz(-1.8119007) q[3];
sx q[3];
rz(-0.74987292) q[3];
sx q[3];
rz(-2.1151224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9146933) q[2];
sx q[2];
rz(-1.2279899) q[2];
sx q[2];
rz(2.4719888) q[2];
rz(-2.2412444) q[3];
sx q[3];
rz(-1.0122296) q[3];
sx q[3];
rz(0.77696925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25903073) q[0];
sx q[0];
rz(-2.7029523) q[0];
sx q[0];
rz(-0.90748179) q[0];
rz(0.59440815) q[1];
sx q[1];
rz(-2.2120357) q[1];
sx q[1];
rz(-2.0118735) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66624712) q[0];
sx q[0];
rz(-0.94303688) q[0];
sx q[0];
rz(-0.23601377) q[0];
x q[1];
rz(-1.8215034) q[2];
sx q[2];
rz(-2.5959229) q[2];
sx q[2];
rz(0.16599338) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.65986827) q[1];
sx q[1];
rz(-2.4818083) q[1];
sx q[1];
rz(-2.6901812) q[1];
rz(2.8606671) q[3];
sx q[3];
rz(-1.8531898) q[3];
sx q[3];
rz(2.1641233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0323459) q[2];
sx q[2];
rz(-1.3999908) q[2];
sx q[2];
rz(-2.4692811) q[2];
rz(0.0191056) q[3];
sx q[3];
rz(-0.34135434) q[3];
sx q[3];
rz(-3.0066971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8371209) q[0];
sx q[0];
rz(-0.79455513) q[0];
sx q[0];
rz(2.9735612) q[0];
rz(-1.2270323) q[1];
sx q[1];
rz(-1.9214168) q[1];
sx q[1];
rz(-0.29774818) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3963602) q[0];
sx q[0];
rz(-2.0848635) q[0];
sx q[0];
rz(-0.6020418) q[0];
x q[1];
rz(-2.3565294) q[2];
sx q[2];
rz(-0.74399555) q[2];
sx q[2];
rz(-2.1659578) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1923026) q[1];
sx q[1];
rz(-1.5459035) q[1];
sx q[1];
rz(-0.20583238) q[1];
rz(-pi) q[2];
rz(-1.566542) q[3];
sx q[3];
rz(-2.9828072) q[3];
sx q[3];
rz(-2.6095853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9510522) q[2];
sx q[2];
rz(-1.0995883) q[2];
sx q[2];
rz(2.5229048) q[2];
rz(-2.3330073) q[3];
sx q[3];
rz(-2.6638668) q[3];
sx q[3];
rz(-1.3396858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65281463) q[0];
sx q[0];
rz(-1.6365016) q[0];
sx q[0];
rz(2.6763647) q[0];
rz(-2.6564927) q[1];
sx q[1];
rz(-2.7310889) q[1];
sx q[1];
rz(0.088265158) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8387091) q[0];
sx q[0];
rz(-1.1195991) q[0];
sx q[0];
rz(2.124435) q[0];
rz(1.5563929) q[2];
sx q[2];
rz(-0.72897899) q[2];
sx q[2];
rz(-0.62504238) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8264721) q[1];
sx q[1];
rz(-1.406028) q[1];
sx q[1];
rz(-0.75854782) q[1];
rz(2.8507455) q[3];
sx q[3];
rz(-2.8418782) q[3];
sx q[3];
rz(-0.3464454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7469067) q[2];
sx q[2];
rz(-0.8553108) q[2];
sx q[2];
rz(-0.49760231) q[2];
rz(0.48007128) q[3];
sx q[3];
rz(-2.2205455) q[3];
sx q[3];
rz(-0.068643071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8467872) q[0];
sx q[0];
rz(-0.46600309) q[0];
sx q[0];
rz(-1.4303327) q[0];
rz(-1.3149186) q[1];
sx q[1];
rz(-2.7480835) q[1];
sx q[1];
rz(1.5865145) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85077943) q[0];
sx q[0];
rz(-1.4516136) q[0];
sx q[0];
rz(-1.9656154) q[0];
x q[1];
rz(2.2261593) q[2];
sx q[2];
rz(-0.829773) q[2];
sx q[2];
rz(2.1754153) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3027398) q[1];
sx q[1];
rz(-1.4697971) q[1];
sx q[1];
rz(-0.98827004) q[1];
x q[2];
rz(2.2022543) q[3];
sx q[3];
rz(-1.5087481) q[3];
sx q[3];
rz(0.59988672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2985208) q[2];
sx q[2];
rz(-0.71121794) q[2];
sx q[2];
rz(-2.2141875) q[2];
rz(1.1585506) q[3];
sx q[3];
rz(-2.4535593) q[3];
sx q[3];
rz(0.096175171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7189099) q[0];
sx q[0];
rz(-0.89253187) q[0];
sx q[0];
rz(-2.6655777) q[0];
rz(-2.1503275) q[1];
sx q[1];
rz(-0.74949336) q[1];
sx q[1];
rz(-3.057726) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1779685) q[0];
sx q[0];
rz(-1.9289079) q[0];
sx q[0];
rz(1.0492508) q[0];
rz(0.65292439) q[2];
sx q[2];
rz(-2.5819375) q[2];
sx q[2];
rz(1.4687572) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2460275) q[1];
sx q[1];
rz(-1.152134) q[1];
sx q[1];
rz(0.093105969) q[1];
rz(1.7937035) q[3];
sx q[3];
rz(-2.768424) q[3];
sx q[3];
rz(2.580433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.91741651) q[2];
sx q[2];
rz(-0.26599628) q[2];
sx q[2];
rz(-0.64845294) q[2];
rz(2.2367541) q[3];
sx q[3];
rz(-1.4584533) q[3];
sx q[3];
rz(-0.31074935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34173486) q[0];
sx q[0];
rz(-0.71001995) q[0];
sx q[0];
rz(0.56890666) q[0];
rz(-2.3078602) q[1];
sx q[1];
rz(-1.2833475) q[1];
sx q[1];
rz(2.1368829) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92623989) q[0];
sx q[0];
rz(-0.66074255) q[0];
sx q[0];
rz(2.8244429) q[0];
rz(-pi) q[1];
rz(-1.1962074) q[2];
sx q[2];
rz(-2.5718712) q[2];
sx q[2];
rz(1.8853055) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1596327) q[1];
sx q[1];
rz(-1.4985604) q[1];
sx q[1];
rz(2.2700929) q[1];
x q[2];
rz(0.58752693) q[3];
sx q[3];
rz(-1.6391132) q[3];
sx q[3];
rz(0.63623601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4363165) q[2];
sx q[2];
rz(-0.87203163) q[2];
sx q[2];
rz(-0.3479859) q[2];
rz(-2.3157388) q[3];
sx q[3];
rz(-0.21819849) q[3];
sx q[3];
rz(-0.60514778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0843049) q[0];
sx q[0];
rz(-2.075752) q[0];
sx q[0];
rz(2.9344946) q[0];
rz(-0.53889489) q[1];
sx q[1];
rz(-2.5205044) q[1];
sx q[1];
rz(-0.60212392) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3823366) q[0];
sx q[0];
rz(-0.70092541) q[0];
sx q[0];
rz(-0.78030326) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5536898) q[2];
sx q[2];
rz(-1.9155115) q[2];
sx q[2];
rz(1.6893139) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94127266) q[1];
sx q[1];
rz(-1.6459904) q[1];
sx q[1];
rz(-1.1635416) q[1];
rz(-pi) q[2];
rz(-1.243351) q[3];
sx q[3];
rz(-2.2909938) q[3];
sx q[3];
rz(2.7011288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.78486097) q[2];
sx q[2];
rz(-0.95180231) q[2];
sx q[2];
rz(1.7968563) q[2];
rz(-0.52062672) q[3];
sx q[3];
rz(-2.8713363) q[3];
sx q[3];
rz(0.80210137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6945334) q[0];
sx q[0];
rz(-2.4166528) q[0];
sx q[0];
rz(1.7550533) q[0];
rz(-2.3583892) q[1];
sx q[1];
rz(-1.422311) q[1];
sx q[1];
rz(1.5625988) q[1];
rz(2.798525) q[2];
sx q[2];
rz(-1.82556) q[2];
sx q[2];
rz(-1.0221046) q[2];
rz(0.24012031) q[3];
sx q[3];
rz(-1.5821531) q[3];
sx q[3];
rz(1.5021363) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
