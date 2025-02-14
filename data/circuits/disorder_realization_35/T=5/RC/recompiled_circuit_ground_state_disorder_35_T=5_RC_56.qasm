OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72260296) q[0];
sx q[0];
rz(-2.498772) q[0];
sx q[0];
rz(8.2052054) q[0];
rz(3.976877) q[1];
sx q[1];
rz(4.4983954) q[1];
sx q[1];
rz(8.5228336) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25489461) q[0];
sx q[0];
rz(-0.6318081) q[0];
sx q[0];
rz(-2.8796701) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7576394) q[2];
sx q[2];
rz(-1.9451491) q[2];
sx q[2];
rz(0.49546212) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.343051) q[1];
sx q[1];
rz(-1.7447116) q[1];
sx q[1];
rz(1.2760722) q[1];
rz(0.91633032) q[3];
sx q[3];
rz(-0.95312762) q[3];
sx q[3];
rz(2.319782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1787662) q[2];
sx q[2];
rz(-2.1361394) q[2];
sx q[2];
rz(2.315305) q[2];
rz(1.4055584) q[3];
sx q[3];
rz(-0.30705753) q[3];
sx q[3];
rz(-2.3981986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58251441) q[0];
sx q[0];
rz(-2.7358416) q[0];
sx q[0];
rz(-1.7412809) q[0];
rz(1.1236313) q[1];
sx q[1];
rz(-0.39389899) q[1];
sx q[1];
rz(-1.0981015) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0535677) q[0];
sx q[0];
rz(-0.19101772) q[0];
sx q[0];
rz(0.6121401) q[0];
x q[1];
rz(2.3674521) q[2];
sx q[2];
rz(-0.4565276) q[2];
sx q[2];
rz(0.67195934) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.84867618) q[1];
sx q[1];
rz(-1.4956988) q[1];
sx q[1];
rz(0.1212705) q[1];
x q[2];
rz(1.9722749) q[3];
sx q[3];
rz(-0.40630925) q[3];
sx q[3];
rz(-0.10029785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0850247) q[2];
sx q[2];
rz(-2.8242064) q[2];
sx q[2];
rz(-1.8918096) q[2];
rz(1.362484) q[3];
sx q[3];
rz(-1.0008078) q[3];
sx q[3];
rz(1.6574297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5054841) q[0];
sx q[0];
rz(-2.4478069) q[0];
sx q[0];
rz(-0.27534494) q[0];
rz(-0.45922008) q[1];
sx q[1];
rz(-0.95454916) q[1];
sx q[1];
rz(3.088248) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5551852) q[0];
sx q[0];
rz(-1.4925028) q[0];
sx q[0];
rz(2.9281062) q[0];
x q[1];
rz(-3.1231327) q[2];
sx q[2];
rz(-2.0316796) q[2];
sx q[2];
rz(0.26319749) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0867696) q[1];
sx q[1];
rz(-2.5882698) q[1];
sx q[1];
rz(-0.077484681) q[1];
rz(-pi) q[2];
rz(0.20823222) q[3];
sx q[3];
rz(-1.8931812) q[3];
sx q[3];
rz(-0.65402189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3205545) q[2];
sx q[2];
rz(-0.36131636) q[2];
sx q[2];
rz(1.6374755) q[2];
rz(-1.5004246) q[3];
sx q[3];
rz(-2.0913561) q[3];
sx q[3];
rz(0.27568451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.7212873) q[0];
sx q[0];
rz(-0.5834226) q[0];
sx q[0];
rz(-0.48459184) q[0];
rz(1.2424319) q[1];
sx q[1];
rz(-0.74598765) q[1];
sx q[1];
rz(-0.34908435) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32767998) q[0];
sx q[0];
rz(-2.591344) q[0];
sx q[0];
rz(1.489086) q[0];
rz(0.075320638) q[2];
sx q[2];
rz(-1.2284797) q[2];
sx q[2];
rz(1.3665733) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.610747) q[1];
sx q[1];
rz(-1.7602923) q[1];
sx q[1];
rz(-2.9086472) q[1];
rz(-0.99557568) q[3];
sx q[3];
rz(-0.68505008) q[3];
sx q[3];
rz(-3.0127061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7829973) q[2];
sx q[2];
rz(-1.9390257) q[2];
sx q[2];
rz(-2.6915468) q[2];
rz(2.3027244) q[3];
sx q[3];
rz(-1.5476371) q[3];
sx q[3];
rz(2.94256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.49804509) q[0];
sx q[0];
rz(-1.4692551) q[0];
sx q[0];
rz(-3.041748) q[0];
rz(1.8210583) q[1];
sx q[1];
rz(-2.1456199) q[1];
sx q[1];
rz(2.5340396) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71220416) q[0];
sx q[0];
rz(-3.0297602) q[0];
sx q[0];
rz(-2.4080316) q[0];
x q[1];
rz(-0.25600146) q[2];
sx q[2];
rz(-2.5814272) q[2];
sx q[2];
rz(-1.3427116) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.95275003) q[1];
sx q[1];
rz(-0.61988587) q[1];
sx q[1];
rz(2.7623738) q[1];
rz(-1.2149548) q[3];
sx q[3];
rz(-2.8666593) q[3];
sx q[3];
rz(0.41809088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9472092) q[2];
sx q[2];
rz(-2.3781229) q[2];
sx q[2];
rz(-0.53059951) q[2];
rz(2.7001906) q[3];
sx q[3];
rz(-2.1680021) q[3];
sx q[3];
rz(-2.1921564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43668231) q[0];
sx q[0];
rz(-1.7969776) q[0];
sx q[0];
rz(2.5675024) q[0];
rz(-0.87999815) q[1];
sx q[1];
rz(-1.1209542) q[1];
sx q[1];
rz(-1.7990187) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3898252) q[0];
sx q[0];
rz(-2.6265567) q[0];
sx q[0];
rz(-1.9866605) q[0];
rz(-0.43812315) q[2];
sx q[2];
rz(-1.9480507) q[2];
sx q[2];
rz(-0.86710189) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.33730971) q[1];
sx q[1];
rz(-1.334467) q[1];
sx q[1];
rz(-1.1996891) q[1];
x q[2];
rz(1.430949) q[3];
sx q[3];
rz(-0.91250832) q[3];
sx q[3];
rz(1.1254532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4673956) q[2];
sx q[2];
rz(-0.28417045) q[2];
sx q[2];
rz(-2.6944842) q[2];
rz(2.1111264) q[3];
sx q[3];
rz(-1.6298529) q[3];
sx q[3];
rz(-0.38465056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7717188) q[0];
sx q[0];
rz(-1.8444703) q[0];
sx q[0];
rz(-0.41279992) q[0];
rz(1.075047) q[1];
sx q[1];
rz(-1.1378891) q[1];
sx q[1];
rz(0.80361754) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5301054) q[0];
sx q[0];
rz(-2.3854333) q[0];
sx q[0];
rz(2.3458781) q[0];
rz(2.2718524) q[2];
sx q[2];
rz(-0.77574965) q[2];
sx q[2];
rz(-2.9091331) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0688096) q[1];
sx q[1];
rz(-2.3549665) q[1];
sx q[1];
rz(1.5508316) q[1];
rz(-2.7429306) q[3];
sx q[3];
rz(-0.19426647) q[3];
sx q[3];
rz(0.27856058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2794118) q[2];
sx q[2];
rz(-1.2181686) q[2];
sx q[2];
rz(1.5009521) q[2];
rz(2.5189279) q[3];
sx q[3];
rz(-0.77087918) q[3];
sx q[3];
rz(-1.6525432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.415446) q[0];
sx q[0];
rz(-1.6212689) q[0];
sx q[0];
rz(-0.55794445) q[0];
rz(-0.56124148) q[1];
sx q[1];
rz(-1.7702425) q[1];
sx q[1];
rz(1.7254613) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19085267) q[0];
sx q[0];
rz(-1.2878364) q[0];
sx q[0];
rz(-2.0377732) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2165856) q[2];
sx q[2];
rz(-2.6753798) q[2];
sx q[2];
rz(0.031161873) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.76745236) q[1];
sx q[1];
rz(-1.472964) q[1];
sx q[1];
rz(1.2409452) q[1];
rz(1.775557) q[3];
sx q[3];
rz(-0.8258709) q[3];
sx q[3];
rz(1.5457758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.8884362) q[2];
sx q[2];
rz(-2.3641868) q[2];
sx q[2];
rz(-1.0003132) q[2];
rz(0.12217626) q[3];
sx q[3];
rz(-1.5001985) q[3];
sx q[3];
rz(-0.15672556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4633453) q[0];
sx q[0];
rz(-1.1868917) q[0];
sx q[0];
rz(-0.004322411) q[0];
rz(-0.16383544) q[1];
sx q[1];
rz(-2.7163353) q[1];
sx q[1];
rz(1.3660627) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3177494) q[0];
sx q[0];
rz(-1.5429521) q[0];
sx q[0];
rz(2.2180024) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4539656) q[2];
sx q[2];
rz(-2.673871) q[2];
sx q[2];
rz(-2.2224838) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8722265) q[1];
sx q[1];
rz(-2.8878643) q[1];
sx q[1];
rz(2.9969508) q[1];
x q[2];
rz(-2.4716464) q[3];
sx q[3];
rz(-2.1996227) q[3];
sx q[3];
rz(2.3184702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1395448) q[2];
sx q[2];
rz(-2.1740972) q[2];
sx q[2];
rz(0.61919332) q[2];
rz(-2.5891417) q[3];
sx q[3];
rz(-1.3055472) q[3];
sx q[3];
rz(2.9964871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.3510975) q[0];
sx q[0];
rz(-2.9572697) q[0];
sx q[0];
rz(-0.47875324) q[0];
rz(1.152773) q[1];
sx q[1];
rz(-2.8515127) q[1];
sx q[1];
rz(-2.1943888) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21175948) q[0];
sx q[0];
rz(-2.9538028) q[0];
sx q[0];
rz(-1.7540356) q[0];
rz(-pi) q[1];
rz(-0.98299111) q[2];
sx q[2];
rz(-1.3653367) q[2];
sx q[2];
rz(0.062177156) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.80120197) q[1];
sx q[1];
rz(-1.2653192) q[1];
sx q[1];
rz(-0.82251541) q[1];
rz(-1.8197038) q[3];
sx q[3];
rz(-1.9134023) q[3];
sx q[3];
rz(3.1150888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7572215) q[2];
sx q[2];
rz(-0.20558509) q[2];
sx q[2];
rz(0.26665404) q[2];
rz(-2.0742778) q[3];
sx q[3];
rz(-1.9284733) q[3];
sx q[3];
rz(-0.96323577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(3.0530479) q[0];
sx q[0];
rz(-1.5374669) q[0];
sx q[0];
rz(1.5446825) q[0];
rz(-2.0883941) q[1];
sx q[1];
rz(-1.9269301) q[1];
sx q[1];
rz(-2.3763837) q[1];
rz(1.5736035) q[2];
sx q[2];
rz(-1.8127828) q[2];
sx q[2];
rz(0.91290963) q[2];
rz(2.8070634) q[3];
sx q[3];
rz(-0.60548895) q[3];
sx q[3];
rz(2.0169712) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
