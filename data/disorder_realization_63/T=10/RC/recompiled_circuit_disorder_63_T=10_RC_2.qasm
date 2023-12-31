OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1749984) q[0];
sx q[0];
rz(-0.35342616) q[0];
sx q[0];
rz(1.0647635) q[0];
rz(3.937768) q[1];
sx q[1];
rz(1.9328971) q[1];
sx q[1];
rz(9.9608496) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6807841) q[0];
sx q[0];
rz(-0.55235282) q[0];
sx q[0];
rz(-2.0244563) q[0];
rz(-3.0739003) q[2];
sx q[2];
rz(-2.2798685) q[2];
sx q[2];
rz(0.46082218) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2823892) q[1];
sx q[1];
rz(-2.5232362) q[1];
sx q[1];
rz(-2.8137141) q[1];
rz(-pi) q[2];
rz(1.0814813) q[3];
sx q[3];
rz(-0.89852528) q[3];
sx q[3];
rz(1.0244474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7321695) q[2];
sx q[2];
rz(-2.7004341) q[2];
sx q[2];
rz(-2.1146963) q[2];
rz(-0.25201592) q[3];
sx q[3];
rz(-1.1427294) q[3];
sx q[3];
rz(2.2759329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082829647) q[0];
sx q[0];
rz(-2.7936462) q[0];
sx q[0];
rz(1.7513562) q[0];
rz(0.83579666) q[1];
sx q[1];
rz(-2.4048769) q[1];
sx q[1];
rz(0.70835152) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0396597) q[0];
sx q[0];
rz(-1.205737) q[0];
sx q[0];
rz(1.2714766) q[0];
rz(-pi) q[1];
rz(1.7827665) q[2];
sx q[2];
rz(-2.5296629) q[2];
sx q[2];
rz(1.3509392) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0185512) q[1];
sx q[1];
rz(-1.0804783) q[1];
sx q[1];
rz(-0.80177387) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3362818) q[3];
sx q[3];
rz(-0.66262965) q[3];
sx q[3];
rz(0.90028541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.11671242) q[2];
sx q[2];
rz(-0.79170266) q[2];
sx q[2];
rz(-2.8498245) q[2];
rz(0.10270384) q[3];
sx q[3];
rz(-1.4029968) q[3];
sx q[3];
rz(-1.6170988) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8125238) q[0];
sx q[0];
rz(-3.0561495) q[0];
sx q[0];
rz(-2.8318751) q[0];
rz(-1.4801056) q[1];
sx q[1];
rz(-1.3269576) q[1];
sx q[1];
rz(0.57166878) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088061995) q[0];
sx q[0];
rz(-2.5998839) q[0];
sx q[0];
rz(0.17266973) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0057726) q[2];
sx q[2];
rz(-1.2425213) q[2];
sx q[2];
rz(-0.091718397) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7844312) q[1];
sx q[1];
rz(-2.2964587) q[1];
sx q[1];
rz(1.7482589) q[1];
rz(-pi) q[2];
rz(1.709135) q[3];
sx q[3];
rz(-1.1225015) q[3];
sx q[3];
rz(-1.5642779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.16033515) q[2];
sx q[2];
rz(-2.2312639) q[2];
sx q[2];
rz(2.4327915) q[2];
rz(-2.839084) q[3];
sx q[3];
rz(-1.442028) q[3];
sx q[3];
rz(-1.9992874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(-0.70808327) q[0];
sx q[0];
rz(-1.1129365) q[0];
sx q[0];
rz(-2.3420912) q[0];
rz(-3.0917621) q[1];
sx q[1];
rz(-0.89490503) q[1];
sx q[1];
rz(0.18049151) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7364396) q[0];
sx q[0];
rz(-0.3973876) q[0];
sx q[0];
rz(2.9050164) q[0];
rz(-1.7454342) q[2];
sx q[2];
rz(-1.9348113) q[2];
sx q[2];
rz(0.82496914) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6199477) q[1];
sx q[1];
rz(-1.4536899) q[1];
sx q[1];
rz(-2.820693) q[1];
rz(-1.0749712) q[3];
sx q[3];
rz(-1.5815445) q[3];
sx q[3];
rz(1.4533952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.089347) q[2];
sx q[2];
rz(-1.1648488) q[2];
sx q[2];
rz(-1.5578516) q[2];
rz(1.3752939) q[3];
sx q[3];
rz(-1.3170653) q[3];
sx q[3];
rz(-2.2089675) q[3];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3809526) q[0];
sx q[0];
rz(-2.8044658) q[0];
sx q[0];
rz(-2.1550762) q[0];
rz(1.9793234) q[1];
sx q[1];
rz(-1.9245851) q[1];
sx q[1];
rz(-2.8900237) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4261599) q[0];
sx q[0];
rz(-0.71007219) q[0];
sx q[0];
rz(-1.8415113) q[0];
rz(-2.9676874) q[2];
sx q[2];
rz(-2.675563) q[2];
sx q[2];
rz(-2.5582061) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0136452) q[1];
sx q[1];
rz(-0.19731678) q[1];
sx q[1];
rz(-0.13983388) q[1];
rz(-pi) q[2];
rz(2.1523347) q[3];
sx q[3];
rz(-0.79205081) q[3];
sx q[3];
rz(1.1359515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.183737) q[2];
sx q[2];
rz(-0.541406) q[2];
sx q[2];
rz(-1.0305369) q[2];
rz(1.7306227) q[3];
sx q[3];
rz(-2.1822699) q[3];
sx q[3];
rz(-2.5103536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75491607) q[0];
sx q[0];
rz(-2.6658391) q[0];
sx q[0];
rz(-2.2914698) q[0];
rz(1.1823581) q[1];
sx q[1];
rz(-1.9071002) q[1];
sx q[1];
rz(0.39302557) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61009208) q[0];
sx q[0];
rz(-1.2628265) q[0];
sx q[0];
rz(0.17225762) q[0];
rz(-pi) q[1];
rz(0.23586258) q[2];
sx q[2];
rz(-2.3013407) q[2];
sx q[2];
rz(-2.0116531) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.45822083) q[1];
sx q[1];
rz(-2.4903989) q[1];
sx q[1];
rz(1.5797735) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7418467) q[3];
sx q[3];
rz(-2.7478455) q[3];
sx q[3];
rz(1.1175931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7810016) q[2];
sx q[2];
rz(-2.0157308) q[2];
sx q[2];
rz(0.97314107) q[2];
rz(-0.94758236) q[3];
sx q[3];
rz(-2.0411453) q[3];
sx q[3];
rz(1.9472286) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4457552) q[0];
sx q[0];
rz(-2.7212302) q[0];
sx q[0];
rz(-1.6181035) q[0];
rz(0.37480005) q[1];
sx q[1];
rz(-1.260489) q[1];
sx q[1];
rz(-3.135625) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8864266) q[0];
sx q[0];
rz(-0.92029858) q[0];
sx q[0];
rz(-2.7889473) q[0];
rz(0.75255021) q[2];
sx q[2];
rz(-1.9117022) q[2];
sx q[2];
rz(2.7698851) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5827427) q[1];
sx q[1];
rz(-2.5834281) q[1];
sx q[1];
rz(1.4504257) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0217516) q[3];
sx q[3];
rz(-0.37444515) q[3];
sx q[3];
rz(-3.0344809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.57198793) q[2];
sx q[2];
rz(-0.54509744) q[2];
sx q[2];
rz(-2.9515284) q[2];
rz(-2.7251785) q[3];
sx q[3];
rz(-2.1834686) q[3];
sx q[3];
rz(2.4068508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1104601) q[0];
sx q[0];
rz(-1.0460331) q[0];
sx q[0];
rz(0.53623143) q[0];
rz(0.93332943) q[1];
sx q[1];
rz(-1.5678762) q[1];
sx q[1];
rz(2.1933864) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3535239) q[0];
sx q[0];
rz(-0.81809645) q[0];
sx q[0];
rz(1.3226932) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6981632) q[2];
sx q[2];
rz(-2.0460874) q[2];
sx q[2];
rz(-2.7623917) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.58029786) q[1];
sx q[1];
rz(-1.5082238) q[1];
sx q[1];
rz(2.9194174) q[1];
rz(-pi) q[2];
rz(-2.7847399) q[3];
sx q[3];
rz(-1.4992504) q[3];
sx q[3];
rz(2.1944291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5006717) q[2];
sx q[2];
rz(-1.6225953) q[2];
sx q[2];
rz(-2.9013157) q[2];
rz(0.37929532) q[3];
sx q[3];
rz(-2.1160188) q[3];
sx q[3];
rz(-1.3482288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3786479) q[0];
sx q[0];
rz(-2.8084016) q[0];
sx q[0];
rz(2.8906524) q[0];
rz(2.8885686) q[1];
sx q[1];
rz(-1.7605942) q[1];
sx q[1];
rz(-2.7889263) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1044554) q[0];
sx q[0];
rz(-1.4482575) q[0];
sx q[0];
rz(3.0260968) q[0];
x q[1];
rz(-2.7025931) q[2];
sx q[2];
rz(-2.8121431) q[2];
sx q[2];
rz(-1.1169369) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.36665146) q[1];
sx q[1];
rz(-1.0171486) q[1];
sx q[1];
rz(1.8371546) q[1];
rz(2.6412233) q[3];
sx q[3];
rz(-0.89958588) q[3];
sx q[3];
rz(2.3994115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3744366) q[2];
sx q[2];
rz(-1.1194976) q[2];
sx q[2];
rz(-0.68332589) q[2];
rz(-1.654489) q[3];
sx q[3];
rz(-2.7023102) q[3];
sx q[3];
rz(1.9911511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3764573) q[0];
sx q[0];
rz(-0.52381223) q[0];
sx q[0];
rz(-1.3002243) q[0];
rz(-0.70872778) q[1];
sx q[1];
rz(-2.6649902) q[1];
sx q[1];
rz(2.2050819) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6694326) q[0];
sx q[0];
rz(-2.5494808) q[0];
sx q[0];
rz(-2.6955312) q[0];
rz(-pi) q[1];
rz(-2.1328743) q[2];
sx q[2];
rz(-0.33318168) q[2];
sx q[2];
rz(0.21466941) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.16441241) q[1];
sx q[1];
rz(-2.3986448) q[1];
sx q[1];
rz(-1.2765902) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1533222) q[3];
sx q[3];
rz(-1.3068145) q[3];
sx q[3];
rz(-1.9895944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5986754) q[2];
sx q[2];
rz(-0.3391372) q[2];
sx q[2];
rz(3.1266406) q[2];
rz(0.22710083) q[3];
sx q[3];
rz(-2.1588219) q[3];
sx q[3];
rz(1.2623513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99075714) q[0];
sx q[0];
rz(-0.80934722) q[0];
sx q[0];
rz(1.9833175) q[0];
rz(-1.5630209) q[1];
sx q[1];
rz(-0.75571267) q[1];
sx q[1];
rz(2.8797348) q[1];
rz(2.0832534) q[2];
sx q[2];
rz(-2.8556311) q[2];
sx q[2];
rz(-0.43559504) q[2];
rz(1.8347918) q[3];
sx q[3];
rz(-2.0859857) q[3];
sx q[3];
rz(1.2386238) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
