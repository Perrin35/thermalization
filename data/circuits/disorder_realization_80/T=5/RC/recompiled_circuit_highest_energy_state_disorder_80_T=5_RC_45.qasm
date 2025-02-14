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
rz(3.5340478) q[0];
sx q[0];
rz(3.4253149) q[0];
sx q[0];
rz(14.520159) q[0];
rz(-0.44161931) q[1];
sx q[1];
rz(-1.8752357) q[1];
sx q[1];
rz(1.1511572) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8609223) q[0];
sx q[0];
rz(-0.09083561) q[0];
sx q[0];
rz(1.8152186) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5106836) q[2];
sx q[2];
rz(-0.5782457) q[2];
sx q[2];
rz(-1.1743682) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2888803) q[1];
sx q[1];
rz(-1.4433644) q[1];
sx q[1];
rz(2.7136449) q[1];
rz(2.092716) q[3];
sx q[3];
rz(-1.3112231) q[3];
sx q[3];
rz(-2.5775016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89897951) q[2];
sx q[2];
rz(-2.2448764) q[2];
sx q[2];
rz(2.8420281) q[2];
rz(0.35563955) q[3];
sx q[3];
rz(-0.99323291) q[3];
sx q[3];
rz(-0.43660823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75293175) q[0];
sx q[0];
rz(-0.54838538) q[0];
sx q[0];
rz(0.27632982) q[0];
rz(0.07946864) q[1];
sx q[1];
rz(-2.6092968) q[1];
sx q[1];
rz(0.4963378) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4547098) q[0];
sx q[0];
rz(-1.0609896) q[0];
sx q[0];
rz(-2.8532993) q[0];
rz(-pi) q[1];
rz(-0.33832835) q[2];
sx q[2];
rz(-1.4979432) q[2];
sx q[2];
rz(-0.05073994) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6092094) q[1];
sx q[1];
rz(-2.2953798) q[1];
sx q[1];
rz(-2.0128472) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1383508) q[3];
sx q[3];
rz(-1.7326983) q[3];
sx q[3];
rz(1.1461794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3695099) q[2];
sx q[2];
rz(-1.778506) q[2];
sx q[2];
rz(-0.17862865) q[2];
rz(1.5313088) q[3];
sx q[3];
rz(-2.2154) q[3];
sx q[3];
rz(-0.76242751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4302706) q[0];
sx q[0];
rz(-2.0549759) q[0];
sx q[0];
rz(1.7311199) q[0];
rz(1.7544282) q[1];
sx q[1];
rz(-1.4924563) q[1];
sx q[1];
rz(1.656146) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1324155) q[0];
sx q[0];
rz(-0.88736358) q[0];
sx q[0];
rz(1.6130436) q[0];
rz(-pi) q[1];
rz(-2.335821) q[2];
sx q[2];
rz(-1.8745595) q[2];
sx q[2];
rz(-0.31651503) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8285813) q[1];
sx q[1];
rz(-1.269956) q[1];
sx q[1];
rz(-1.1238696) q[1];
x q[2];
rz(-1.1165421) q[3];
sx q[3];
rz(-2.4245533) q[3];
sx q[3];
rz(-1.7771135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2014655) q[2];
sx q[2];
rz(-2.1374233) q[2];
sx q[2];
rz(0.49261576) q[2];
rz(-2.7228739) q[3];
sx q[3];
rz(-1.0227572) q[3];
sx q[3];
rz(2.54971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(3.1241322) q[0];
sx q[0];
rz(-3.0191665) q[0];
sx q[0];
rz(0.19805743) q[0];
rz(-0.65912229) q[1];
sx q[1];
rz(-2.5810869) q[1];
sx q[1];
rz(-2.9445599) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2360508) q[0];
sx q[0];
rz(-1.5371548) q[0];
sx q[0];
rz(2.6108526) q[0];
rz(-2.2617499) q[2];
sx q[2];
rz(-0.36565271) q[2];
sx q[2];
rz(-1.3297538) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2391165) q[1];
sx q[1];
rz(-2.4817339) q[1];
sx q[1];
rz(-2.541594) q[1];
x q[2];
rz(-1.0258338) q[3];
sx q[3];
rz(-1.539388) q[3];
sx q[3];
rz(1.9592913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3531602) q[2];
sx q[2];
rz(-1.3015231) q[2];
sx q[2];
rz(-1.8678467) q[2];
rz(-1.5431917) q[3];
sx q[3];
rz(-1.9485995) q[3];
sx q[3];
rz(0.25632349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0565599) q[0];
sx q[0];
rz(-1.4075449) q[0];
sx q[0];
rz(-3.0388167) q[0];
rz(-3.00434) q[1];
sx q[1];
rz(-0.44932258) q[1];
sx q[1];
rz(1.7030254) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.611127) q[0];
sx q[0];
rz(-2.1411863) q[0];
sx q[0];
rz(-2.3378951) q[0];
x q[1];
rz(2.6784727) q[2];
sx q[2];
rz(-1.5296522) q[2];
sx q[2];
rz(0.1938614) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2762982) q[1];
sx q[1];
rz(-1.7182847) q[1];
sx q[1];
rz(-1.3091722) q[1];
rz(-2.4047671) q[3];
sx q[3];
rz(-2.9758436) q[3];
sx q[3];
rz(-2.152657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.30456257) q[2];
sx q[2];
rz(-1.6670767) q[2];
sx q[2];
rz(2.4414818) q[2];
rz(-2.2434798) q[3];
sx q[3];
rz(-2.5105295) q[3];
sx q[3];
rz(1.496544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9333618) q[0];
sx q[0];
rz(-2.5550186) q[0];
sx q[0];
rz(1.5775648) q[0];
rz(2.481752) q[1];
sx q[1];
rz(-0.78958646) q[1];
sx q[1];
rz(-2.1671364) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3269269) q[0];
sx q[0];
rz(-1.5788933) q[0];
sx q[0];
rz(-0.011122313) q[0];
rz(-pi) q[1];
rz(0.31807138) q[2];
sx q[2];
rz(-0.50919767) q[2];
sx q[2];
rz(-1.731002) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.240307) q[1];
sx q[1];
rz(-2.0723596) q[1];
sx q[1];
rz(1.5585414) q[1];
x q[2];
rz(1.5100689) q[3];
sx q[3];
rz(-1.0332881) q[3];
sx q[3];
rz(0.93858657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.66880354) q[2];
sx q[2];
rz(-2.4432224) q[2];
sx q[2];
rz(-2.2051956) q[2];
rz(0.89013731) q[3];
sx q[3];
rz(-1.7809159) q[3];
sx q[3];
rz(-2.9357125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5278006) q[0];
sx q[0];
rz(-0.19659909) q[0];
sx q[0];
rz(3.0278681) q[0];
rz(1.3971036) q[1];
sx q[1];
rz(-1.0662268) q[1];
sx q[1];
rz(0.024959175) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67201383) q[0];
sx q[0];
rz(-0.22525283) q[0];
sx q[0];
rz(1.8173056) q[0];
x q[1];
rz(-1.8309085) q[2];
sx q[2];
rz(-0.94365135) q[2];
sx q[2];
rz(2.8994034) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0164707) q[1];
sx q[1];
rz(-2.2844587) q[1];
sx q[1];
rz(3.0550856) q[1];
rz(1.5146144) q[3];
sx q[3];
rz(-1.3732519) q[3];
sx q[3];
rz(2.2584525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.52594841) q[2];
sx q[2];
rz(-1.4333466) q[2];
sx q[2];
rz(2.136266) q[2];
rz(-0.89056906) q[3];
sx q[3];
rz(-1.7361448) q[3];
sx q[3];
rz(-2.4864206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5263851) q[0];
sx q[0];
rz(-0.40281519) q[0];
sx q[0];
rz(2.030754) q[0];
rz(0.088518294) q[1];
sx q[1];
rz(-2.7169777) q[1];
sx q[1];
rz(1.7875338) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.272404) q[0];
sx q[0];
rz(-1.4581632) q[0];
sx q[0];
rz(2.3478481) q[0];
x q[1];
rz(1.0646787) q[2];
sx q[2];
rz(-0.30914834) q[2];
sx q[2];
rz(2.791601) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.718487) q[1];
sx q[1];
rz(-2.8837648) q[1];
sx q[1];
rz(0.30155525) q[1];
rz(-pi) q[2];
rz(-0.76137162) q[3];
sx q[3];
rz(-1.1283481) q[3];
sx q[3];
rz(-0.82514742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2985349) q[2];
sx q[2];
rz(-0.83710805) q[2];
sx q[2];
rz(2.7561772) q[2];
rz(0.75000969) q[3];
sx q[3];
rz(-2.1842726) q[3];
sx q[3];
rz(0.065940417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0176004) q[0];
sx q[0];
rz(-2.0662859) q[0];
sx q[0];
rz(-0.83220926) q[0];
rz(2.0596793) q[1];
sx q[1];
rz(-2.3257906) q[1];
sx q[1];
rz(1.5218511) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6523217) q[0];
sx q[0];
rz(-1.1025363) q[0];
sx q[0];
rz(0.44083198) q[0];
rz(-pi) q[1];
rz(-0.25427962) q[2];
sx q[2];
rz(-1.8449718) q[2];
sx q[2];
rz(2.5667896) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.79241133) q[1];
sx q[1];
rz(-2.1508625) q[1];
sx q[1];
rz(0.68934727) q[1];
rz(2.2275739) q[3];
sx q[3];
rz(-1.7631712) q[3];
sx q[3];
rz(0.094948204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1443783) q[2];
sx q[2];
rz(-1.245456) q[2];
sx q[2];
rz(-1.9027556) q[2];
rz(1.291412) q[3];
sx q[3];
rz(-1.8700799) q[3];
sx q[3];
rz(0.13135697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0111897) q[0];
sx q[0];
rz(-0.953453) q[0];
sx q[0];
rz(2.3977872) q[0];
rz(3.0939843) q[1];
sx q[1];
rz(-1.1046537) q[1];
sx q[1];
rz(2.0000752) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6358546) q[0];
sx q[0];
rz(-2.4199652) q[0];
sx q[0];
rz(-2.3139335) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.52454) q[2];
sx q[2];
rz(-1.8931555) q[2];
sx q[2];
rz(-2.6094391) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0298668) q[1];
sx q[1];
rz(-0.85241417) q[1];
sx q[1];
rz(-0.54775441) q[1];
rz(-2.4968653) q[3];
sx q[3];
rz(-1.4253741) q[3];
sx q[3];
rz(-2.4032088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6514674) q[2];
sx q[2];
rz(-0.23568025) q[2];
sx q[2];
rz(-0.30533314) q[2];
rz(-1.7134282) q[3];
sx q[3];
rz(-1.7686663) q[3];
sx q[3];
rz(1.8491687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74844985) q[0];
sx q[0];
rz(-2.9105757) q[0];
sx q[0];
rz(-2.0052746) q[0];
rz(-2.582386) q[1];
sx q[1];
rz(-1.8103841) q[1];
sx q[1];
rz(0.24787535) q[1];
rz(-2.7087173) q[2];
sx q[2];
rz(-0.71472157) q[2];
sx q[2];
rz(1.0403487) q[2];
rz(1.4356267) q[3];
sx q[3];
rz(-1.4528989) q[3];
sx q[3];
rz(-1.0624878) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
