OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4189897) q[0];
sx q[0];
rz(-0.64282066) q[0];
sx q[0];
rz(1.2195725) q[0];
rz(-2.3063083) q[1];
sx q[1];
rz(-1.78479) q[1];
sx q[1];
rz(-0.9019444) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.565958) q[0];
sx q[0];
rz(-0.96372858) q[0];
sx q[0];
rz(-1.3835039) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44166126) q[2];
sx q[2];
rz(-2.7251789) q[2];
sx q[2];
rz(3.1232782) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8512501) q[1];
sx q[1];
rz(-0.34091967) q[1];
sx q[1];
rz(1.0268289) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2252623) q[3];
sx q[3];
rz(-0.95312762) q[3];
sx q[3];
rz(-0.82181069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.96282643) q[2];
sx q[2];
rz(-1.0054532) q[2];
sx q[2];
rz(0.82628769) q[2];
rz(1.4055584) q[3];
sx q[3];
rz(-2.8345351) q[3];
sx q[3];
rz(2.3981986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58251441) q[0];
sx q[0];
rz(-2.7358416) q[0];
sx q[0];
rz(-1.7412809) q[0];
rz(-2.0179613) q[1];
sx q[1];
rz(-2.7476937) q[1];
sx q[1];
rz(1.0981015) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0535677) q[0];
sx q[0];
rz(-2.9505749) q[0];
sx q[0];
rz(0.6121401) q[0];
x q[1];
rz(-1.9015354) q[2];
sx q[2];
rz(-1.2501226) q[2];
sx q[2];
rz(-1.5000686) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.84867618) q[1];
sx q[1];
rz(-1.4956988) q[1];
sx q[1];
rz(-3.0203222) q[1];
rz(-1.9722749) q[3];
sx q[3];
rz(-2.7352834) q[3];
sx q[3];
rz(-0.10029785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0850247) q[2];
sx q[2];
rz(-2.8242064) q[2];
sx q[2];
rz(-1.8918096) q[2];
rz(-1.7791087) q[3];
sx q[3];
rz(-1.0008078) q[3];
sx q[3];
rz(1.6574297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5054841) q[0];
sx q[0];
rz(-0.69378575) q[0];
sx q[0];
rz(2.8662477) q[0];
rz(0.45922008) q[1];
sx q[1];
rz(-0.95454916) q[1];
sx q[1];
rz(0.053344639) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33063525) q[0];
sx q[0];
rz(-0.22718469) q[0];
sx q[0];
rz(0.35463984) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0317475) q[2];
sx q[2];
rz(-1.58733) q[2];
sx q[2];
rz(1.2993882) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5819323) q[1];
sx q[1];
rz(-1.6114863) q[1];
sx q[1];
rz(2.5896124) q[1];
x q[2];
rz(-2.1249849) q[3];
sx q[3];
rz(-0.38180581) q[3];
sx q[3];
rz(-3.0757381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.82103819) q[2];
sx q[2];
rz(-0.36131636) q[2];
sx q[2];
rz(1.5041171) q[2];
rz(1.5004246) q[3];
sx q[3];
rz(-1.0502366) q[3];
sx q[3];
rz(-2.8659081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42030537) q[0];
sx q[0];
rz(-2.5581701) q[0];
sx q[0];
rz(-2.6570008) q[0];
rz(1.8991607) q[1];
sx q[1];
rz(-2.395605) q[1];
sx q[1];
rz(-0.34908435) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42345966) q[0];
sx q[0];
rz(-1.0225931) q[0];
sx q[0];
rz(-3.0915652) q[0];
rz(-pi) q[1];
rz(1.3626839) q[2];
sx q[2];
rz(-0.35018626) q[2];
sx q[2];
rz(1.1454358) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.71098498) q[1];
sx q[1];
rz(-2.8423956) q[1];
sx q[1];
rz(-2.4483212) q[1];
x q[2];
rz(0.99557568) q[3];
sx q[3];
rz(-2.4565426) q[3];
sx q[3];
rz(0.12888651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.35859534) q[2];
sx q[2];
rz(-1.9390257) q[2];
sx q[2];
rz(-0.45004582) q[2];
rz(-0.83886823) q[3];
sx q[3];
rz(-1.5939555) q[3];
sx q[3];
rz(-2.94256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6435476) q[0];
sx q[0];
rz(-1.4692551) q[0];
sx q[0];
rz(0.099844649) q[0];
rz(1.3205344) q[1];
sx q[1];
rz(-2.1456199) q[1];
sx q[1];
rz(0.60755306) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4488825) q[0];
sx q[0];
rz(-1.4878055) q[0];
sx q[0];
rz(1.4957499) q[0];
x q[1];
rz(-0.25600146) q[2];
sx q[2];
rz(-2.5814272) q[2];
sx q[2];
rz(1.7988811) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.4974571) q[1];
sx q[1];
rz(-2.1408399) q[1];
sx q[1];
rz(-1.312478) q[1];
rz(-pi) q[2];
rz(-1.2149548) q[3];
sx q[3];
rz(-0.27493335) q[3];
sx q[3];
rz(-0.41809088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1943835) q[2];
sx q[2];
rz(-0.76346976) q[2];
sx q[2];
rz(-0.53059951) q[2];
rz(0.44140205) q[3];
sx q[3];
rz(-0.97359052) q[3];
sx q[3];
rz(0.94943625) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
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
rz(-2.0206385) q[1];
sx q[1];
rz(-1.3425739) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.955637) q[0];
sx q[0];
rz(-1.370472) q[0];
sx q[0];
rz(1.093051) q[0];
x q[1];
rz(2.3903177) q[2];
sx q[2];
rz(-0.57005586) q[2];
sx q[2];
rz(1.7715724) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8042829) q[1];
sx q[1];
rz(-1.8071257) q[1];
sx q[1];
rz(1.9419036) q[1];
x q[2];
rz(0.66304147) q[3];
sx q[3];
rz(-1.6812857) q[3];
sx q[3];
rz(-2.7821531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.67419702) q[2];
sx q[2];
rz(-2.8574222) q[2];
sx q[2];
rz(2.6944842) q[2];
rz(-2.1111264) q[3];
sx q[3];
rz(-1.5117398) q[3];
sx q[3];
rz(-0.38465056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(2.7717188) q[0];
sx q[0];
rz(-1.2971224) q[0];
sx q[0];
rz(2.7287927) q[0];
rz(-1.075047) q[1];
sx q[1];
rz(-2.0037035) q[1];
sx q[1];
rz(-2.3379751) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61148723) q[0];
sx q[0];
rz(-0.75615935) q[0];
sx q[0];
rz(-2.3458781) q[0];
rz(-2.2140131) q[2];
sx q[2];
rz(-1.1021492) q[2];
sx q[2];
rz(1.8806632) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6576801) q[1];
sx q[1];
rz(-1.5849304) q[1];
sx q[1];
rz(-0.78426984) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7429306) q[3];
sx q[3];
rz(-0.19426647) q[3];
sx q[3];
rz(-2.8630321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.86218086) q[2];
sx q[2];
rz(-1.2181686) q[2];
sx q[2];
rz(-1.6406406) q[2];
rz(2.5189279) q[3];
sx q[3];
rz(-0.77087918) q[3];
sx q[3];
rz(1.4890495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72614661) q[0];
sx q[0];
rz(-1.5203238) q[0];
sx q[0];
rz(0.55794445) q[0];
rz(2.5803512) q[1];
sx q[1];
rz(-1.3713501) q[1];
sx q[1];
rz(1.4161313) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.95074) q[0];
sx q[0];
rz(-1.8537562) q[0];
sx q[0];
rz(-1.1038194) q[0];
rz(-pi) q[1];
rz(0.92500706) q[2];
sx q[2];
rz(-2.6753798) q[2];
sx q[2];
rz(3.1104308) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3741403) q[1];
sx q[1];
rz(-1.6686286) q[1];
sx q[1];
rz(1.2409452) q[1];
rz(-pi) q[2];
rz(1.3660356) q[3];
sx q[3];
rz(-0.8258709) q[3];
sx q[3];
rz(1.5958169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.8884362) q[2];
sx q[2];
rz(-2.3641868) q[2];
sx q[2];
rz(2.1412795) q[2];
rz(-0.12217626) q[3];
sx q[3];
rz(-1.6413942) q[3];
sx q[3];
rz(2.9848671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67824739) q[0];
sx q[0];
rz(-1.954701) q[0];
sx q[0];
rz(0.004322411) q[0];
rz(0.16383544) q[1];
sx q[1];
rz(-2.7163353) q[1];
sx q[1];
rz(1.77553) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8238433) q[0];
sx q[0];
rz(-1.5986406) q[0];
sx q[0];
rz(-2.2180024) q[0];
x q[1];
rz(0.37213426) q[2];
sx q[2];
rz(-1.8610125) q[2];
sx q[2];
rz(1.2843101) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9802347) q[1];
sx q[1];
rz(-1.606985) q[1];
sx q[1];
rz(-0.25118942) q[1];
x q[2];
rz(-2.2774599) q[3];
sx q[3];
rz(-0.88418761) q[3];
sx q[3];
rz(-0.1089801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1395448) q[2];
sx q[2];
rz(-0.9674955) q[2];
sx q[2];
rz(2.5223993) q[2];
rz(2.5891417) q[3];
sx q[3];
rz(-1.3055472) q[3];
sx q[3];
rz(-2.9964871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3510975) q[0];
sx q[0];
rz(-2.9572697) q[0];
sx q[0];
rz(2.6628394) q[0];
rz(1.9888196) q[1];
sx q[1];
rz(-2.8515127) q[1];
sx q[1];
rz(2.1943888) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9626434) q[0];
sx q[0];
rz(-1.6048204) q[0];
sx q[0];
rz(-1.3860788) q[0];
x q[1];
rz(1.2113153) q[2];
sx q[2];
rz(-0.61865846) q[2];
sx q[2];
rz(-1.930069) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.80120197) q[1];
sx q[1];
rz(-1.2653192) q[1];
sx q[1];
rz(-2.3190772) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5371505) q[3];
sx q[3];
rz(-2.7210208) q[3];
sx q[3];
rz(-0.67422359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3843711) q[2];
sx q[2];
rz(-0.20558509) q[2];
sx q[2];
rz(0.26665404) q[2];
rz(1.0673149) q[3];
sx q[3];
rz(-1.9284733) q[3];
sx q[3];
rz(2.1783569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08854475) q[0];
sx q[0];
rz(-1.5374669) q[0];
sx q[0];
rz(1.5446825) q[0];
rz(-2.0883941) q[1];
sx q[1];
rz(-1.9269301) q[1];
sx q[1];
rz(-2.3763837) q[1];
rz(-3.13022) q[2];
sx q[2];
rz(-0.24200242) q[2];
sx q[2];
rz(0.92462362) q[2];
rz(1.7942747) q[3];
sx q[3];
rz(-1.0031932) q[3];
sx q[3];
rz(-0.72465988) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
