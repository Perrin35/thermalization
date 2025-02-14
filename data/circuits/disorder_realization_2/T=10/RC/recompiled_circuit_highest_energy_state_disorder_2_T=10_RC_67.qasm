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
rz(2.958137) q[0];
sx q[0];
rz(-2.3975211) q[0];
sx q[0];
rz(-2.0896572) q[0];
rz(0.19368859) q[1];
sx q[1];
rz(-0.63313484) q[1];
sx q[1];
rz(0.57300353) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0664803) q[0];
sx q[0];
rz(-1.6329935) q[0];
sx q[0];
rz(-2.1960427) q[0];
x q[1];
rz(-0.51197211) q[2];
sx q[2];
rz(-0.7535156) q[2];
sx q[2];
rz(1.7640424) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.68202735) q[1];
sx q[1];
rz(-0.50217705) q[1];
sx q[1];
rz(-2.9952496) q[1];
x q[2];
rz(-3.0188881) q[3];
sx q[3];
rz(-0.70022455) q[3];
sx q[3];
rz(-0.52470696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.28213349) q[2];
sx q[2];
rz(-0.30974516) q[2];
sx q[2];
rz(-1.0563043) q[2];
rz(0.022196444) q[3];
sx q[3];
rz(-0.76265097) q[3];
sx q[3];
rz(-1.7215884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2317155) q[0];
sx q[0];
rz(-2.660399) q[0];
sx q[0];
rz(0.43014446) q[0];
rz(-3.0138956) q[1];
sx q[1];
rz(-1.9857429) q[1];
sx q[1];
rz(-1.4375623) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6352954) q[0];
sx q[0];
rz(-2.7397635) q[0];
sx q[0];
rz(-0.70962972) q[0];
rz(1.6280074) q[2];
sx q[2];
rz(-1.5199516) q[2];
sx q[2];
rz(1.0982996) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9641407) q[1];
sx q[1];
rz(-2.3712082) q[1];
sx q[1];
rz(-0.48078534) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3540157) q[3];
sx q[3];
rz(-2.1000897) q[3];
sx q[3];
rz(2.2640219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.11640707) q[2];
sx q[2];
rz(-1.6802639) q[2];
sx q[2];
rz(-2.6961668) q[2];
rz(-0.40924254) q[3];
sx q[3];
rz(-2.0425115) q[3];
sx q[3];
rz(-0.87944952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5889848) q[0];
sx q[0];
rz(-0.092985066) q[0];
sx q[0];
rz(0.71480042) q[0];
rz(-2.0843166) q[1];
sx q[1];
rz(-2.7209268) q[1];
sx q[1];
rz(0.32726273) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11571685) q[0];
sx q[0];
rz(-1.030632) q[0];
sx q[0];
rz(-0.19465036) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1959841) q[2];
sx q[2];
rz(-0.66788061) q[2];
sx q[2];
rz(0.57648522) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2043874) q[1];
sx q[1];
rz(-1.722285) q[1];
sx q[1];
rz(0.020843055) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9182358) q[3];
sx q[3];
rz(-1.8104324) q[3];
sx q[3];
rz(-2.6289866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1383692) q[2];
sx q[2];
rz(-0.97347632) q[2];
sx q[2];
rz(-1.5039911) q[2];
rz(-1.9231632) q[3];
sx q[3];
rz(-2.7259493) q[3];
sx q[3];
rz(-2.159582) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0686491) q[0];
sx q[0];
rz(-1.8953841) q[0];
sx q[0];
rz(-0.15596998) q[0];
rz(-2.7929557) q[1];
sx q[1];
rz(-2.5376384) q[1];
sx q[1];
rz(-1.2329996) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017576305) q[0];
sx q[0];
rz(-1.6917545) q[0];
sx q[0];
rz(-1.8328387) q[0];
rz(1.2945588) q[2];
sx q[2];
rz(-0.66670185) q[2];
sx q[2];
rz(-1.897097) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8089227) q[1];
sx q[1];
rz(-0.70007174) q[1];
sx q[1];
rz(-2.0168522) q[1];
rz(0.25147049) q[3];
sx q[3];
rz(-1.5565539) q[3];
sx q[3];
rz(-2.3492658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6984581) q[2];
sx q[2];
rz(-3.105574) q[2];
sx q[2];
rz(-1.6301463) q[2];
rz(-1.1394507) q[3];
sx q[3];
rz(-1.3316493) q[3];
sx q[3];
rz(2.1512234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
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
rz(1.3646506) q[0];
sx q[0];
rz(-2.7237837) q[0];
sx q[0];
rz(-1.1530217) q[0];
rz(-2.5979089) q[1];
sx q[1];
rz(-1.2666356) q[1];
sx q[1];
rz(-0.26184729) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1516929) q[0];
sx q[0];
rz(-1.6542871) q[0];
sx q[0];
rz(3.1079453) q[0];
rz(-pi) q[1];
rz(2.1543617) q[2];
sx q[2];
rz(-1.2886184) q[2];
sx q[2];
rz(0.63167494) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.10298577) q[1];
sx q[1];
rz(-0.9894045) q[1];
sx q[1];
rz(-1.3036779) q[1];
rz(-pi) q[2];
rz(0.7821726) q[3];
sx q[3];
rz(-1.425079) q[3];
sx q[3];
rz(-0.71469939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.61949817) q[2];
sx q[2];
rz(-0.60569373) q[2];
sx q[2];
rz(1.7363133) q[2];
rz(-0.74553472) q[3];
sx q[3];
rz(-1.2657974) q[3];
sx q[3];
rz(0.94329992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0010506823) q[0];
sx q[0];
rz(-0.11740919) q[0];
sx q[0];
rz(2.7897799) q[0];
rz(-2.70128) q[1];
sx q[1];
rz(-1.7761209) q[1];
sx q[1];
rz(2.9702759) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3765613) q[0];
sx q[0];
rz(-0.91463551) q[0];
sx q[0];
rz(-1.879346) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2662884) q[2];
sx q[2];
rz(-2.300548) q[2];
sx q[2];
rz(-2.8685121) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3503987) q[1];
sx q[1];
rz(-1.5499252) q[1];
sx q[1];
rz(2.4430165) q[1];
rz(0.76254179) q[3];
sx q[3];
rz(-2.2665215) q[3];
sx q[3];
rz(2.1660337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.064284023) q[2];
sx q[2];
rz(-1.372154) q[2];
sx q[2];
rz(-2.1997931) q[2];
rz(1.1549548) q[3];
sx q[3];
rz(-0.20808163) q[3];
sx q[3];
rz(-1.340516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6963541) q[0];
sx q[0];
rz(-1.1055163) q[0];
sx q[0];
rz(3.136694) q[0];
rz(-0.44662961) q[1];
sx q[1];
rz(-0.59372562) q[1];
sx q[1];
rz(-0.68797025) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1991581) q[0];
sx q[0];
rz(-0.73478886) q[0];
sx q[0];
rz(2.8834729) q[0];
x q[1];
rz(2.6737988) q[2];
sx q[2];
rz(-2.4681598) q[2];
sx q[2];
rz(1.837567) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6354949) q[1];
sx q[1];
rz(-0.48759547) q[1];
sx q[1];
rz(2.0388842) q[1];
x q[2];
rz(2.2673685) q[3];
sx q[3];
rz(-2.1316574) q[3];
sx q[3];
rz(-1.107533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.61116162) q[2];
sx q[2];
rz(-0.40803424) q[2];
sx q[2];
rz(0.92998695) q[2];
rz(1.2215349) q[3];
sx q[3];
rz(-2.0285716) q[3];
sx q[3];
rz(0.59337029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68194836) q[0];
sx q[0];
rz(-1.821803) q[0];
sx q[0];
rz(-0.090106877) q[0];
rz(1.4247165) q[1];
sx q[1];
rz(-0.63980353) q[1];
sx q[1];
rz(-2.7065014) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97532192) q[0];
sx q[0];
rz(-1.4251801) q[0];
sx q[0];
rz(2.6377489) q[0];
x q[1];
rz(-2.4736011) q[2];
sx q[2];
rz(-1.4282303) q[2];
sx q[2];
rz(3.0270456) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.69999483) q[1];
sx q[1];
rz(-0.95788237) q[1];
sx q[1];
rz(-0.15721486) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5716721) q[3];
sx q[3];
rz(-1.9951576) q[3];
sx q[3];
rz(-0.76155886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6626176) q[2];
sx q[2];
rz(-2.0765897) q[2];
sx q[2];
rz(-2.5153861) q[2];
rz(-2.9446757) q[3];
sx q[3];
rz(-2.1508689) q[3];
sx q[3];
rz(-1.3885434) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58186746) q[0];
sx q[0];
rz(-1.6904866) q[0];
sx q[0];
rz(0.054280601) q[0];
rz(2.718603) q[1];
sx q[1];
rz(-1.3754247) q[1];
sx q[1];
rz(-1.8064226) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2535219) q[0];
sx q[0];
rz(-1.0253064) q[0];
sx q[0];
rz(2.2175199) q[0];
rz(-2.7846856) q[2];
sx q[2];
rz(-1.0400598) q[2];
sx q[2];
rz(-2.4414947) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.74354913) q[1];
sx q[1];
rz(-2.365592) q[1];
sx q[1];
rz(2.3659124) q[1];
rz(-pi) q[2];
rz(-1.7087206) q[3];
sx q[3];
rz(-1.63878) q[3];
sx q[3];
rz(1.8658416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6173031) q[2];
sx q[2];
rz(-1.6733988) q[2];
sx q[2];
rz(-2.7900901) q[2];
rz(-1.7678123) q[3];
sx q[3];
rz(-2.6280554) q[3];
sx q[3];
rz(-0.26941776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7518625) q[0];
sx q[0];
rz(-2.8808012) q[0];
sx q[0];
rz(-2.3824298) q[0];
rz(1.0836541) q[1];
sx q[1];
rz(-1.5307129) q[1];
sx q[1];
rz(2.0738475) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5006951) q[0];
sx q[0];
rz(-1.0984549) q[0];
sx q[0];
rz(0.91820902) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2933838) q[2];
sx q[2];
rz(-1.2870064) q[2];
sx q[2];
rz(1.0671771) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2764923) q[1];
sx q[1];
rz(-0.70521627) q[1];
sx q[1];
rz(2.221285) q[1];
rz(-pi) q[2];
rz(-2.2678947) q[3];
sx q[3];
rz(-1.3109968) q[3];
sx q[3];
rz(1.3316621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7119673) q[2];
sx q[2];
rz(-1.9316614) q[2];
sx q[2];
rz(2.9313226) q[2];
rz(-0.78768864) q[3];
sx q[3];
rz(-1.7588047) q[3];
sx q[3];
rz(2.6113966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44337153) q[0];
sx q[0];
rz(-0.5589232) q[0];
sx q[0];
rz(1.9236175) q[0];
rz(-1.5086077) q[1];
sx q[1];
rz(-1.8764381) q[1];
sx q[1];
rz(1.1988342) q[1];
rz(2.0532578) q[2];
sx q[2];
rz(-2.3934622) q[2];
sx q[2];
rz(-0.29665034) q[2];
rz(-0.15565025) q[3];
sx q[3];
rz(-1.8425103) q[3];
sx q[3];
rz(1.507148) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
