OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.27622142) q[0];
sx q[0];
rz(5.4260317) q[0];
sx q[0];
rz(9.5572588) q[0];
rz(-2.8744856) q[1];
sx q[1];
rz(-2.5565956) q[1];
sx q[1];
rz(-2.4490228) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8938457) q[0];
sx q[0];
rz(-0.5286628) q[0];
sx q[0];
rz(1.7696487) q[0];
rz(2.4970826) q[2];
sx q[2];
rz(-1.95382) q[2];
sx q[2];
rz(0.884998) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.3092279) q[1];
sx q[1];
rz(-1.4432866) q[1];
sx q[1];
rz(3.068919) q[1];
rz(-pi) q[2];
rz(2.0444144) q[3];
sx q[3];
rz(-1.7930822) q[3];
sx q[3];
rz(1.3003295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1074368) q[2];
sx q[2];
rz(-0.52705708) q[2];
sx q[2];
rz(-1.6050603) q[2];
rz(-1.5213373) q[3];
sx q[3];
rz(-1.4884357) q[3];
sx q[3];
rz(3.1055514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81543106) q[0];
sx q[0];
rz(-2.8679929) q[0];
sx q[0];
rz(-1.2492299) q[0];
rz(-2.5800887) q[1];
sx q[1];
rz(-2.3655472) q[1];
sx q[1];
rz(-0.5805648) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8037572) q[0];
sx q[0];
rz(-2.8135186) q[0];
sx q[0];
rz(0.34123811) q[0];
rz(-pi) q[1];
rz(-2.850769) q[2];
sx q[2];
rz(-0.65867701) q[2];
sx q[2];
rz(0.80292279) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.21743449) q[1];
sx q[1];
rz(-2.1719451) q[1];
sx q[1];
rz(-0.8767599) q[1];
rz(-1.1167691) q[3];
sx q[3];
rz(-1.2516216) q[3];
sx q[3];
rz(1.9516731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7291752) q[2];
sx q[2];
rz(-2.9340332) q[2];
sx q[2];
rz(0.87835971) q[2];
rz(-0.39204028) q[3];
sx q[3];
rz(-1.6974028) q[3];
sx q[3];
rz(-2.5382606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47675258) q[0];
sx q[0];
rz(-2.219253) q[0];
sx q[0];
rz(-2.3679249) q[0];
rz(0.0013008612) q[1];
sx q[1];
rz(-1.5258077) q[1];
sx q[1];
rz(-3.1087648) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5531909) q[0];
sx q[0];
rz(-1.8817888) q[0];
sx q[0];
rz(-1.8131959) q[0];
rz(-pi) q[1];
rz(-0.29067729) q[2];
sx q[2];
rz(-2.103984) q[2];
sx q[2];
rz(0.79613396) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9877316) q[1];
sx q[1];
rz(-0.78480936) q[1];
sx q[1];
rz(-1.8042817) q[1];
rz(-pi) q[2];
rz(-2.7955416) q[3];
sx q[3];
rz(-1.0533353) q[3];
sx q[3];
rz(0.82511653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7971928) q[2];
sx q[2];
rz(-2.1160647) q[2];
sx q[2];
rz(0.27734217) q[2];
rz(2.7456361) q[3];
sx q[3];
rz(-1.5405416) q[3];
sx q[3];
rz(-2.4424281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4375777) q[0];
sx q[0];
rz(-0.33575785) q[0];
sx q[0];
rz(2.8787956) q[0];
rz(2.908005) q[1];
sx q[1];
rz(-0.83507744) q[1];
sx q[1];
rz(-2.3707726) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0924776) q[0];
sx q[0];
rz(-0.025509398) q[0];
sx q[0];
rz(2.269948) q[0];
rz(-pi) q[1];
rz(-1.0806482) q[2];
sx q[2];
rz(-1.4785826) q[2];
sx q[2];
rz(-2.8871418) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.091688823) q[1];
sx q[1];
rz(-2.4871475) q[1];
sx q[1];
rz(-0.62033886) q[1];
rz(1.3464438) q[3];
sx q[3];
rz(-1.3849349) q[3];
sx q[3];
rz(0.88691521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9757441) q[2];
sx q[2];
rz(-1.5074915) q[2];
sx q[2];
rz(-0.7129933) q[2];
rz(-1.0130079) q[3];
sx q[3];
rz(-0.37390798) q[3];
sx q[3];
rz(-2.1876984) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35904303) q[0];
sx q[0];
rz(-2.0694216) q[0];
sx q[0];
rz(1.4404526) q[0];
rz(0.094093181) q[1];
sx q[1];
rz(-2.4021939) q[1];
sx q[1];
rz(-0.17000155) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.819034) q[0];
sx q[0];
rz(-1.337262) q[0];
sx q[0];
rz(0.76604953) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.22959392) q[2];
sx q[2];
rz(-2.0320971) q[2];
sx q[2];
rz(-2.4597944) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4449094) q[1];
sx q[1];
rz(-2.141436) q[1];
sx q[1];
rz(-2.6615104) q[1];
x q[2];
rz(0.87336297) q[3];
sx q[3];
rz(-2.1057099) q[3];
sx q[3];
rz(1.3267335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1468982) q[2];
sx q[2];
rz(-2.6165104) q[2];
sx q[2];
rz(-1.4040995) q[2];
rz(-1.5935625) q[3];
sx q[3];
rz(-2.352495) q[3];
sx q[3];
rz(-1.7061957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65574044) q[0];
sx q[0];
rz(-1.9952554) q[0];
sx q[0];
rz(-2.6830542) q[0];
rz(0.25587747) q[1];
sx q[1];
rz(-1.8829388) q[1];
sx q[1];
rz(-2.4564254) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.348649) q[0];
sx q[0];
rz(-0.92175882) q[0];
sx q[0];
rz(0.087260212) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7148758) q[2];
sx q[2];
rz(-0.79877582) q[2];
sx q[2];
rz(0.069552334) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.398715) q[1];
sx q[1];
rz(-1.0814953) q[1];
sx q[1];
rz(-0.87280416) q[1];
x q[2];
rz(-0.98486793) q[3];
sx q[3];
rz(-1.586048) q[3];
sx q[3];
rz(-2.0963984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.3926065) q[2];
sx q[2];
rz(-0.41967732) q[2];
sx q[2];
rz(1.4292599) q[2];
rz(-1.0990934) q[3];
sx q[3];
rz(-0.50656879) q[3];
sx q[3];
rz(2.9523622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(3.1290865) q[0];
sx q[0];
rz(-1.5456454) q[0];
sx q[0];
rz(0.72934735) q[0];
rz(-0.29306456) q[1];
sx q[1];
rz(-0.23935071) q[1];
sx q[1];
rz(1.9940631) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4097737) q[0];
sx q[0];
rz(-1.6454576) q[0];
sx q[0];
rz(0.41418196) q[0];
x q[1];
rz(0.36979923) q[2];
sx q[2];
rz(-0.97019201) q[2];
sx q[2];
rz(-0.72593867) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.16312576) q[1];
sx q[1];
rz(-1.6118057) q[1];
sx q[1];
rz(-0.48467111) q[1];
rz(-2.494874) q[3];
sx q[3];
rz(-0.59520703) q[3];
sx q[3];
rz(-1.3884384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3532233) q[2];
sx q[2];
rz(-1.0228144) q[2];
sx q[2];
rz(-2.3925171) q[2];
rz(-0.64368147) q[3];
sx q[3];
rz(-2.1285074) q[3];
sx q[3];
rz(-0.914004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1483243) q[0];
sx q[0];
rz(-2.0355621) q[0];
sx q[0];
rz(-0.83129445) q[0];
rz(1.7656901) q[1];
sx q[1];
rz(-0.81326905) q[1];
sx q[1];
rz(2.7430699) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7300028) q[0];
sx q[0];
rz(-1.411479) q[0];
sx q[0];
rz(-1.0986064) q[0];
rz(-pi) q[1];
rz(-2.987791) q[2];
sx q[2];
rz(-2.5349344) q[2];
sx q[2];
rz(2.3369044) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.41259137) q[1];
sx q[1];
rz(-1.593643) q[1];
sx q[1];
rz(2.1680135) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6284292) q[3];
sx q[3];
rz(-2.0675821) q[3];
sx q[3];
rz(2.7548807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1901671) q[2];
sx q[2];
rz(-1.7931033) q[2];
sx q[2];
rz(1.8438967) q[2];
rz(-1.1249582) q[3];
sx q[3];
rz(-1.8816032) q[3];
sx q[3];
rz(-2.4979533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1289537) q[0];
sx q[0];
rz(-2.5057827) q[0];
sx q[0];
rz(1.3893611) q[0];
rz(-1.5147491) q[1];
sx q[1];
rz(-1.6747968) q[1];
sx q[1];
rz(1.0983889) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3745981) q[0];
sx q[0];
rz(-0.54531389) q[0];
sx q[0];
rz(1.1101515) q[0];
rz(-pi) q[1];
rz(2.3214925) q[2];
sx q[2];
rz(-1.9165877) q[2];
sx q[2];
rz(2.1095914) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.92022773) q[1];
sx q[1];
rz(-0.5792633) q[1];
sx q[1];
rz(3.1406162) q[1];
rz(-1.3994201) q[3];
sx q[3];
rz(-1.8212574) q[3];
sx q[3];
rz(0.40303883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8998469) q[2];
sx q[2];
rz(-1.2925623) q[2];
sx q[2];
rz(2.6220654) q[2];
rz(1.8064921) q[3];
sx q[3];
rz(-2.3050008) q[3];
sx q[3];
rz(1.3174723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7463503) q[0];
sx q[0];
rz(-2.0210176) q[0];
sx q[0];
rz(-2.675132) q[0];
rz(-0.17164224) q[1];
sx q[1];
rz(-1.2152834) q[1];
sx q[1];
rz(-2.5126273) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91613149) q[0];
sx q[0];
rz(-0.3779419) q[0];
sx q[0];
rz(1.7037237) q[0];
rz(-pi) q[1];
rz(1.2538818) q[2];
sx q[2];
rz(-0.66714087) q[2];
sx q[2];
rz(-2.4184879) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2168857) q[1];
sx q[1];
rz(-2.6994884) q[1];
sx q[1];
rz(-3.1367723) q[1];
rz(-2.7528796) q[3];
sx q[3];
rz(-0.75838381) q[3];
sx q[3];
rz(0.10314108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.24370596) q[2];
sx q[2];
rz(-2.0508998) q[2];
sx q[2];
rz(2.0824599) q[2];
rz(0.6774261) q[3];
sx q[3];
rz(-2.1493561) q[3];
sx q[3];
rz(-0.89390755) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8582936) q[0];
sx q[0];
rz(-2.0170006) q[0];
sx q[0];
rz(2.2977805) q[0];
rz(2.8876866) q[1];
sx q[1];
rz(-2.0575247) q[1];
sx q[1];
rz(2.5622096) q[1];
rz(0.049747808) q[2];
sx q[2];
rz(-0.72409734) q[2];
sx q[2];
rz(-3.004336) q[2];
rz(3.0729978) q[3];
sx q[3];
rz(-2.1566236) q[3];
sx q[3];
rz(-1.9038283) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
