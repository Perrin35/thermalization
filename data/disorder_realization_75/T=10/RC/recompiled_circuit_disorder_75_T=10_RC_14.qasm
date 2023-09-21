OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5443213) q[0];
sx q[0];
rz(3.5452329) q[0];
sx q[0];
rz(9.7950254) q[0];
rz(-2.9397842) q[1];
sx q[1];
rz(-1.8887853) q[1];
sx q[1];
rz(-1.7226146) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0231853) q[0];
sx q[0];
rz(-1.4010795) q[0];
sx q[0];
rz(-0.24766185) q[0];
x q[1];
rz(-0.51214829) q[2];
sx q[2];
rz(-1.3379659) q[2];
sx q[2];
rz(-1.6793959) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.943489) q[1];
sx q[1];
rz(-1.6232821) q[1];
sx q[1];
rz(2.1211038) q[1];
x q[2];
rz(0.3317823) q[3];
sx q[3];
rz(-2.3506769) q[3];
sx q[3];
rz(0.055981759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3966763) q[2];
sx q[2];
rz(-1.6884721) q[2];
sx q[2];
rz(-0.35787004) q[2];
rz(-0.19168028) q[3];
sx q[3];
rz(-2.7087757) q[3];
sx q[3];
rz(-2.5884957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1382004) q[0];
sx q[0];
rz(-1.1509742) q[0];
sx q[0];
rz(-0.068280846) q[0];
rz(1.0066907) q[1];
sx q[1];
rz(-0.12193646) q[1];
sx q[1];
rz(-0.65111792) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3119445) q[0];
sx q[0];
rz(-0.94046578) q[0];
sx q[0];
rz(-1.3208566) q[0];
x q[1];
rz(-3.100349) q[2];
sx q[2];
rz(-2.6154499) q[2];
sx q[2];
rz(2.0522842) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5755641) q[1];
sx q[1];
rz(-2.063942) q[1];
sx q[1];
rz(-0.50599392) q[1];
rz(-1.7884368) q[3];
sx q[3];
rz(-1.0895035) q[3];
sx q[3];
rz(-2.7841115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6887001) q[2];
sx q[2];
rz(-0.15658997) q[2];
sx q[2];
rz(-2.2382656) q[2];
rz(0.7545169) q[3];
sx q[3];
rz(-1.4947596) q[3];
sx q[3];
rz(0.27404684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2925401) q[0];
sx q[0];
rz(-0.64340574) q[0];
sx q[0];
rz(-2.6480411) q[0];
rz(-1.6429398) q[1];
sx q[1];
rz(-2.7354) q[1];
sx q[1];
rz(1.0292056) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2632336) q[0];
sx q[0];
rz(-1.2326816) q[0];
sx q[0];
rz(-0.46237207) q[0];
x q[1];
rz(-2.2649293) q[2];
sx q[2];
rz(-0.52774094) q[2];
sx q[2];
rz(-1.1610247) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.610565) q[1];
sx q[1];
rz(-0.37279168) q[1];
sx q[1];
rz(2.5337266) q[1];
x q[2];
rz(3.0164099) q[3];
sx q[3];
rz(-2.1117961) q[3];
sx q[3];
rz(2.560905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9323953) q[2];
sx q[2];
rz(-0.63964996) q[2];
sx q[2];
rz(-1.2970682) q[2];
rz(-2.5629937) q[3];
sx q[3];
rz(-1.9208627) q[3];
sx q[3];
rz(0.67563081) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9349174) q[0];
sx q[0];
rz(-0.6466372) q[0];
sx q[0];
rz(-0.40859616) q[0];
rz(-1.4273377) q[1];
sx q[1];
rz(-1.4988377) q[1];
sx q[1];
rz(0.88159195) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13215412) q[0];
sx q[0];
rz(-2.3883925) q[0];
sx q[0];
rz(2.8892345) q[0];
rz(-pi) q[1];
rz(-0.84076968) q[2];
sx q[2];
rz(-0.86849125) q[2];
sx q[2];
rz(-2.4384987) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9095163) q[1];
sx q[1];
rz(-1.0788003) q[1];
sx q[1];
rz(-0.75922774) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9922585) q[3];
sx q[3];
rz(-0.42156005) q[3];
sx q[3];
rz(1.644852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0980229) q[2];
sx q[2];
rz(-1.6514401) q[2];
sx q[2];
rz(2.0289452) q[2];
rz(2.5860795) q[3];
sx q[3];
rz(-1.7881309) q[3];
sx q[3];
rz(1.5295193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9920138) q[0];
sx q[0];
rz(-1.9861789) q[0];
sx q[0];
rz(1.3647112) q[0];
rz(-0.31750202) q[1];
sx q[1];
rz(-2.1798539) q[1];
sx q[1];
rz(0.11725765) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2909572) q[0];
sx q[0];
rz(-2.1521749) q[0];
sx q[0];
rz(-1.4145538) q[0];
rz(-pi) q[1];
rz(2.8357382) q[2];
sx q[2];
rz(-0.62539414) q[2];
sx q[2];
rz(1.2139699) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9980364) q[1];
sx q[1];
rz(-1.7245502) q[1];
sx q[1];
rz(2.9424332) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1579671) q[3];
sx q[3];
rz(-1.5067325) q[3];
sx q[3];
rz(-2.4487045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.34510288) q[2];
sx q[2];
rz(-1.779665) q[2];
sx q[2];
rz(1.3797181) q[2];
rz(-1.951925) q[3];
sx q[3];
rz(-0.15933557) q[3];
sx q[3];
rz(0.074507944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9538486) q[0];
sx q[0];
rz(-1.501361) q[0];
sx q[0];
rz(0.70621079) q[0];
rz(-1.1114936) q[1];
sx q[1];
rz(-0.62875426) q[1];
sx q[1];
rz(0.10791735) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7100922) q[0];
sx q[0];
rz(-1.8653231) q[0];
sx q[0];
rz(0.72580238) q[0];
rz(2.8651587) q[2];
sx q[2];
rz(-3.099953) q[2];
sx q[2];
rz(2.4561938) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.38938746) q[1];
sx q[1];
rz(-2.2332193) q[1];
sx q[1];
rz(-1.9325158) q[1];
rz(-pi) q[2];
rz(2.2599254) q[3];
sx q[3];
rz(-1.1881184) q[3];
sx q[3];
rz(1.4985639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.32249054) q[2];
sx q[2];
rz(-0.97444797) q[2];
sx q[2];
rz(-2.0325913) q[2];
rz(1.8479944) q[3];
sx q[3];
rz(-1.3556017) q[3];
sx q[3];
rz(-0.10722815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.368822) q[0];
sx q[0];
rz(-1.6596376) q[0];
sx q[0];
rz(-3.1266881) q[0];
rz(-2.7203454) q[1];
sx q[1];
rz(-1.0528456) q[1];
sx q[1];
rz(-0.79963911) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0945471) q[0];
sx q[0];
rz(-2.4106815) q[0];
sx q[0];
rz(-2.7546309) q[0];
rz(-pi) q[1];
rz(-1.2866576) q[2];
sx q[2];
rz(-1.26789) q[2];
sx q[2];
rz(2.6877833) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.423973) q[1];
sx q[1];
rz(-1.7160001) q[1];
sx q[1];
rz(-1.5145626) q[1];
rz(1.7855438) q[3];
sx q[3];
rz(-2.5657006) q[3];
sx q[3];
rz(1.4240595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5853167) q[2];
sx q[2];
rz(-0.63627807) q[2];
sx q[2];
rz(-2.2195393) q[2];
rz(-1.8317892) q[3];
sx q[3];
rz(-1.9366879) q[3];
sx q[3];
rz(-0.95782763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17469445) q[0];
sx q[0];
rz(-2.0860724) q[0];
sx q[0];
rz(-2.877537) q[0];
rz(-1.4008201) q[1];
sx q[1];
rz(-1.4512647) q[1];
sx q[1];
rz(0.31731269) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70586328) q[0];
sx q[0];
rz(-1.1840491) q[0];
sx q[0];
rz(-0.078588967) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1605026) q[2];
sx q[2];
rz(-1.2264894) q[2];
sx q[2];
rz(-2.3366994) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4999769) q[1];
sx q[1];
rz(-2.9950905) q[1];
sx q[1];
rz(0.67540692) q[1];
x q[2];
rz(-0.21934261) q[3];
sx q[3];
rz(-2.4845124) q[3];
sx q[3];
rz(2.5430162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.43549609) q[2];
sx q[2];
rz(-2.2075682) q[2];
sx q[2];
rz(-1.502011) q[2];
rz(0.25029415) q[3];
sx q[3];
rz(-1.4151662) q[3];
sx q[3];
rz(1.1423473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2450927) q[0];
sx q[0];
rz(-0.30830202) q[0];
sx q[0];
rz(-1.7171575) q[0];
rz(0.4793438) q[1];
sx q[1];
rz(-1.665325) q[1];
sx q[1];
rz(-0.11553484) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4452267) q[0];
sx q[0];
rz(-1.0123024) q[0];
sx q[0];
rz(-0.88748705) q[0];
rz(-pi) q[1];
rz(-2.8250474) q[2];
sx q[2];
rz(-0.74005055) q[2];
sx q[2];
rz(-2.9209602) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2781196) q[1];
sx q[1];
rz(-1.9895456) q[1];
sx q[1];
rz(0.31660415) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69620903) q[3];
sx q[3];
rz(-1.6632102) q[3];
sx q[3];
rz(-2.948451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.66403786) q[2];
sx q[2];
rz(-1.1096191) q[2];
sx q[2];
rz(-1.7158562) q[2];
rz(1.4005631) q[3];
sx q[3];
rz(-1.0481342) q[3];
sx q[3];
rz(-0.28276309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35587674) q[0];
sx q[0];
rz(-0.73369217) q[0];
sx q[0];
rz(2.8724331) q[0];
rz(1.0956988) q[1];
sx q[1];
rz(-2.2311189) q[1];
sx q[1];
rz(1.3795308) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3740765) q[0];
sx q[0];
rz(-0.19664581) q[0];
sx q[0];
rz(-1.0023414) q[0];
rz(-1.7494406) q[2];
sx q[2];
rz(-2.7182455) q[2];
sx q[2];
rz(1.64738) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.97314944) q[1];
sx q[1];
rz(-1.397555) q[1];
sx q[1];
rz(-0.030862191) q[1];
x q[2];
rz(1.2606603) q[3];
sx q[3];
rz(-1.2536067) q[3];
sx q[3];
rz(-0.96084259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0925838) q[2];
sx q[2];
rz(-0.2139341) q[2];
sx q[2];
rz(1.5157549) q[2];
rz(1.9745291) q[3];
sx q[3];
rz(-1.5853106) q[3];
sx q[3];
rz(-1.0313755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9243069) q[0];
sx q[0];
rz(-1.8368245) q[0];
sx q[0];
rz(0.40689847) q[0];
rz(2.6869607) q[1];
sx q[1];
rz(-1.1063207) q[1];
sx q[1];
rz(2.8938821) q[1];
rz(1.0977279) q[2];
sx q[2];
rz(-1.4164783) q[2];
sx q[2];
rz(2.867792) q[2];
rz(1.3439988) q[3];
sx q[3];
rz(-0.26567017) q[3];
sx q[3];
rz(-2.9220823) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];