OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2296978) q[0];
sx q[0];
rz(-1.8954281) q[0];
sx q[0];
rz(1.5204313) q[0];
rz(0.35960943) q[1];
sx q[1];
rz(-2.8728027) q[1];
sx q[1];
rz(-0.51528817) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7069076) q[0];
sx q[0];
rz(-1.3680172) q[0];
sx q[0];
rz(-2.7903656) q[0];
x q[1];
rz(-1.2688387) q[2];
sx q[2];
rz(-1.0450604) q[2];
sx q[2];
rz(2.5765004) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88564116) q[1];
sx q[1];
rz(-0.4547387) q[1];
sx q[1];
rz(-0.498147) q[1];
rz(0.73909594) q[3];
sx q[3];
rz(-2.1198273) q[3];
sx q[3];
rz(-2.9573553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9109853) q[2];
sx q[2];
rz(-1.5585941) q[2];
sx q[2];
rz(0.67414635) q[2];
rz(-2.9437183) q[3];
sx q[3];
rz(-1.2400235) q[3];
sx q[3];
rz(1.5052634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3182217) q[0];
sx q[0];
rz(-0.30105337) q[0];
sx q[0];
rz(-2.9252885) q[0];
rz(-0.11953106) q[1];
sx q[1];
rz(-1.1263589) q[1];
sx q[1];
rz(-1.4215887) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57935539) q[0];
sx q[0];
rz(-0.52412141) q[0];
sx q[0];
rz(1.7299132) q[0];
x q[1];
rz(1.7825215) q[2];
sx q[2];
rz(-1.1781409) q[2];
sx q[2];
rz(-0.77048877) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.51829631) q[1];
sx q[1];
rz(-2.5824248) q[1];
sx q[1];
rz(-2.1914047) q[1];
rz(-pi) q[2];
rz(-2.414538) q[3];
sx q[3];
rz(-2.0306572) q[3];
sx q[3];
rz(2.1279631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8372832) q[2];
sx q[2];
rz(-1.195636) q[2];
sx q[2];
rz(2.1514814) q[2];
rz(-0.75634161) q[3];
sx q[3];
rz(-0.6476616) q[3];
sx q[3];
rz(2.9653449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15683098) q[0];
sx q[0];
rz(-2.3909843) q[0];
sx q[0];
rz(-2.005715) q[0];
rz(-1.0095949) q[1];
sx q[1];
rz(-1.986809) q[1];
sx q[1];
rz(0.44581595) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71190055) q[0];
sx q[0];
rz(-1.4662379) q[0];
sx q[0];
rz(0.069653102) q[0];
rz(1.5405212) q[2];
sx q[2];
rz(-0.86756268) q[2];
sx q[2];
rz(1.3959259) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.34796068) q[1];
sx q[1];
rz(-2.1588481) q[1];
sx q[1];
rz(-2.6403202) q[1];
rz(-1.5169083) q[3];
sx q[3];
rz(-1.1240376) q[3];
sx q[3];
rz(-1.7487546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3710215) q[2];
sx q[2];
rz(-1.0022481) q[2];
sx q[2];
rz(0.64024964) q[2];
rz(-0.69784969) q[3];
sx q[3];
rz(-1.4638487) q[3];
sx q[3];
rz(-0.33997926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3014389) q[0];
sx q[0];
rz(-0.88742632) q[0];
sx q[0];
rz(-1.812717) q[0];
rz(-1.9245194) q[1];
sx q[1];
rz(-2.4220059) q[1];
sx q[1];
rz(0.77635366) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8912131) q[0];
sx q[0];
rz(-2.4107412) q[0];
sx q[0];
rz(-2.5688085) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59036915) q[2];
sx q[2];
rz(-1.9861097) q[2];
sx q[2];
rz(1.6990952) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1963819) q[1];
sx q[1];
rz(-0.27881611) q[1];
sx q[1];
rz(1.6120738) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68459529) q[3];
sx q[3];
rz(-1.5696313) q[3];
sx q[3];
rz(-2.5161285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2886469) q[2];
sx q[2];
rz(-2.8198346) q[2];
sx q[2];
rz(-1.9256437) q[2];
rz(-1.1207885) q[3];
sx q[3];
rz(-1.8782764) q[3];
sx q[3];
rz(-2.2585675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23585606) q[0];
sx q[0];
rz(-2.164542) q[0];
sx q[0];
rz(2.6781154) q[0];
rz(2.0526759) q[1];
sx q[1];
rz(-1.2605647) q[1];
sx q[1];
rz(-1.3190528) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80281335) q[0];
sx q[0];
rz(-1.3283214) q[0];
sx q[0];
rz(0.26870299) q[0];
rz(-pi) q[1];
rz(1.7402788) q[2];
sx q[2];
rz(-0.78022829) q[2];
sx q[2];
rz(1.4223157) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.12015776) q[1];
sx q[1];
rz(-1.48367) q[1];
sx q[1];
rz(2.1539262) q[1];
rz(-pi) q[2];
rz(1.1518258) q[3];
sx q[3];
rz(-1.4052466) q[3];
sx q[3];
rz(-2.3922684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9702381) q[2];
sx q[2];
rz(-0.78533185) q[2];
sx q[2];
rz(0.63924092) q[2];
rz(0.81131896) q[3];
sx q[3];
rz(-2.6526484) q[3];
sx q[3];
rz(-1.5343687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92523471) q[0];
sx q[0];
rz(-0.41330591) q[0];
sx q[0];
rz(-2.9226724) q[0];
rz(-1.9256598) q[1];
sx q[1];
rz(-0.464012) q[1];
sx q[1];
rz(-1.6485515) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3603148) q[0];
sx q[0];
rz(-2.6584315) q[0];
sx q[0];
rz(-0.29668087) q[0];
x q[1];
rz(-2.7240407) q[2];
sx q[2];
rz(-1.5680015) q[2];
sx q[2];
rz(-0.33139834) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.56144372) q[1];
sx q[1];
rz(-0.83858788) q[1];
sx q[1];
rz(1.5484018) q[1];
x q[2];
rz(-1.5739176) q[3];
sx q[3];
rz(-1.471648) q[3];
sx q[3];
rz(-1.329545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1047989) q[2];
sx q[2];
rz(-1.3497738) q[2];
sx q[2];
rz(1.8297423) q[2];
rz(-1.5051684) q[3];
sx q[3];
rz(-1.2994095) q[3];
sx q[3];
rz(-2.6817491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5055607) q[0];
sx q[0];
rz(-0.4902896) q[0];
sx q[0];
rz(-1.3469411) q[0];
rz(-1.3268283) q[1];
sx q[1];
rz(-2.3074) q[1];
sx q[1];
rz(-0.38536513) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13305328) q[0];
sx q[0];
rz(-0.36046777) q[0];
sx q[0];
rz(0.089715606) q[0];
x q[1];
rz(-0.48882882) q[2];
sx q[2];
rz(-1.0932473) q[2];
sx q[2];
rz(1.0350641) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7984385) q[1];
sx q[1];
rz(-2.0166409) q[1];
sx q[1];
rz(-1.6002161) q[1];
rz(-pi) q[2];
rz(-2.3919098) q[3];
sx q[3];
rz(-2.5307641) q[3];
sx q[3];
rz(1.2884097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4555326) q[2];
sx q[2];
rz(-2.3251688) q[2];
sx q[2];
rz(-1.653999) q[2];
rz(0.16573302) q[3];
sx q[3];
rz(-2.3159852) q[3];
sx q[3];
rz(2.9674271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1626749) q[0];
sx q[0];
rz(-2.5161777) q[0];
sx q[0];
rz(-0.22853525) q[0];
rz(-2.8288815) q[1];
sx q[1];
rz(-0.87930185) q[1];
sx q[1];
rz(1.762134) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40423597) q[0];
sx q[0];
rz(-1.9725058) q[0];
sx q[0];
rz(-1.4866923) q[0];
rz(-1.3993456) q[2];
sx q[2];
rz(-1.2510692) q[2];
sx q[2];
rz(0.4188183) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5818565) q[1];
sx q[1];
rz(-0.35071555) q[1];
sx q[1];
rz(-1.4697671) q[1];
rz(2.5614235) q[3];
sx q[3];
rz(-2.3456367) q[3];
sx q[3];
rz(0.43886504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.10869965) q[2];
sx q[2];
rz(-1.0445107) q[2];
sx q[2];
rz(-1.0408939) q[2];
rz(-2.6563472) q[3];
sx q[3];
rz(-1.3121366) q[3];
sx q[3];
rz(0.41788873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27337209) q[0];
sx q[0];
rz(-2.1594248) q[0];
sx q[0];
rz(-2.3305273) q[0];
rz(-1.2366933) q[1];
sx q[1];
rz(-0.86634723) q[1];
sx q[1];
rz(2.9467357) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4351866) q[0];
sx q[0];
rz(-1.331745) q[0];
sx q[0];
rz(-1.8349706) q[0];
rz(0.92020184) q[2];
sx q[2];
rz(-2.0363931) q[2];
sx q[2];
rz(-2.1289785) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5734308) q[1];
sx q[1];
rz(-0.57933148) q[1];
sx q[1];
rz(2.2327077) q[1];
rz(2.2556188) q[3];
sx q[3];
rz(-0.77903226) q[3];
sx q[3];
rz(-2.1942735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1607232) q[2];
sx q[2];
rz(-1.2474493) q[2];
sx q[2];
rz(-0.54171872) q[2];
rz(-1.3228275) q[3];
sx q[3];
rz(-1.8018689) q[3];
sx q[3];
rz(2.8790348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5498891) q[0];
sx q[0];
rz(-1.6240969) q[0];
sx q[0];
rz(0.24895915) q[0];
rz(1.2173563) q[1];
sx q[1];
rz(-1.267642) q[1];
sx q[1];
rz(-2.731146) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9684943) q[0];
sx q[0];
rz(-1.9763034) q[0];
sx q[0];
rz(-1.9963229) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4798231) q[2];
sx q[2];
rz(-1.5814648) q[2];
sx q[2];
rz(-1.7144698) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.96406781) q[1];
sx q[1];
rz(-0.52492889) q[1];
sx q[1];
rz(-1.0704416) q[1];
rz(1.0656359) q[3];
sx q[3];
rz(-2.1200051) q[3];
sx q[3];
rz(2.4545936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2534788) q[2];
sx q[2];
rz(-1.9733182) q[2];
sx q[2];
rz(-1.1178364) q[2];
rz(0.11876336) q[3];
sx q[3];
rz(-0.6898841) q[3];
sx q[3];
rz(1.9411055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4111907) q[0];
sx q[0];
rz(-1.406519) q[0];
sx q[0];
rz(2.7997959) q[0];
rz(0.047601184) q[1];
sx q[1];
rz(-1.2536512) q[1];
sx q[1];
rz(1.8258078) q[1];
rz(-0.81472266) q[2];
sx q[2];
rz(-1.2263032) q[2];
sx q[2];
rz(-2.6487614) q[2];
rz(-2.4763465) q[3];
sx q[3];
rz(-0.86626296) q[3];
sx q[3];
rz(-2.5588425) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
