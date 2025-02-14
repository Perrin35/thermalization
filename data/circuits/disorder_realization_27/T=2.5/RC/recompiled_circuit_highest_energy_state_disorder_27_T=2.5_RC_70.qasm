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
rz(-1.6838411) q[0];
sx q[0];
rz(-0.10310752) q[0];
sx q[0];
rz(-2.4007894) q[0];
rz(2.9840772) q[1];
sx q[1];
rz(-0.73556757) q[1];
sx q[1];
rz(0.29677376) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5723853) q[0];
sx q[0];
rz(-0.57842904) q[0];
sx q[0];
rz(-2.3142772) q[0];
x q[1];
rz(0.27539092) q[2];
sx q[2];
rz(-0.73705929) q[2];
sx q[2];
rz(-0.136497) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.43208968) q[1];
sx q[1];
rz(-2.5910834) q[1];
sx q[1];
rz(0.014660346) q[1];
rz(1.2052676) q[3];
sx q[3];
rz(-2.6726279) q[3];
sx q[3];
rz(-2.7597357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.890581) q[2];
sx q[2];
rz(-3.0457532) q[2];
sx q[2];
rz(1.4061692) q[2];
rz(-0.74875325) q[3];
sx q[3];
rz(-0.65776062) q[3];
sx q[3];
rz(-2.9296618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9894079) q[0];
sx q[0];
rz(-2.2693372) q[0];
sx q[0];
rz(2.8023791) q[0];
rz(-0.57505125) q[1];
sx q[1];
rz(-1.9564068) q[1];
sx q[1];
rz(0.38160479) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1763827) q[0];
sx q[0];
rz(-1.3383703) q[0];
sx q[0];
rz(-1.9984972) q[0];
x q[1];
rz(0.26338844) q[2];
sx q[2];
rz(-1.8991422) q[2];
sx q[2];
rz(-1.8957317) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4424628) q[1];
sx q[1];
rz(-2.7425296) q[1];
sx q[1];
rz(1.5956053) q[1];
rz(-pi) q[2];
rz(-1.0142542) q[3];
sx q[3];
rz(-0.84662163) q[3];
sx q[3];
rz(-0.70412815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7316458) q[2];
sx q[2];
rz(-2.4531328) q[2];
sx q[2];
rz(0.48639578) q[2];
rz(2.5977123) q[3];
sx q[3];
rz(-2.0982274) q[3];
sx q[3];
rz(2.9774408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51089066) q[0];
sx q[0];
rz(-0.4158026) q[0];
sx q[0];
rz(-0.99935943) q[0];
rz(-2.4103145) q[1];
sx q[1];
rz(-0.46842289) q[1];
sx q[1];
rz(2.9670002) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.098487735) q[0];
sx q[0];
rz(-0.44666651) q[0];
sx q[0];
rz(1.5652324) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9810437) q[2];
sx q[2];
rz(-1.3957275) q[2];
sx q[2];
rz(-1.0700981) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.120003) q[1];
sx q[1];
rz(-1.5183416) q[1];
sx q[1];
rz(-1.4823556) q[1];
rz(-pi) q[2];
rz(-2.0295983) q[3];
sx q[3];
rz(-1.7672886) q[3];
sx q[3];
rz(-0.3199164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.38873765) q[2];
sx q[2];
rz(-2.2896705) q[2];
sx q[2];
rz(-1.277415) q[2];
rz(-0.82469213) q[3];
sx q[3];
rz(-0.84452355) q[3];
sx q[3];
rz(-2.9741014) q[3];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81553066) q[0];
sx q[0];
rz(-0.46634665) q[0];
sx q[0];
rz(-1.9933568) q[0];
rz(-2.8094021) q[1];
sx q[1];
rz(-2.5416608) q[1];
sx q[1];
rz(-2.0567599) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6116219) q[0];
sx q[0];
rz(-2.8154439) q[0];
sx q[0];
rz(1.1205313) q[0];
rz(-pi) q[1];
rz(0.81697322) q[2];
sx q[2];
rz(-1.8228056) q[2];
sx q[2];
rz(-1.1622045) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0177637) q[1];
sx q[1];
rz(-2.4014086) q[1];
sx q[1];
rz(-0.077880903) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74777909) q[3];
sx q[3];
rz(-0.59855312) q[3];
sx q[3];
rz(-2.7876496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.79160249) q[2];
sx q[2];
rz(-2.377066) q[2];
sx q[2];
rz(-0.0066268607) q[2];
rz(2.918112) q[3];
sx q[3];
rz(-1.0379182) q[3];
sx q[3];
rz(-2.3124783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4357736) q[0];
sx q[0];
rz(-2.3888102) q[0];
sx q[0];
rz(1.1821795) q[0];
rz(-1.301282) q[1];
sx q[1];
rz(-2.2186406) q[1];
sx q[1];
rz(-2.9822947) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1263016) q[0];
sx q[0];
rz(-1.8582004) q[0];
sx q[0];
rz(-1.5319583) q[0];
x q[1];
rz(-1.8865239) q[2];
sx q[2];
rz(-1.3906533) q[2];
sx q[2];
rz(0.9100998) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9827798) q[1];
sx q[1];
rz(-2.4430664) q[1];
sx q[1];
rz(-2.34312) q[1];
rz(-pi) q[2];
rz(2.0924657) q[3];
sx q[3];
rz(-2.4043466) q[3];
sx q[3];
rz(2.3741219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2948239) q[2];
sx q[2];
rz(-0.20192768) q[2];
sx q[2];
rz(-1.3429886) q[2];
rz(2.790847) q[3];
sx q[3];
rz(-1.1387768) q[3];
sx q[3];
rz(0.32575592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(0.89966929) q[0];
sx q[0];
rz(-2.3133008) q[0];
sx q[0];
rz(2.9747466) q[0];
rz(0.22653656) q[1];
sx q[1];
rz(-1.3275361) q[1];
sx q[1];
rz(2.66364) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95889752) q[0];
sx q[0];
rz(-1.0144466) q[0];
sx q[0];
rz(2.7378847) q[0];
rz(-pi) q[1];
rz(-1.9327546) q[2];
sx q[2];
rz(-0.88506341) q[2];
sx q[2];
rz(2.509523) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9354269) q[1];
sx q[1];
rz(-0.7784673) q[1];
sx q[1];
rz(2.4099518) q[1];
rz(-2.7109954) q[3];
sx q[3];
rz(-1.1599419) q[3];
sx q[3];
rz(1.2263067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6106674) q[2];
sx q[2];
rz(-1.3767367) q[2];
sx q[2];
rz(-1.8492071) q[2];
rz(-2.2409706) q[3];
sx q[3];
rz(-2.3714122) q[3];
sx q[3];
rz(1.5463411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6677299) q[0];
sx q[0];
rz(-0.97309363) q[0];
sx q[0];
rz(1.5245755) q[0];
rz(1.698311) q[1];
sx q[1];
rz(-2.2459005) q[1];
sx q[1];
rz(2.6406094) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.058280073) q[0];
sx q[0];
rz(-2.4078363) q[0];
sx q[0];
rz(-2.1408129) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47246859) q[2];
sx q[2];
rz(-0.61206619) q[2];
sx q[2];
rz(-1.2304359) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0308548) q[1];
sx q[1];
rz(-0.99910611) q[1];
sx q[1];
rz(-2.5377889) q[1];
rz(-0.37338169) q[3];
sx q[3];
rz(-2.8115559) q[3];
sx q[3];
rz(-0.16597834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.78918004) q[2];
sx q[2];
rz(-0.12426201) q[2];
sx q[2];
rz(2.6782356) q[2];
rz(-3.0739259) q[3];
sx q[3];
rz(-1.3031518) q[3];
sx q[3];
rz(2.9786003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8891334) q[0];
sx q[0];
rz(-1.8661789) q[0];
sx q[0];
rz(0.5109936) q[0];
rz(-0.64396089) q[1];
sx q[1];
rz(-2.0168346) q[1];
sx q[1];
rz(0.28800979) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0181684) q[0];
sx q[0];
rz(-2.9459758) q[0];
sx q[0];
rz(0.19851144) q[0];
rz(-pi) q[1];
rz(0.87784572) q[2];
sx q[2];
rz(-1.0509509) q[2];
sx q[2];
rz(-0.23990384) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0499514) q[1];
sx q[1];
rz(-2.5041951) q[1];
sx q[1];
rz(2.6693488) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7631825) q[3];
sx q[3];
rz(-1.6749951) q[3];
sx q[3];
rz(-0.7507594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1945232) q[2];
sx q[2];
rz(-1.1627407) q[2];
sx q[2];
rz(0.21491773) q[2];
rz(-1.8042709) q[3];
sx q[3];
rz(-2.603172) q[3];
sx q[3];
rz(-2.9609093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27338481) q[0];
sx q[0];
rz(-0.93623638) q[0];
sx q[0];
rz(-1.0816164) q[0];
rz(2.1514905) q[1];
sx q[1];
rz(-1.5847881) q[1];
sx q[1];
rz(0.33499151) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7107527) q[0];
sx q[0];
rz(-0.3246322) q[0];
sx q[0];
rz(1.6765094) q[0];
x q[1];
rz(1.5849466) q[2];
sx q[2];
rz(-0.42193896) q[2];
sx q[2];
rz(2.5463605) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44030467) q[1];
sx q[1];
rz(-3.0106643) q[1];
sx q[1];
rz(0.71016117) q[1];
rz(-pi) q[2];
rz(1.5906628) q[3];
sx q[3];
rz(-2.7405313) q[3];
sx q[3];
rz(0.10894081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.058502402) q[2];
sx q[2];
rz(-0.52512705) q[2];
sx q[2];
rz(-2.7527909) q[2];
rz(2.7723516) q[3];
sx q[3];
rz(-2.8856314) q[3];
sx q[3];
rz(2.2970439) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19287547) q[0];
sx q[0];
rz(-2.188864) q[0];
sx q[0];
rz(-0.70190758) q[0];
rz(-1.5319872) q[1];
sx q[1];
rz(-0.67370266) q[1];
sx q[1];
rz(-2.7598377) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88696357) q[0];
sx q[0];
rz(-1.6338631) q[0];
sx q[0];
rz(-3.0596759) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6959305) q[2];
sx q[2];
rz(-1.1627253) q[2];
sx q[2];
rz(3.0681075) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1211135) q[1];
sx q[1];
rz(-1.8217297) q[1];
sx q[1];
rz(-2.8931055) q[1];
rz(-pi) q[2];
rz(0.66017229) q[3];
sx q[3];
rz(-0.91554612) q[3];
sx q[3];
rz(1.2869664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6751042) q[2];
sx q[2];
rz(-2.1110822) q[2];
sx q[2];
rz(2.4934736) q[2];
rz(1.0068007) q[3];
sx q[3];
rz(-1.1454134) q[3];
sx q[3];
rz(-3.0692611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61160144) q[0];
sx q[0];
rz(-1.4016822) q[0];
sx q[0];
rz(0.66521426) q[0];
rz(-2.8499659) q[1];
sx q[1];
rz(-1.3529774) q[1];
sx q[1];
rz(-1.1381961) q[1];
rz(0.57403471) q[2];
sx q[2];
rz(-2.5435401) q[2];
sx q[2];
rz(0.5141573) q[2];
rz(1.9907822) q[3];
sx q[3];
rz(-1.715015) q[3];
sx q[3];
rz(0.7072995) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
