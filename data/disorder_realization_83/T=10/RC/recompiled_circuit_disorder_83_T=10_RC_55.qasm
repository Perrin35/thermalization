OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5873592) q[0];
sx q[0];
rz(-0.98266196) q[0];
sx q[0];
rz(-0.76173705) q[0];
rz(0.76454437) q[1];
sx q[1];
rz(-2.0643056) q[1];
sx q[1];
rz(0.74365562) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24510358) q[0];
sx q[0];
rz(-1.5050383) q[0];
sx q[0];
rz(-1.3301646) q[0];
rz(-pi) q[1];
rz(0.6526297) q[2];
sx q[2];
rz(-1.666781) q[2];
sx q[2];
rz(1.8362311) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.88120645) q[1];
sx q[1];
rz(-1.2938061) q[1];
sx q[1];
rz(1.1587515) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7573962) q[3];
sx q[3];
rz(-0.36986923) q[3];
sx q[3];
rz(2.6922525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4347697) q[2];
sx q[2];
rz(-2.7396024) q[2];
sx q[2];
rz(-3.0337231) q[2];
rz(-0.14262959) q[3];
sx q[3];
rz(-1.7248036) q[3];
sx q[3];
rz(2.4690348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0649439) q[0];
sx q[0];
rz(-2.3363484) q[0];
sx q[0];
rz(0.23072492) q[0];
rz(-1.8143066) q[1];
sx q[1];
rz(-2.4704411) q[1];
sx q[1];
rz(0.040963106) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3687392) q[0];
sx q[0];
rz(-0.29631796) q[0];
sx q[0];
rz(-0.98451891) q[0];
rz(-2.3054753) q[2];
sx q[2];
rz(-0.40514075) q[2];
sx q[2];
rz(-0.9678313) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7663824) q[1];
sx q[1];
rz(-0.76645215) q[1];
sx q[1];
rz(2.4541928) q[1];
rz(-pi) q[2];
rz(-0.81539865) q[3];
sx q[3];
rz(-2.6147463) q[3];
sx q[3];
rz(-2.6032053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9791947) q[2];
sx q[2];
rz(-1.6684063) q[2];
sx q[2];
rz(-2.5734148) q[2];
rz(0.5125106) q[3];
sx q[3];
rz(-2.5458702) q[3];
sx q[3];
rz(0.56186831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42999643) q[0];
sx q[0];
rz(-1.5189518) q[0];
sx q[0];
rz(0.45021737) q[0];
rz(-1.2954767) q[1];
sx q[1];
rz(-1.1643012) q[1];
sx q[1];
rz(0.67726642) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5325461) q[0];
sx q[0];
rz(-2.9653774) q[0];
sx q[0];
rz(2.8410068) q[0];
rz(-pi) q[1];
rz(0.1707503) q[2];
sx q[2];
rz(-3.1351334) q[2];
sx q[2];
rz(3.1250172) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3213768) q[1];
sx q[1];
rz(-1.4396588) q[1];
sx q[1];
rz(1.7507491) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47795313) q[3];
sx q[3];
rz(-1.2198997) q[3];
sx q[3];
rz(-1.4357476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9230817) q[2];
sx q[2];
rz(-1.9873025) q[2];
sx q[2];
rz(-3.0947321) q[2];
rz(0.81165195) q[3];
sx q[3];
rz(-2.4620158) q[3];
sx q[3];
rz(2.9714382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99455225) q[0];
sx q[0];
rz(-3.0299598) q[0];
sx q[0];
rz(3.0901093) q[0];
rz(2.6507846) q[1];
sx q[1];
rz(-2.1665116) q[1];
sx q[1];
rz(-1.1725918) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8797982) q[0];
sx q[0];
rz(-0.9167295) q[0];
sx q[0];
rz(1.4021224) q[0];
rz(-pi) q[1];
rz(-1.5863717) q[2];
sx q[2];
rz(-2.1767463) q[2];
sx q[2];
rz(1.7473999) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.51845156) q[1];
sx q[1];
rz(-1.8703096) q[1];
sx q[1];
rz(1.9407942) q[1];
rz(-pi) q[2];
rz(2.0026607) q[3];
sx q[3];
rz(-1.6383088) q[3];
sx q[3];
rz(-2.1907012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.88901687) q[2];
sx q[2];
rz(-2.2106407) q[2];
sx q[2];
rz(2.6546997) q[2];
rz(1.1934818) q[3];
sx q[3];
rz(-2.8119757) q[3];
sx q[3];
rz(-0.025432767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30381969) q[0];
sx q[0];
rz(-1.0205512) q[0];
sx q[0];
rz(1.8160965) q[0];
rz(0.59108132) q[1];
sx q[1];
rz(-2.4609844) q[1];
sx q[1];
rz(0.65471929) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7514873) q[0];
sx q[0];
rz(-2.562398) q[0];
sx q[0];
rz(0.67436995) q[0];
rz(-pi) q[1];
rz(1.3833369) q[2];
sx q[2];
rz(-1.9831295) q[2];
sx q[2];
rz(-2.2802441) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3760738) q[1];
sx q[1];
rz(-2.3507833) q[1];
sx q[1];
rz(2.5942624) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4773024) q[3];
sx q[3];
rz(-1.1870541) q[3];
sx q[3];
rz(1.547471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6902265) q[2];
sx q[2];
rz(-1.9407242) q[2];
sx q[2];
rz(-2.5704685) q[2];
rz(2.5518104) q[3];
sx q[3];
rz(-2.6815806) q[3];
sx q[3];
rz(0.067972876) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1968483) q[0];
sx q[0];
rz(-1.1905043) q[0];
sx q[0];
rz(-2.8425472) q[0];
rz(-1.3202745) q[1];
sx q[1];
rz(-0.25779217) q[1];
sx q[1];
rz(1.6437795) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70443557) q[0];
sx q[0];
rz(-2.2309982) q[0];
sx q[0];
rz(-3.0755088) q[0];
rz(-pi) q[1];
rz(-2.5619179) q[2];
sx q[2];
rz(-2.0849166) q[2];
sx q[2];
rz(-0.74619734) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.99243473) q[1];
sx q[1];
rz(-0.5545534) q[1];
sx q[1];
rz(-2.1402332) q[1];
x q[2];
rz(2.8619814) q[3];
sx q[3];
rz(-1.8108484) q[3];
sx q[3];
rz(0.9595426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.51612878) q[2];
sx q[2];
rz(-1.7136145) q[2];
sx q[2];
rz(-0.027475474) q[2];
rz(2.6190858) q[3];
sx q[3];
rz(-2.3404739) q[3];
sx q[3];
rz(2.3251422) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7291173) q[0];
sx q[0];
rz(-1.1183879) q[0];
sx q[0];
rz(-3.0274042) q[0];
rz(2.1633637) q[1];
sx q[1];
rz(-2.5949635) q[1];
sx q[1];
rz(-0.79089975) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40490155) q[0];
sx q[0];
rz(-1.6653403) q[0];
sx q[0];
rz(2.1928284) q[0];
x q[1];
rz(2.3979264) q[2];
sx q[2];
rz(-1.3197834) q[2];
sx q[2];
rz(0.92168346) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.447688) q[1];
sx q[1];
rz(-1.4920007) q[1];
sx q[1];
rz(-1.2055956) q[1];
rz(-2.0049719) q[3];
sx q[3];
rz(-2.4063111) q[3];
sx q[3];
rz(0.47664627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.87166446) q[2];
sx q[2];
rz(-2.0285138) q[2];
sx q[2];
rz(-2.2765735) q[2];
rz(-0.67251742) q[3];
sx q[3];
rz(-2.0130242) q[3];
sx q[3];
rz(-0.095656693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.4928116) q[0];
sx q[0];
rz(-0.5195986) q[0];
sx q[0];
rz(-2.6742324) q[0];
rz(-2.6043747) q[1];
sx q[1];
rz(-0.98058128) q[1];
sx q[1];
rz(0.25407243) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2652889) q[0];
sx q[0];
rz(-0.2067925) q[0];
sx q[0];
rz(-0.32949038) q[0];
x q[1];
rz(-3.1363856) q[2];
sx q[2];
rz(-1.5723096) q[2];
sx q[2];
rz(3.0031799) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.961812) q[1];
sx q[1];
rz(-1.6976377) q[1];
sx q[1];
rz(2.8282053) q[1];
rz(-pi) q[2];
rz(-1.6510321) q[3];
sx q[3];
rz(-2.1685765) q[3];
sx q[3];
rz(2.762592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.49999985) q[2];
sx q[2];
rz(-0.77074146) q[2];
sx q[2];
rz(0.33561486) q[2];
rz(2.8149758) q[3];
sx q[3];
rz(-0.89544046) q[3];
sx q[3];
rz(1.4183104) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083387233) q[0];
sx q[0];
rz(-0.21033062) q[0];
sx q[0];
rz(-1.0539508) q[0];
rz(-0.15696934) q[1];
sx q[1];
rz(-1.418768) q[1];
sx q[1];
rz(-2.1597247) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0197524) q[0];
sx q[0];
rz(-1.6822018) q[0];
sx q[0];
rz(2.7287672) q[0];
x q[1];
rz(-3.0332546) q[2];
sx q[2];
rz(-1.0045369) q[2];
sx q[2];
rz(-0.94891753) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2541483) q[1];
sx q[1];
rz(-1.5132628) q[1];
sx q[1];
rz(-1.385034) q[1];
x q[2];
rz(-1.8198265) q[3];
sx q[3];
rz(-0.8753652) q[3];
sx q[3];
rz(1.4679366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92423576) q[2];
sx q[2];
rz(-2.0841667) q[2];
sx q[2];
rz(-2.1441933) q[2];
rz(-2.635397) q[3];
sx q[3];
rz(-2.180407) q[3];
sx q[3];
rz(0.034328073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2845594) q[0];
sx q[0];
rz(-0.71145809) q[0];
sx q[0];
rz(-0.51666623) q[0];
rz(2.754028) q[1];
sx q[1];
rz(-2.0811847) q[1];
sx q[1];
rz(0.26836747) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3049406) q[0];
sx q[0];
rz(-0.9965082) q[0];
sx q[0];
rz(-0.71026295) q[0];
rz(2.1568314) q[2];
sx q[2];
rz(-1.9433937) q[2];
sx q[2];
rz(-0.44024703) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.25071535) q[1];
sx q[1];
rz(-2.4836575) q[1];
sx q[1];
rz(0.4472181) q[1];
x q[2];
rz(-2.8769365) q[3];
sx q[3];
rz(-1.0088682) q[3];
sx q[3];
rz(0.86268007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2051852) q[2];
sx q[2];
rz(-0.78339094) q[2];
sx q[2];
rz(0.56383413) q[2];
rz(-1.9514203) q[3];
sx q[3];
rz(-2.1824013) q[3];
sx q[3];
rz(2.4998375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.664809) q[0];
sx q[0];
rz(-1.2510779) q[0];
sx q[0];
rz(-1.0673987) q[0];
rz(1.339636) q[1];
sx q[1];
rz(-1.4383153) q[1];
sx q[1];
rz(-1.7972606) q[1];
rz(-1.2673169) q[2];
sx q[2];
rz(-1.0721285) q[2];
sx q[2];
rz(-0.68821651) q[2];
rz(-1.7715122) q[3];
sx q[3];
rz(-1.8772535) q[3];
sx q[3];
rz(1.5542961) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
