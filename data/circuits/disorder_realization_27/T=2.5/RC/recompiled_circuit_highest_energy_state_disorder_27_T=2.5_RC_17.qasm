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
rz(0.74080324) q[0];
rz(-0.15751547) q[1];
sx q[1];
rz(3.8771602) q[1];
sx q[1];
rz(9.1280042) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5723853) q[0];
sx q[0];
rz(-0.57842904) q[0];
sx q[0];
rz(-0.82731547) q[0];
x q[1];
rz(-2.4236083) q[2];
sx q[2];
rz(-1.7545926) q[2];
sx q[2];
rz(-1.2280457) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.41488838) q[1];
sx q[1];
rz(-1.0203531) q[1];
sx q[1];
rz(-1.5617983) q[1];
rz(-pi) q[2];
rz(-1.2052676) q[3];
sx q[3];
rz(-0.46896471) q[3];
sx q[3];
rz(0.381857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.2510117) q[2];
sx q[2];
rz(-3.0457532) q[2];
sx q[2];
rz(-1.4061692) q[2];
rz(-0.74875325) q[3];
sx q[3];
rz(-0.65776062) q[3];
sx q[3];
rz(-2.9296618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9894079) q[0];
sx q[0];
rz(-2.2693372) q[0];
sx q[0];
rz(-0.33921355) q[0];
rz(0.57505125) q[1];
sx q[1];
rz(-1.1851858) q[1];
sx q[1];
rz(0.38160479) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71020102) q[0];
sx q[0];
rz(-1.1553191) q[0];
sx q[0];
rz(0.25449591) q[0];
rz(2.223338) q[2];
sx q[2];
rz(-2.7236669) q[2];
sx q[2];
rz(-1.9422836) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4424628) q[1];
sx q[1];
rz(-2.7425296) q[1];
sx q[1];
rz(-1.5459874) q[1];
rz(-pi) q[2];
rz(0.53839243) q[3];
sx q[3];
rz(-0.88148553) q[3];
sx q[3];
rz(-1.6834843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7316458) q[2];
sx q[2];
rz(-0.68845981) q[2];
sx q[2];
rz(-0.48639578) q[2];
rz(-2.5977123) q[3];
sx q[3];
rz(-1.0433652) q[3];
sx q[3];
rz(-0.16415183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.630702) q[0];
sx q[0];
rz(-2.7257901) q[0];
sx q[0];
rz(2.1422332) q[0];
rz(-0.73127812) q[1];
sx q[1];
rz(-2.6731698) q[1];
sx q[1];
rz(2.9670002) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.664266) q[0];
sx q[0];
rz(-1.5731997) q[0];
sx q[0];
rz(-1.1241358) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9810437) q[2];
sx q[2];
rz(-1.3957275) q[2];
sx q[2];
rz(-2.0714945) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5445575) q[1];
sx q[1];
rz(-1.4824776) q[1];
sx q[1];
rz(-3.0889325) q[1];
rz(-1.1483795) q[3];
sx q[3];
rz(-0.49634051) q[3];
sx q[3];
rz(1.6272735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.38873765) q[2];
sx q[2];
rz(-0.85192215) q[2];
sx q[2];
rz(-1.277415) q[2];
rz(2.3169005) q[3];
sx q[3];
rz(-0.84452355) q[3];
sx q[3];
rz(-2.9741014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.326062) q[0];
sx q[0];
rz(-2.675246) q[0];
sx q[0];
rz(1.1482358) q[0];
rz(2.8094021) q[1];
sx q[1];
rz(-0.59993184) q[1];
sx q[1];
rz(-2.0567599) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0834616) q[0];
sx q[0];
rz(-1.2781743) q[0];
sx q[0];
rz(-0.14614848) q[0];
rz(1.930619) q[2];
sx q[2];
rz(-2.3546989) q[2];
sx q[2];
rz(-0.66821276) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22913223) q[1];
sx q[1];
rz(-0.83337754) q[1];
sx q[1];
rz(1.6417437) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3938136) q[3];
sx q[3];
rz(-2.5430395) q[3];
sx q[3];
rz(-0.35394305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.79160249) q[2];
sx q[2];
rz(-0.76452667) q[2];
sx q[2];
rz(-3.1349658) q[2];
rz(0.22348063) q[3];
sx q[3];
rz(-1.0379182) q[3];
sx q[3];
rz(2.3124783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4357736) q[0];
sx q[0];
rz(-0.75278246) q[0];
sx q[0];
rz(-1.1821795) q[0];
rz(-1.8403107) q[1];
sx q[1];
rz(-0.92295206) q[1];
sx q[1];
rz(-2.9822947) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2625339) q[0];
sx q[0];
rz(-2.8516483) q[0];
sx q[0];
rz(-3.0109809) q[0];
rz(-pi) q[1];
rz(-2.9522991) q[2];
sx q[2];
rz(-1.2603501) q[2];
sx q[2];
rz(-0.60223641) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9827798) q[1];
sx q[1];
rz(-0.69852622) q[1];
sx q[1];
rz(0.79847269) q[1];
x q[2];
rz(2.2377272) q[3];
sx q[3];
rz(-1.9124036) q[3];
sx q[3];
rz(-1.9359534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.84676877) q[2];
sx q[2];
rz(-0.20192768) q[2];
sx q[2];
rz(-1.798604) q[2];
rz(2.790847) q[3];
sx q[3];
rz(-1.1387768) q[3];
sx q[3];
rz(0.32575592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2419234) q[0];
sx q[0];
rz(-2.3133008) q[0];
sx q[0];
rz(-0.1668461) q[0];
rz(-0.22653656) q[1];
sx q[1];
rz(-1.8140565) q[1];
sx q[1];
rz(2.66364) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8628381) q[0];
sx q[0];
rz(-2.4669381) q[0];
sx q[0];
rz(1.007403) q[0];
x q[1];
rz(-2.4228335) q[2];
sx q[2];
rz(-1.2931817) q[2];
sx q[2];
rz(1.1740545) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8357827) q[1];
sx q[1];
rz(-2.1205582) q[1];
sx q[1];
rz(2.1533986) q[1];
x q[2];
rz(2.7109954) q[3];
sx q[3];
rz(-1.1599419) q[3];
sx q[3];
rz(-1.2263067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.53092521) q[2];
sx q[2];
rz(-1.764856) q[2];
sx q[2];
rz(1.8492071) q[2];
rz(2.2409706) q[3];
sx q[3];
rz(-2.3714122) q[3];
sx q[3];
rz(-1.5463411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47386277) q[0];
sx q[0];
rz(-0.97309363) q[0];
sx q[0];
rz(-1.6170172) q[0];
rz(1.4432817) q[1];
sx q[1];
rz(-2.2459005) q[1];
sx q[1];
rz(-2.6406094) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4878236) q[0];
sx q[0];
rz(-2.1697576) q[0];
sx q[0];
rz(2.6887211) q[0];
rz(-2.5829264) q[2];
sx q[2];
rz(-1.8353454) q[2];
sx q[2];
rz(-3.0857744) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9364719) q[1];
sx q[1];
rz(-2.3354482) q[1];
sx q[1];
rz(-2.2938726) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30876183) q[3];
sx q[3];
rz(-1.4523066) q[3];
sx q[3];
rz(-2.091696) q[3];
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
rz(-1.8384408) q[3];
sx q[3];
rz(0.1629924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25245923) q[0];
sx q[0];
rz(-1.8661789) q[0];
sx q[0];
rz(-0.5109936) q[0];
rz(-2.4976318) q[1];
sx q[1];
rz(-1.124758) q[1];
sx q[1];
rz(-2.8535829) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.889042) q[0];
sx q[0];
rz(-1.5324549) q[0];
sx q[0];
rz(-2.949723) q[0];
rz(-pi) q[1];
rz(0.87784572) q[2];
sx q[2];
rz(-1.0509509) q[2];
sx q[2];
rz(-0.23990384) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0102699) q[1];
sx q[1];
rz(-1.2966709) q[1];
sx q[1];
rz(-2.5585973) q[1];
rz(-pi) q[2];
rz(2.8657416) q[3];
sx q[3];
rz(-0.39182651) q[3];
sx q[3];
rz(-1.0759169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1945232) q[2];
sx q[2];
rz(-1.9788519) q[2];
sx q[2];
rz(-0.21491773) q[2];
rz(1.8042709) q[3];
sx q[3];
rz(-2.603172) q[3];
sx q[3];
rz(2.9609093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27338481) q[0];
sx q[0];
rz(-2.2053563) q[0];
sx q[0];
rz(-2.0599763) q[0];
rz(-2.1514905) q[1];
sx q[1];
rz(-1.5847881) q[1];
sx q[1];
rz(-0.33499151) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2401857) q[0];
sx q[0];
rz(-1.5371345) q[0];
sx q[0];
rz(-1.2478527) q[0];
rz(-1.556646) q[2];
sx q[2];
rz(-2.7196537) q[2];
sx q[2];
rz(0.59523216) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8674842) q[1];
sx q[1];
rz(-1.47164) q[1];
sx q[1];
rz(1.4851557) q[1];
rz(1.5509299) q[3];
sx q[3];
rz(-2.7405313) q[3];
sx q[3];
rz(3.0326518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0830903) q[2];
sx q[2];
rz(-2.6164656) q[2];
sx q[2];
rz(-0.38880175) q[2];
rz(0.36924103) q[3];
sx q[3];
rz(-2.8856314) q[3];
sx q[3];
rz(-2.2970439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19287547) q[0];
sx q[0];
rz(-0.95272869) q[0];
sx q[0];
rz(-2.4396851) q[0];
rz(1.6096055) q[1];
sx q[1];
rz(-2.46789) q[1];
sx q[1];
rz(-0.38175499) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88696357) q[0];
sx q[0];
rz(-1.5077295) q[0];
sx q[0];
rz(3.0596759) q[0];
rz(-pi) q[1];
rz(-2.8605531) q[2];
sx q[2];
rz(-0.42579415) q[2];
sx q[2];
rz(0.38044924) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4873969) q[1];
sx q[1];
rz(-1.3302478) q[1];
sx q[1];
rz(1.3122561) q[1];
rz(-0.89729805) q[3];
sx q[3];
rz(-2.2475911) q[3];
sx q[3];
rz(0.38161665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6751042) q[2];
sx q[2];
rz(-1.0305104) q[2];
sx q[2];
rz(-0.64811903) q[2];
rz(-2.134792) q[3];
sx q[3];
rz(-1.9961793) q[3];
sx q[3];
rz(-0.072331585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61160144) q[0];
sx q[0];
rz(-1.4016822) q[0];
sx q[0];
rz(0.66521426) q[0];
rz(2.8499659) q[1];
sx q[1];
rz(-1.7886152) q[1];
sx q[1];
rz(2.0033966) q[1];
rz(1.2164581) q[2];
sx q[2];
rz(-2.0632498) q[2];
sx q[2];
rz(-1.9634631) q[2];
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
