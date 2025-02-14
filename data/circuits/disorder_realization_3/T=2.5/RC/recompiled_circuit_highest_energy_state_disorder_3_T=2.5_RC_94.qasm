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
rz(-2.3961177) q[0];
sx q[0];
rz(-2.5479654) q[0];
sx q[0];
rz(2.797085) q[0];
rz(1.1747107) q[1];
sx q[1];
rz(-2.7732958) q[1];
sx q[1];
rz(-2.5370497) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9637361) q[0];
sx q[0];
rz(-2.3174598) q[0];
sx q[0];
rz(-0.080029052) q[0];
rz(-pi) q[1];
x q[1];
rz(3.072158) q[2];
sx q[2];
rz(-1.3926818) q[2];
sx q[2];
rz(1.3645862) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4010594) q[1];
sx q[1];
rz(-1.9895253) q[1];
sx q[1];
rz(1.5685521) q[1];
rz(-pi) q[2];
rz(2.3326268) q[3];
sx q[3];
rz(-2.3936317) q[3];
sx q[3];
rz(-2.7873519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.979226) q[2];
sx q[2];
rz(-1.484885) q[2];
sx q[2];
rz(1.7171198) q[2];
rz(-3.034397) q[3];
sx q[3];
rz(-0.19459477) q[3];
sx q[3];
rz(0.9596107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2694038) q[0];
sx q[0];
rz(-0.93282455) q[0];
sx q[0];
rz(1.3099571) q[0];
rz(0.45920363) q[1];
sx q[1];
rz(-1.2745067) q[1];
sx q[1];
rz(-2.8244663) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.223359) q[0];
sx q[0];
rz(-0.14592136) q[0];
sx q[0];
rz(-0.33182611) q[0];
rz(1.5731847) q[2];
sx q[2];
rz(-1.6206073) q[2];
sx q[2];
rz(-0.43149155) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9307946) q[1];
sx q[1];
rz(-2.5902777) q[1];
sx q[1];
rz(-3.0153794) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7789654) q[3];
sx q[3];
rz(-2.31593) q[3];
sx q[3];
rz(-2.7221808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.7121048) q[2];
sx q[2];
rz(-1.2965344) q[2];
sx q[2];
rz(-2.3340161) q[2];
rz(0.56337041) q[3];
sx q[3];
rz(-2.6475776) q[3];
sx q[3];
rz(-1.0529244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6323557) q[0];
sx q[0];
rz(-0.18993264) q[0];
sx q[0];
rz(0.83455363) q[0];
rz(-2.6368311) q[1];
sx q[1];
rz(-0.36634645) q[1];
sx q[1];
rz(-1.6827481) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8887688) q[0];
sx q[0];
rz(-1.3712168) q[0];
sx q[0];
rz(0.56364023) q[0];
x q[1];
rz(0.81291109) q[2];
sx q[2];
rz(-0.92070075) q[2];
sx q[2];
rz(1.7342099) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9533938) q[1];
sx q[1];
rz(-2.2911304) q[1];
sx q[1];
rz(-1.6745643) q[1];
rz(0.97455131) q[3];
sx q[3];
rz(-0.85874346) q[3];
sx q[3];
rz(-0.94030583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0316281) q[2];
sx q[2];
rz(-1.6692903) q[2];
sx q[2];
rz(-2.7785981) q[2];
rz(-1.1686769) q[3];
sx q[3];
rz(-0.76255885) q[3];
sx q[3];
rz(-0.76538435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7211683) q[0];
sx q[0];
rz(-0.52614251) q[0];
sx q[0];
rz(2.5508733) q[0];
rz(-2.4144454) q[1];
sx q[1];
rz(-2.7666481) q[1];
sx q[1];
rz(2.4730543) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1626133) q[0];
sx q[0];
rz(-0.091446459) q[0];
sx q[0];
rz(-1.0981111) q[0];
rz(-pi) q[1];
rz(-2.3756251) q[2];
sx q[2];
rz(-1.2184288) q[2];
sx q[2];
rz(-2.3841928) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.82610735) q[1];
sx q[1];
rz(-0.8902227) q[1];
sx q[1];
rz(2.1904551) q[1];
x q[2];
rz(-2.5741124) q[3];
sx q[3];
rz(-0.82616185) q[3];
sx q[3];
rz(-1.7989649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0380402) q[2];
sx q[2];
rz(-1.2710287) q[2];
sx q[2];
rz(-1.5070149) q[2];
rz(-0.42190894) q[3];
sx q[3];
rz(-0.76015893) q[3];
sx q[3];
rz(-0.71934492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6307395) q[0];
sx q[0];
rz(-0.35950867) q[0];
sx q[0];
rz(1.0605633) q[0];
rz(-3.0781436) q[1];
sx q[1];
rz(-0.63065204) q[1];
sx q[1];
rz(-0.55387703) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3097274) q[0];
sx q[0];
rz(-0.99156443) q[0];
sx q[0];
rz(-2.6088891) q[0];
rz(-pi) q[1];
rz(-2.2138344) q[2];
sx q[2];
rz(-1.0500488) q[2];
sx q[2];
rz(-0.82633229) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0530776) q[1];
sx q[1];
rz(-1.5309528) q[1];
sx q[1];
rz(0.30327176) q[1];
x q[2];
rz(1.2451671) q[3];
sx q[3];
rz(-1.4417329) q[3];
sx q[3];
rz(1.8092524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8744897) q[2];
sx q[2];
rz(-0.59659448) q[2];
sx q[2];
rz(0.60274094) q[2];
rz(-3.0722669) q[3];
sx q[3];
rz(-1.5452789) q[3];
sx q[3];
rz(-0.16665211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.199274) q[0];
sx q[0];
rz(-0.60180226) q[0];
sx q[0];
rz(2.6875575) q[0];
rz(-1.8858689) q[1];
sx q[1];
rz(-2.1818826) q[1];
sx q[1];
rz(1.7032636) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3011264) q[0];
sx q[0];
rz(-0.92341237) q[0];
sx q[0];
rz(0.50917888) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82346114) q[2];
sx q[2];
rz(-2.744906) q[2];
sx q[2];
rz(1.4422063) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0918947) q[1];
sx q[1];
rz(-1.0842853) q[1];
sx q[1];
rz(-0.92384932) q[1];
x q[2];
rz(2.4289729) q[3];
sx q[3];
rz(-1.9639059) q[3];
sx q[3];
rz(-2.676323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.85256514) q[2];
sx q[2];
rz(-1.270741) q[2];
sx q[2];
rz(-1.3555869) q[2];
rz(-1.2191314) q[3];
sx q[3];
rz(-1.4798603) q[3];
sx q[3];
rz(1.5662955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3179625) q[0];
sx q[0];
rz(-0.83331236) q[0];
sx q[0];
rz(1.3939567) q[0];
rz(2.6849003) q[1];
sx q[1];
rz(-0.87867457) q[1];
sx q[1];
rz(-1.3880233) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41075024) q[0];
sx q[0];
rz(-1.9534982) q[0];
sx q[0];
rz(0.14234538) q[0];
rz(-pi) q[1];
x q[1];
rz(0.20095436) q[2];
sx q[2];
rz(-1.6840874) q[2];
sx q[2];
rz(0.73071161) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9364096) q[1];
sx q[1];
rz(-0.35297063) q[1];
sx q[1];
rz(2.754753) q[1];
x q[2];
rz(2.0783992) q[3];
sx q[3];
rz(-1.0143544) q[3];
sx q[3];
rz(2.0991652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8689416) q[2];
sx q[2];
rz(-0.77241263) q[2];
sx q[2];
rz(0.19115494) q[2];
rz(1.8027421) q[3];
sx q[3];
rz(-0.50722417) q[3];
sx q[3];
rz(-0.52201456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.987748) q[0];
sx q[0];
rz(-2.4697883) q[0];
sx q[0];
rz(1.2220569) q[0];
rz(-2.1568495) q[1];
sx q[1];
rz(-0.9877111) q[1];
sx q[1];
rz(-0.15880671) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75750763) q[0];
sx q[0];
rz(-2.473978) q[0];
sx q[0];
rz(-0.13170858) q[0];
rz(-1.5222015) q[2];
sx q[2];
rz(-1.144295) q[2];
sx q[2];
rz(1.0792062) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1013704) q[1];
sx q[1];
rz(-1.985038) q[1];
sx q[1];
rz(2.657402) q[1];
rz(1.7114337) q[3];
sx q[3];
rz(-1.780073) q[3];
sx q[3];
rz(1.9491553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.020393546) q[2];
sx q[2];
rz(-2.8195916) q[2];
sx q[2];
rz(-2.2484692) q[2];
rz(2.761306) q[3];
sx q[3];
rz(-1.2023353) q[3];
sx q[3];
rz(1.7072385) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073008386) q[0];
sx q[0];
rz(-0.46982729) q[0];
sx q[0];
rz(-1.481886) q[0];
rz(2.1549554) q[1];
sx q[1];
rz(-1.6529129) q[1];
sx q[1];
rz(0.557244) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.754234) q[0];
sx q[0];
rz(-0.12436238) q[0];
sx q[0];
rz(-0.4496622) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87961228) q[2];
sx q[2];
rz(-0.45347255) q[2];
sx q[2];
rz(-1.536608) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4367974) q[1];
sx q[1];
rz(-2.3782502) q[1];
sx q[1];
rz(-0.82398681) q[1];
x q[2];
rz(-1.5596296) q[3];
sx q[3];
rz(-1.5876932) q[3];
sx q[3];
rz(-3.0143085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3096932) q[2];
sx q[2];
rz(-2.5688186) q[2];
sx q[2];
rz(2.0834303) q[2];
rz(-1.7582827) q[3];
sx q[3];
rz(-1.0388831) q[3];
sx q[3];
rz(0.96341187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6365373) q[0];
sx q[0];
rz(-1.9872682) q[0];
sx q[0];
rz(2.7099047) q[0];
rz(2.441326) q[1];
sx q[1];
rz(-1.2338748) q[1];
sx q[1];
rz(0.69469992) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85633959) q[0];
sx q[0];
rz(-2.4972557) q[0];
sx q[0];
rz(0.9327234) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0371713) q[2];
sx q[2];
rz(-2.0215144) q[2];
sx q[2];
rz(1.8397651) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8011915) q[1];
sx q[1];
rz(-2.4266198) q[1];
sx q[1];
rz(-2.2919876) q[1];
rz(-pi) q[2];
rz(-2.9630345) q[3];
sx q[3];
rz(-1.6135912) q[3];
sx q[3];
rz(-0.6346441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6557287) q[2];
sx q[2];
rz(-2.8317917) q[2];
sx q[2];
rz(0.41717213) q[2];
rz(-2.396865) q[3];
sx q[3];
rz(-1.797902) q[3];
sx q[3];
rz(-2.1630796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3485296) q[0];
sx q[0];
rz(-1.8804258) q[0];
sx q[0];
rz(-1.6649167) q[0];
rz(-1.9229802) q[1];
sx q[1];
rz(-1.1398133) q[1];
sx q[1];
rz(-2.555991) q[1];
rz(-0.89305604) q[2];
sx q[2];
rz(-2.0231267) q[2];
sx q[2];
rz(-0.33742661) q[2];
rz(-2.811609) q[3];
sx q[3];
rz(-0.97616227) q[3];
sx q[3];
rz(-0.57310692) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
