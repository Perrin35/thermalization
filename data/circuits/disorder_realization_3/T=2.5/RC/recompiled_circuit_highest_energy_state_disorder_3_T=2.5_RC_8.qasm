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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4473576) q[0];
sx q[0];
rz(-1.5120872) q[0];
sx q[0];
rz(-2.319058) q[0];
x q[1];
rz(-0.069434631) q[2];
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
rz(pi/2) q[0];
sx q[0];
rz(0.82935059) q[1];
sx q[1];
rz(-1.568746) q[1];
sx q[1];
rz(-0.41872989) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3326268) q[3];
sx q[3];
rz(-2.3936317) q[3];
sx q[3];
rz(0.35424074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1623666) q[2];
sx q[2];
rz(-1.6567076) q[2];
sx q[2];
rz(-1.4244728) q[2];
rz(-3.034397) q[3];
sx q[3];
rz(-0.19459477) q[3];
sx q[3];
rz(-2.181982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87218881) q[0];
sx q[0];
rz(-2.2087681) q[0];
sx q[0];
rz(-1.8316356) q[0];
rz(-0.45920363) q[1];
sx q[1];
rz(-1.2745067) q[1];
sx q[1];
rz(-0.31712636) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91823365) q[0];
sx q[0];
rz(-0.14592136) q[0];
sx q[0];
rz(0.33182611) q[0];
rz(-3.0917815) q[2];
sx q[2];
rz(-1.5684109) q[2];
sx q[2];
rz(-2.002169) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0786653) q[1];
sx q[1];
rz(-2.1172273) q[1];
sx q[1];
rz(1.648047) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22038118) q[3];
sx q[3];
rz(-2.3733449) q[3];
sx q[3];
rz(0.11742442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.7121048) q[2];
sx q[2];
rz(-1.2965344) q[2];
sx q[2];
rz(0.8075766) q[2];
rz(0.56337041) q[3];
sx q[3];
rz(-2.6475776) q[3];
sx q[3];
rz(2.0886683) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5092369) q[0];
sx q[0];
rz(-0.18993264) q[0];
sx q[0];
rz(-0.83455363) q[0];
rz(-2.6368311) q[1];
sx q[1];
rz(-2.7752462) q[1];
sx q[1];
rz(1.6827481) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.014011009) q[0];
sx q[0];
rz(-2.547285) q[0];
sx q[0];
rz(2.7796714) q[0];
rz(2.3332852) q[2];
sx q[2];
rz(-0.9918074) q[2];
sx q[2];
rz(2.4579687) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9533938) q[1];
sx q[1];
rz(-2.2911304) q[1];
sx q[1];
rz(1.6745643) q[1];
rz(-pi) q[2];
rz(0.80648544) q[3];
sx q[3];
rz(-2.0098662) q[3];
sx q[3];
rz(2.0937378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.10996455) q[2];
sx q[2];
rz(-1.6692903) q[2];
sx q[2];
rz(-2.7785981) q[2];
rz(1.1686769) q[3];
sx q[3];
rz(-2.3790338) q[3];
sx q[3];
rz(2.3762083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4204243) q[0];
sx q[0];
rz(-0.52614251) q[0];
sx q[0];
rz(2.5508733) q[0];
rz(2.4144454) q[1];
sx q[1];
rz(-0.37494451) q[1];
sx q[1];
rz(2.4730543) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97897935) q[0];
sx q[0];
rz(-3.0501462) q[0];
sx q[0];
rz(-2.0434815) q[0];
rz(-pi) q[1];
rz(-2.6538969) q[2];
sx q[2];
rz(-0.82789153) q[2];
sx q[2];
rz(2.6726892) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4674899) q[1];
sx q[1];
rz(-2.2558442) q[1];
sx q[1];
rz(-2.5193307) q[1];
rz(-0.7410243) q[3];
sx q[3];
rz(-1.1644568) q[3];
sx q[3];
rz(-2.5056553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.10355243) q[2];
sx q[2];
rz(-1.870564) q[2];
sx q[2];
rz(1.5070149) q[2];
rz(-2.7196837) q[3];
sx q[3];
rz(-2.3814337) q[3];
sx q[3];
rz(2.4222477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6307395) q[0];
sx q[0];
rz(-0.35950867) q[0];
sx q[0];
rz(2.0810293) q[0];
rz(0.063449055) q[1];
sx q[1];
rz(-0.63065204) q[1];
sx q[1];
rz(2.5877156) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3097274) q[0];
sx q[0];
rz(-2.1500282) q[0];
sx q[0];
rz(0.53270355) q[0];
rz(-0.62184288) q[2];
sx q[2];
rz(-2.1178195) q[2];
sx q[2];
rz(-2.0402997) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.088515048) q[1];
sx q[1];
rz(-1.6106399) q[1];
sx q[1];
rz(-0.30327176) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13613607) q[3];
sx q[3];
rz(-1.247974) q[3];
sx q[3];
rz(-2.8597067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.267103) q[2];
sx q[2];
rz(-0.59659448) q[2];
sx q[2];
rz(0.60274094) q[2];
rz(0.069325773) q[3];
sx q[3];
rz(-1.5452789) q[3];
sx q[3];
rz(2.9749405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.199274) q[0];
sx q[0];
rz(-2.5397904) q[0];
sx q[0];
rz(-0.45403516) q[0];
rz(-1.2557238) q[1];
sx q[1];
rz(-0.95971003) q[1];
sx q[1];
rz(-1.4383291) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8404663) q[0];
sx q[0];
rz(-2.2181803) q[0];
sx q[0];
rz(-0.50917888) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27737646) q[2];
sx q[2];
rz(-1.8581333) q[2];
sx q[2];
rz(2.4871021) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.033129902) q[1];
sx q[1];
rz(-2.3537043) q[1];
sx q[1];
rz(2.2910396) q[1];
rz(-2.4289729) q[3];
sx q[3];
rz(-1.9639059) q[3];
sx q[3];
rz(2.676323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2890275) q[2];
sx q[2];
rz(-1.8708517) q[2];
sx q[2];
rz(1.7860058) q[2];
rz(1.2191314) q[3];
sx q[3];
rz(-1.6617323) q[3];
sx q[3];
rz(-1.5752972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3179625) q[0];
sx q[0];
rz(-2.3082803) q[0];
sx q[0];
rz(-1.3939567) q[0];
rz(-0.45669237) q[1];
sx q[1];
rz(-2.2629181) q[1];
sx q[1];
rz(1.3880233) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41075024) q[0];
sx q[0];
rz(-1.9534982) q[0];
sx q[0];
rz(0.14234538) q[0];
x q[1];
rz(-0.20095436) q[2];
sx q[2];
rz(-1.4575053) q[2];
sx q[2];
rz(0.73071161) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4108666) q[1];
sx q[1];
rz(-1.4400088) q[1];
sx q[1];
rz(-2.8128002) q[1];
rz(-2.5230266) q[3];
sx q[3];
rz(-1.1453218) q[3];
sx q[3];
rz(2.8989603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8689416) q[2];
sx q[2];
rz(-0.77241263) q[2];
sx q[2];
rz(0.19115494) q[2];
rz(1.3388505) q[3];
sx q[3];
rz(-0.50722417) q[3];
sx q[3];
rz(0.52201456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1538447) q[0];
sx q[0];
rz(-0.67180434) q[0];
sx q[0];
rz(-1.2220569) q[0];
rz(-0.98474312) q[1];
sx q[1];
rz(-2.1538815) q[1];
sx q[1];
rz(2.9827859) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5903871) q[0];
sx q[0];
rz(-2.231601) q[0];
sx q[0];
rz(-1.673965) q[0];
rz(-pi) q[1];
rz(-3.0350948) q[2];
sx q[2];
rz(-2.7125008) q[2];
sx q[2];
rz(0.96218357) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1837801) q[1];
sx q[1];
rz(-2.5152825) q[1];
sx q[1];
rz(-0.7568936) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7114337) q[3];
sx q[3];
rz(-1.3615196) q[3];
sx q[3];
rz(1.9491553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.020393546) q[2];
sx q[2];
rz(-0.32200107) q[2];
sx q[2];
rz(0.89312345) q[2];
rz(-0.3802866) q[3];
sx q[3];
rz(-1.2023353) q[3];
sx q[3];
rz(-1.4343542) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073008386) q[0];
sx q[0];
rz(-2.6717654) q[0];
sx q[0];
rz(-1.481886) q[0];
rz(-2.1549554) q[1];
sx q[1];
rz(-1.4886798) q[1];
sx q[1];
rz(0.557244) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6300723) q[0];
sx q[0];
rz(-1.5168539) q[0];
sx q[0];
rz(3.0294837) q[0];
rz(-pi) q[1];
rz(2.2619804) q[2];
sx q[2];
rz(-0.45347255) q[2];
sx q[2];
rz(-1.6049847) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4367974) q[1];
sx q[1];
rz(-2.3782502) q[1];
sx q[1];
rz(-0.82398681) q[1];
rz(-pi) q[2];
x q[2];
rz(2.557665) q[3];
sx q[3];
rz(-0.020253094) q[3];
sx q[3];
rz(2.4302866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8318994) q[2];
sx q[2];
rz(-2.5688186) q[2];
sx q[2];
rz(-1.0581623) q[2];
rz(-1.38331) q[3];
sx q[3];
rz(-1.0388831) q[3];
sx q[3];
rz(2.1781808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50505534) q[0];
sx q[0];
rz(-1.9872682) q[0];
sx q[0];
rz(2.7099047) q[0];
rz(2.441326) q[1];
sx q[1];
rz(-1.9077178) q[1];
sx q[1];
rz(2.4468927) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10852495) q[0];
sx q[0];
rz(-2.0742848) q[0];
sx q[0];
rz(-0.42078542) q[0];
rz(-pi) q[1];
rz(-1.117933) q[2];
sx q[2];
rz(-1.6647571) q[2];
sx q[2];
rz(2.8270023) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8164639) q[1];
sx q[1];
rz(-1.123114) q[1];
sx q[1];
rz(0.992985) q[1];
rz(-pi) q[2];
rz(1.6142817) q[3];
sx q[3];
rz(-1.3924034) q[3];
sx q[3];
rz(-2.1977193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.485864) q[2];
sx q[2];
rz(-2.8317917) q[2];
sx q[2];
rz(-0.41717213) q[2];
rz(0.74472767) q[3];
sx q[3];
rz(-1.3436907) q[3];
sx q[3];
rz(-0.97851306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7930631) q[0];
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
rz(2.1914235) q[3];
sx q[3];
rz(-1.8425377) q[3];
sx q[3];
rz(0.80811926) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
