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
rz(0.87297451) q[0];
sx q[0];
rz(-0.31384808) q[0];
sx q[0];
rz(2.0883972) q[0];
rz(1.6246417) q[1];
sx q[1];
rz(2.8487974) q[1];
sx q[1];
rz(8.4881633) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7760843) q[0];
sx q[0];
rz(-2.0788686) q[0];
sx q[0];
rz(-3.1209141) q[0];
rz(-pi) q[1];
rz(-2.6153271) q[2];
sx q[2];
rz(-0.91134763) q[2];
sx q[2];
rz(1.1629205) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0055252) q[1];
sx q[1];
rz(-0.64856358) q[1];
sx q[1];
rz(-0.4737718) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7009007) q[3];
sx q[3];
rz(-0.42817906) q[3];
sx q[3];
rz(2.6994914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8006353) q[2];
sx q[2];
rz(-0.95658797) q[2];
sx q[2];
rz(-2.0802278) q[2];
rz(-0.64570767) q[3];
sx q[3];
rz(-0.3568477) q[3];
sx q[3];
rz(2.1319353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89479947) q[0];
sx q[0];
rz(-2.6428887) q[0];
sx q[0];
rz(0.45951581) q[0];
rz(1.7543606) q[1];
sx q[1];
rz(-0.51001716) q[1];
sx q[1];
rz(2.0769108) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.014039847) q[0];
sx q[0];
rz(-1.8888374) q[0];
sx q[0];
rz(-1.4038588) q[0];
x q[1];
rz(0.063195066) q[2];
sx q[2];
rz(-2.8738351) q[2];
sx q[2];
rz(-1.7187985) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.766651) q[1];
sx q[1];
rz(-1.8901575) q[1];
sx q[1];
rz(-1.6691895) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2336938) q[3];
sx q[3];
rz(-1.3114942) q[3];
sx q[3];
rz(1.8725267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3620944) q[2];
sx q[2];
rz(-0.67216122) q[2];
sx q[2];
rz(-2.6008495) q[2];
rz(-2.8703459) q[3];
sx q[3];
rz(-1.8999148) q[3];
sx q[3];
rz(2.0103683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3830477) q[0];
sx q[0];
rz(-2.8309412) q[0];
sx q[0];
rz(0.3824105) q[0];
rz(-2.8657939) q[1];
sx q[1];
rz(-2.661992) q[1];
sx q[1];
rz(-2.3013505) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0226757) q[0];
sx q[0];
rz(-1.9249601) q[0];
sx q[0];
rz(0.35665956) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98767583) q[2];
sx q[2];
rz(-2.7842369) q[2];
sx q[2];
rz(-2.5492956) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0178595) q[1];
sx q[1];
rz(-0.85461939) q[1];
sx q[1];
rz(2.2804428) q[1];
rz(-pi) q[2];
rz(-0.48583416) q[3];
sx q[3];
rz(-1.6260901) q[3];
sx q[3];
rz(3.0818617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8648839) q[2];
sx q[2];
rz(-2.101892) q[2];
sx q[2];
rz(3.0806105) q[2];
rz(-1.3649155) q[3];
sx q[3];
rz(-2.958332) q[3];
sx q[3];
rz(1.1158367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10114305) q[0];
sx q[0];
rz(-2.1231066) q[0];
sx q[0];
rz(2.2818991) q[0];
rz(-2.1172093) q[1];
sx q[1];
rz(-2.5442217) q[1];
sx q[1];
rz(0.76622564) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9767155) q[0];
sx q[0];
rz(-0.055792965) q[0];
sx q[0];
rz(-0.68875046) q[0];
rz(1.8623975) q[2];
sx q[2];
rz(-0.5480786) q[2];
sx q[2];
rz(1.6432135) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1658881) q[1];
sx q[1];
rz(-0.59671003) q[1];
sx q[1];
rz(-2.6590682) q[1];
x q[2];
rz(-0.071841876) q[3];
sx q[3];
rz(-2.2719682) q[3];
sx q[3];
rz(1.5909242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.031294558) q[2];
sx q[2];
rz(-2.2790907) q[2];
sx q[2];
rz(-0.026570126) q[2];
rz(2.8734112) q[3];
sx q[3];
rz(-0.53332204) q[3];
sx q[3];
rz(-2.4466483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-0.37102315) q[0];
sx q[0];
rz(-0.71582782) q[0];
sx q[0];
rz(-2.7243966) q[0];
rz(-0.5873276) q[1];
sx q[1];
rz(-0.53343499) q[1];
sx q[1];
rz(0.74648285) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70770812) q[0];
sx q[0];
rz(-1.2217772) q[0];
sx q[0];
rz(1.2322578) q[0];
x q[1];
rz(2.3801489) q[2];
sx q[2];
rz(-2.1816744) q[2];
sx q[2];
rz(2.2100984) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.58140671) q[1];
sx q[1];
rz(-0.65196067) q[1];
sx q[1];
rz(-0.65309872) q[1];
rz(0.77961214) q[3];
sx q[3];
rz(-2.5031708) q[3];
sx q[3];
rz(2.8299931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8314961) q[2];
sx q[2];
rz(-0.3336755) q[2];
sx q[2];
rz(0.92158544) q[2];
rz(2.0690252) q[3];
sx q[3];
rz(-1.8662063) q[3];
sx q[3];
rz(2.5785562) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6801341) q[0];
sx q[0];
rz(-0.38668329) q[0];
sx q[0];
rz(0.65163809) q[0];
rz(-0.98546511) q[1];
sx q[1];
rz(-2.5990504) q[1];
sx q[1];
rz(-0.14269565) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14772077) q[0];
sx q[0];
rz(-1.8447478) q[0];
sx q[0];
rz(3.0763294) q[0];
x q[1];
rz(-2.8593117) q[2];
sx q[2];
rz(-2.9476894) q[2];
sx q[2];
rz(-1.0218203) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7276926) q[1];
sx q[1];
rz(-1.6928612) q[1];
sx q[1];
rz(-2.6660566) q[1];
rz(0.13655314) q[3];
sx q[3];
rz(-1.1205744) q[3];
sx q[3];
rz(-3.0409851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5695213) q[2];
sx q[2];
rz(-2.5321952) q[2];
sx q[2];
rz(0.93192464) q[2];
rz(2.7052687) q[3];
sx q[3];
rz(-0.47567979) q[3];
sx q[3];
rz(-1.9816192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6465004) q[0];
sx q[0];
rz(-2.2252872) q[0];
sx q[0];
rz(1.5402933) q[0];
rz(-2.1401999) q[1];
sx q[1];
rz(-1.5903558) q[1];
sx q[1];
rz(2.1839949) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1100028) q[0];
sx q[0];
rz(-0.45781198) q[0];
sx q[0];
rz(-0.97719595) q[0];
rz(1.4740491) q[2];
sx q[2];
rz(-1.3203586) q[2];
sx q[2];
rz(-1.4055523) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0066788) q[1];
sx q[1];
rz(-1.5777262) q[1];
sx q[1];
rz(-2.0459941) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0882629) q[3];
sx q[3];
rz(-2.5825273) q[3];
sx q[3];
rz(-2.6757765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6340948) q[2];
sx q[2];
rz(-1.8562506) q[2];
sx q[2];
rz(1.1320587) q[2];
rz(-0.13103983) q[3];
sx q[3];
rz(-1.1538426) q[3];
sx q[3];
rz(1.2426144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2496654) q[0];
sx q[0];
rz(-0.51814336) q[0];
sx q[0];
rz(3.0666572) q[0];
rz(-0.19649188) q[1];
sx q[1];
rz(-2.7939929) q[1];
sx q[1];
rz(2.4620655) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9411036) q[0];
sx q[0];
rz(-1.3692432) q[0];
sx q[0];
rz(-2.5747712) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9747693) q[2];
sx q[2];
rz(-1.543664) q[2];
sx q[2];
rz(2.2441442) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.841015) q[1];
sx q[1];
rz(-1.1493582) q[1];
sx q[1];
rz(2.8309062) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1155247) q[3];
sx q[3];
rz(-0.51122087) q[3];
sx q[3];
rz(-2.8554684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3139265) q[2];
sx q[2];
rz(-0.9845261) q[2];
sx q[2];
rz(2.0619681) q[2];
rz(0.45352724) q[3];
sx q[3];
rz(-0.87117666) q[3];
sx q[3];
rz(-2.1760904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.51665783) q[0];
sx q[0];
rz(-2.9624532) q[0];
sx q[0];
rz(3.0058885) q[0];
rz(-0.16595674) q[1];
sx q[1];
rz(-2.143492) q[1];
sx q[1];
rz(0.96806324) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57763273) q[0];
sx q[0];
rz(-1.5819307) q[0];
sx q[0];
rz(2.0861113) q[0];
rz(2.9918475) q[2];
sx q[2];
rz(-2.1080842) q[2];
sx q[2];
rz(1.4531236) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1727454) q[1];
sx q[1];
rz(-2.0683534) q[1];
sx q[1];
rz(-1.8109591) q[1];
x q[2];
rz(-2.5510213) q[3];
sx q[3];
rz(-2.6988523) q[3];
sx q[3];
rz(2.8273432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.82219899) q[2];
sx q[2];
rz(-2.2182756) q[2];
sx q[2];
rz(2.609002) q[2];
rz(2.3432664) q[3];
sx q[3];
rz(-2.7440378) q[3];
sx q[3];
rz(-0.09440162) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7867197) q[0];
sx q[0];
rz(-2.4479471) q[0];
sx q[0];
rz(-2.4548446) q[0];
rz(-1.2314388) q[1];
sx q[1];
rz(-1.310692) q[1];
sx q[1];
rz(3.0795857) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69294482) q[0];
sx q[0];
rz(-0.8668859) q[0];
sx q[0];
rz(1.6567898) q[0];
rz(-pi) q[1];
rz(-1.643963) q[2];
sx q[2];
rz(-1.3806297) q[2];
sx q[2];
rz(-2.721618) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0289284) q[1];
sx q[1];
rz(-0.24953574) q[1];
sx q[1];
rz(0.94503944) q[1];
x q[2];
rz(2.8436321) q[3];
sx q[3];
rz(-1.584867) q[3];
sx q[3];
rz(1.640682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9318781) q[2];
sx q[2];
rz(-0.96499062) q[2];
sx q[2];
rz(-0.17267257) q[2];
rz(-2.4058345) q[3];
sx q[3];
rz(-0.57307214) q[3];
sx q[3];
rz(0.47634038) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.944735) q[0];
sx q[0];
rz(-1.6898962) q[0];
sx q[0];
rz(1.663399) q[0];
rz(2.2810777) q[1];
sx q[1];
rz(-1.1970701) q[1];
sx q[1];
rz(1.7477716) q[1];
rz(1.6780268) q[2];
sx q[2];
rz(-1.4489531) q[2];
sx q[2];
rz(1.0482241) q[2];
rz(2.2279578) q[3];
sx q[3];
rz(-1.1362094) q[3];
sx q[3];
rz(2.1528578) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
