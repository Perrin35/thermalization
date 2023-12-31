OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.94937593) q[0];
sx q[0];
rz(5.2360143) q[0];
sx q[0];
rz(9.4935023) q[0];
rz(1.7460495) q[1];
sx q[1];
rz(4.6739251) q[1];
sx q[1];
rz(8.2164017) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87586227) q[0];
sx q[0];
rz(-0.10673545) q[0];
sx q[0];
rz(2.021832) q[0];
rz(3/(10*pi)) q[2];
sx q[2];
rz(-0.12636939) q[2];
sx q[2];
rz(-2.7356844) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2627416) q[1];
sx q[1];
rz(-0.95572119) q[1];
sx q[1];
rz(-1.7064852) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0923907) q[3];
sx q[3];
rz(-1.0983101) q[3];
sx q[3];
rz(-2.9977968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8157114) q[2];
sx q[2];
rz(-1.7408966) q[2];
sx q[2];
rz(1.4665843) q[2];
rz(-0.69774929) q[3];
sx q[3];
rz(-2.0402699) q[3];
sx q[3];
rz(2.3944323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0242457) q[0];
sx q[0];
rz(-2.0139366) q[0];
sx q[0];
rz(1.1741937) q[0];
rz(0.17114561) q[1];
sx q[1];
rz(-2.0967963) q[1];
sx q[1];
rz(0.29719621) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3127808) q[0];
sx q[0];
rz(-0.096185616) q[0];
sx q[0];
rz(-2.0975153) q[0];
rz(-pi) q[1];
rz(2.0446288) q[2];
sx q[2];
rz(-0.91544881) q[2];
sx q[2];
rz(1.6608134) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9651523) q[1];
sx q[1];
rz(-2.0808176) q[1];
sx q[1];
rz(0.39217197) q[1];
x q[2];
rz(1.9678715) q[3];
sx q[3];
rz(-2.2928772) q[3];
sx q[3];
rz(1.3305566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.979636) q[2];
sx q[2];
rz(-1.1844144) q[2];
sx q[2];
rz(2.5276108) q[2];
rz(0.87614122) q[3];
sx q[3];
rz(-2.4738779) q[3];
sx q[3];
rz(0.63703018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.0959594) q[0];
sx q[0];
rz(-1.2852083) q[0];
sx q[0];
rz(-2.8339548) q[0];
rz(-0.74854198) q[1];
sx q[1];
rz(-2.8100439) q[1];
sx q[1];
rz(-0.83980733) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6162286) q[0];
sx q[0];
rz(-1.9814241) q[0];
sx q[0];
rz(1.8812268) q[0];
rz(-pi) q[1];
rz(-2.8022021) q[2];
sx q[2];
rz(-1.6600142) q[2];
sx q[2];
rz(0.81775507) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7786583) q[1];
sx q[1];
rz(-2.0211453) q[1];
sx q[1];
rz(1.2308916) q[1];
rz(-pi) q[2];
rz(2.2037376) q[3];
sx q[3];
rz(-1.7114534) q[3];
sx q[3];
rz(-1.7733639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2741189) q[2];
sx q[2];
rz(-1.3511191) q[2];
sx q[2];
rz(-1.8236558) q[2];
rz(1.2157724) q[3];
sx q[3];
rz(-2.7850745) q[3];
sx q[3];
rz(-1.5095476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3777305) q[0];
sx q[0];
rz(-1.3803991) q[0];
sx q[0];
rz(0.44556251) q[0];
rz(0.62082779) q[1];
sx q[1];
rz(-1.8955684) q[1];
sx q[1];
rz(0.96558085) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25528204) q[0];
sx q[0];
rz(-0.83034407) q[0];
sx q[0];
rz(-1.9489954) q[0];
rz(0.031164073) q[2];
sx q[2];
rz(-0.84261299) q[2];
sx q[2];
rz(-2.503501) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.021918745) q[1];
sx q[1];
rz(-1.6026755) q[1];
sx q[1];
rz(-1.8244513) q[1];
rz(-pi) q[2];
rz(1.4569034) q[3];
sx q[3];
rz(-1.4741065) q[3];
sx q[3];
rz(-0.25988042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.821637) q[2];
sx q[2];
rz(-0.76239061) q[2];
sx q[2];
rz(-0.92932534) q[2];
rz(-2.4980513) q[3];
sx q[3];
rz(-2.1054335) q[3];
sx q[3];
rz(-2.1381366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4716361) q[0];
sx q[0];
rz(-1.5595373) q[0];
sx q[0];
rz(-1.2840282) q[0];
rz(-0.28981003) q[1];
sx q[1];
rz(-0.73957864) q[1];
sx q[1];
rz(-1.0481542) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3583583) q[0];
sx q[0];
rz(-2.2336707) q[0];
sx q[0];
rz(2.5213084) q[0];
x q[1];
rz(-0.58446144) q[2];
sx q[2];
rz(-2.2197154) q[2];
sx q[2];
rz(-1.4635758) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2939799) q[1];
sx q[1];
rz(-1.249053) q[1];
sx q[1];
rz(-1.4732248) q[1];
rz(-0.9976451) q[3];
sx q[3];
rz(-0.76612681) q[3];
sx q[3];
rz(-2.0875974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8706878) q[2];
sx q[2];
rz(-2.1257766) q[2];
sx q[2];
rz(2.2407545) q[2];
rz(-1.0926931) q[3];
sx q[3];
rz(-2.1381502) q[3];
sx q[3];
rz(-1.2341011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1822405) q[0];
sx q[0];
rz(-1.2077967) q[0];
sx q[0];
rz(2.545488) q[0];
rz(1.4959363) q[1];
sx q[1];
rz(-2.2438965) q[1];
sx q[1];
rz(1.8966819) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67958528) q[0];
sx q[0];
rz(-1.2512659) q[0];
sx q[0];
rz(2.1035478) q[0];
x q[1];
rz(2.5858324) q[2];
sx q[2];
rz(-2.5854977) q[2];
sx q[2];
rz(-0.13242002) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4391172) q[1];
sx q[1];
rz(-2.2687952) q[1];
sx q[1];
rz(2.4707787) q[1];
rz(0.75002589) q[3];
sx q[3];
rz(-2.8083028) q[3];
sx q[3];
rz(0.94543524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.231455) q[2];
sx q[2];
rz(-0.69880501) q[2];
sx q[2];
rz(-2.9166252) q[2];
rz(3.0531626) q[3];
sx q[3];
rz(-1.4513353) q[3];
sx q[3];
rz(-2.6627873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46463075) q[0];
sx q[0];
rz(-1.2372274) q[0];
sx q[0];
rz(2.8616469) q[0];
rz(-1.6784558) q[1];
sx q[1];
rz(-1.2660374) q[1];
sx q[1];
rz(2.8889012) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42851105) q[0];
sx q[0];
rz(-1.5230852) q[0];
sx q[0];
rz(-1.2776572) q[0];
x q[1];
rz(1.5451317) q[2];
sx q[2];
rz(-2.5070094) q[2];
sx q[2];
rz(-2.6013825) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6934769) q[1];
sx q[1];
rz(-2.2709393) q[1];
sx q[1];
rz(2.4561988) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85315506) q[3];
sx q[3];
rz(-0.95082885) q[3];
sx q[3];
rz(-2.6593047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.32020405) q[2];
sx q[2];
rz(-2.2479222) q[2];
sx q[2];
rz(-1.5318711) q[2];
rz(-1.1931217) q[3];
sx q[3];
rz(-2.2333998) q[3];
sx q[3];
rz(-1.0866722) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10483345) q[0];
sx q[0];
rz(-1.7585254) q[0];
sx q[0];
rz(-1.6554792) q[0];
rz(-2.8727818) q[1];
sx q[1];
rz(-1.1299645) q[1];
sx q[1];
rz(-2.862646) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19764087) q[0];
sx q[0];
rz(-0.949172) q[0];
sx q[0];
rz(0.62244121) q[0];
rz(-pi) q[1];
rz(-2.0729162) q[2];
sx q[2];
rz(-1.4197822) q[2];
sx q[2];
rz(-0.47889027) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.61356269) q[1];
sx q[1];
rz(-0.87606214) q[1];
sx q[1];
rz(-1.0747521) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99854462) q[3];
sx q[3];
rz(-1.4456985) q[3];
sx q[3];
rz(0.78748122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3387317) q[2];
sx q[2];
rz(-1.677745) q[2];
sx q[2];
rz(1.5489244) q[2];
rz(1.0507978) q[3];
sx q[3];
rz(-1.934634) q[3];
sx q[3];
rz(-1.9416434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6779697) q[0];
sx q[0];
rz(-0.93358731) q[0];
sx q[0];
rz(1.2532225) q[0];
rz(-0.62250096) q[1];
sx q[1];
rz(-1.6815192) q[1];
sx q[1];
rz(1.1463096) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63908731) q[0];
sx q[0];
rz(-0.82477942) q[0];
sx q[0];
rz(-1.4750925) q[0];
x q[1];
rz(-0.81252819) q[2];
sx q[2];
rz(-0.84156893) q[2];
sx q[2];
rz(1.3844045) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.031530023) q[1];
sx q[1];
rz(-1.6441802) q[1];
sx q[1];
rz(0.88295464) q[1];
x q[2];
rz(0.63486926) q[3];
sx q[3];
rz(-1.2253309) q[3];
sx q[3];
rz(2.8431901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.15381947) q[2];
sx q[2];
rz(-1.035707) q[2];
sx q[2];
rz(-1.0158319) q[2];
rz(0.9097957) q[3];
sx q[3];
rz(-0.67009059) q[3];
sx q[3];
rz(-0.36809665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1290454) q[0];
sx q[0];
rz(-2.5279901) q[0];
sx q[0];
rz(3.1273499) q[0];
rz(-0.8447389) q[1];
sx q[1];
rz(-2.1886487) q[1];
sx q[1];
rz(1.3815809) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89109126) q[0];
sx q[0];
rz(-1.7133461) q[0];
sx q[0];
rz(0.010682627) q[0];
rz(-pi) q[1];
rz(-0.064247473) q[2];
sx q[2];
rz(-0.99284961) q[2];
sx q[2];
rz(-2.3648175) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6784918) q[1];
sx q[1];
rz(-1.4271724) q[1];
sx q[1];
rz(1.0113869) q[1];
x q[2];
rz(-0.53346177) q[3];
sx q[3];
rz(-1.3812314) q[3];
sx q[3];
rz(1.5272527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0766729) q[2];
sx q[2];
rz(-1.939247) q[2];
sx q[2];
rz(-0.98999611) q[2];
rz(-2.2475217) q[3];
sx q[3];
rz(-2.6529513) q[3];
sx q[3];
rz(2.9161684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0902696) q[0];
sx q[0];
rz(-2.4933503) q[0];
sx q[0];
rz(-1.1052263) q[0];
rz(-1.3399711) q[1];
sx q[1];
rz(-2.5201288) q[1];
sx q[1];
rz(-2.7531243) q[1];
rz(-2.4253035) q[2];
sx q[2];
rz(-1.204797) q[2];
sx q[2];
rz(2.7439678) q[2];
rz(-0.43511012) q[3];
sx q[3];
rz(-0.995244) q[3];
sx q[3];
rz(1.4959195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
