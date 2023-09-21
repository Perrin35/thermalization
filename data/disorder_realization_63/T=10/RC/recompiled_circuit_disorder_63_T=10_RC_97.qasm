OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9665943) q[0];
sx q[0];
rz(-2.7881665) q[0];
sx q[0];
rz(2.0768291) q[0];
rz(-2.3454173) q[1];
sx q[1];
rz(-1.2086955) q[1];
sx q[1];
rz(2.605521) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.160643) q[0];
sx q[0];
rz(-1.0796709) q[0];
sx q[0];
rz(0.26382291) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0739003) q[2];
sx q[2];
rz(-2.2798685) q[2];
sx q[2];
rz(-0.46082218) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.55878996) q[1];
sx q[1];
rz(-1.7585808) q[1];
sx q[1];
rz(-0.59273984) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0814813) q[3];
sx q[3];
rz(-2.2430674) q[3];
sx q[3];
rz(2.1171452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4094231) q[2];
sx q[2];
rz(-2.7004341) q[2];
sx q[2];
rz(-1.0268964) q[2];
rz(-2.8895767) q[3];
sx q[3];
rz(-1.9988632) q[3];
sx q[3];
rz(-0.86565971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.058763) q[0];
sx q[0];
rz(-0.34794647) q[0];
sx q[0];
rz(1.7513562) q[0];
rz(-2.305796) q[1];
sx q[1];
rz(-0.73671571) q[1];
sx q[1];
rz(-0.70835152) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7824552) q[0];
sx q[0];
rz(-1.8498427) q[0];
sx q[0];
rz(-2.7611087) q[0];
x q[1];
rz(-2.1721241) q[2];
sx q[2];
rz(-1.4496441) q[2];
sx q[2];
rz(3.0960992) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.12304141) q[1];
sx q[1];
rz(-2.0611144) q[1];
sx q[1];
rz(0.80177387) q[1];
rz(-0.80531081) q[3];
sx q[3];
rz(-0.66262965) q[3];
sx q[3];
rz(-0.90028541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0248802) q[2];
sx q[2];
rz(-2.34989) q[2];
sx q[2];
rz(2.8498245) q[2];
rz(-0.10270384) q[3];
sx q[3];
rz(-1.4029968) q[3];
sx q[3];
rz(1.6170988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3290688) q[0];
sx q[0];
rz(-0.085443184) q[0];
sx q[0];
rz(-0.30971757) q[0];
rz(1.6614871) q[1];
sx q[1];
rz(-1.3269576) q[1];
sx q[1];
rz(-2.5699239) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11274352) q[0];
sx q[0];
rz(-1.0380121) q[0];
sx q[0];
rz(1.4677731) q[0];
rz(-pi) q[1];
rz(1.13582) q[2];
sx q[2];
rz(-1.8990714) q[2];
sx q[2];
rz(-0.091718397) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.09517041) q[1];
sx q[1];
rz(-1.4383525) q[1];
sx q[1];
rz(2.4080647) q[1];
x q[2];
rz(0.27922697) q[3];
sx q[3];
rz(-2.673827) q[3];
sx q[3];
rz(1.2665018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.16033515) q[2];
sx q[2];
rz(-0.91032878) q[2];
sx q[2];
rz(-0.70880115) q[2];
rz(-2.839084) q[3];
sx q[3];
rz(-1.6995647) q[3];
sx q[3];
rz(-1.1423053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70808327) q[0];
sx q[0];
rz(-2.0286562) q[0];
sx q[0];
rz(-2.3420912) q[0];
rz(-3.0917621) q[1];
sx q[1];
rz(-0.89490503) q[1];
sx q[1];
rz(0.18049151) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40515306) q[0];
sx q[0];
rz(-0.3973876) q[0];
sx q[0];
rz(-0.23657628) q[0];
rz(-pi) q[1];
rz(-1.7454342) q[2];
sx q[2];
rz(-1.9348113) q[2];
sx q[2];
rz(0.82496914) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.71112878) q[1];
sx q[1];
rz(-0.34090484) q[1];
sx q[1];
rz(-2.7845963) q[1];
rz(-1.5482076) q[3];
sx q[3];
rz(-2.6456607) q[3];
sx q[3];
rz(3.0043234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.089347) q[2];
sx q[2];
rz(-1.1648488) q[2];
sx q[2];
rz(1.583741) q[2];
rz(-1.7662988) q[3];
sx q[3];
rz(-1.3170653) q[3];
sx q[3];
rz(-2.2089675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76064008) q[0];
sx q[0];
rz(-0.33712688) q[0];
sx q[0];
rz(2.1550762) q[0];
rz(-1.1622693) q[1];
sx q[1];
rz(-1.9245851) q[1];
sx q[1];
rz(0.25156897) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0662721) q[0];
sx q[0];
rz(-2.2499654) q[0];
sx q[0];
rz(-2.9156296) q[0];
rz(0.45995633) q[2];
sx q[2];
rz(-1.492968) q[2];
sx q[2];
rz(-2.3098582) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5800036) q[1];
sx q[1];
rz(-1.5981234) q[1];
sx q[1];
rz(-0.1954397) q[1];
x q[2];
rz(-2.2734363) q[3];
sx q[3];
rz(-1.1690825) q[3];
sx q[3];
rz(0.86740869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.183737) q[2];
sx q[2];
rz(-0.541406) q[2];
sx q[2];
rz(-2.1110558) q[2];
rz(1.41097) q[3];
sx q[3];
rz(-2.1822699) q[3];
sx q[3];
rz(2.5103536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3866766) q[0];
sx q[0];
rz(-0.47575352) q[0];
sx q[0];
rz(2.2914698) q[0];
rz(-1.1823581) q[1];
sx q[1];
rz(-1.2344924) q[1];
sx q[1];
rz(0.39302557) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.128199) q[0];
sx q[0];
rz(-1.4067187) q[0];
sx q[0];
rz(-1.2584932) q[0];
rz(2.9057301) q[2];
sx q[2];
rz(-2.3013407) q[2];
sx q[2];
rz(-1.1299396) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.46950754) q[1];
sx q[1];
rz(-2.2219594) q[1];
sx q[1];
rz(0.0068411946) q[1];
x q[2];
rz(2.399746) q[3];
sx q[3];
rz(-0.39374712) q[3];
sx q[3];
rz(1.1175931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7810016) q[2];
sx q[2];
rz(-2.0157308) q[2];
sx q[2];
rz(2.1684516) q[2];
rz(-0.94758236) q[3];
sx q[3];
rz(-2.0411453) q[3];
sx q[3];
rz(-1.1943641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4457552) q[0];
sx q[0];
rz(-0.42036244) q[0];
sx q[0];
rz(1.5234891) q[0];
rz(0.37480005) q[1];
sx q[1];
rz(-1.8811036) q[1];
sx q[1];
rz(-0.0059676776) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80124827) q[0];
sx q[0];
rz(-2.4140515) q[0];
sx q[0];
rz(-1.9968541) q[0];
rz(-pi) q[1];
rz(2.3890424) q[2];
sx q[2];
rz(-1.9117022) q[2];
sx q[2];
rz(0.37170751) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.55885) q[1];
sx q[1];
rz(-0.55816459) q[1];
sx q[1];
rz(-1.4504257) q[1];
rz(-1.5238477) q[3];
sx q[3];
rz(-1.1991683) q[3];
sx q[3];
rz(0.23577984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5696047) q[2];
sx q[2];
rz(-0.54509744) q[2];
sx q[2];
rz(0.19006426) q[2];
rz(-2.7251785) q[3];
sx q[3];
rz(-0.9581241) q[3];
sx q[3];
rz(0.73474187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1104601) q[0];
sx q[0];
rz(-2.0955595) q[0];
sx q[0];
rz(2.6053612) q[0];
rz(0.93332943) q[1];
sx q[1];
rz(-1.5678762) q[1];
sx q[1];
rz(-0.94820625) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1873916) q[0];
sx q[0];
rz(-1.3905977) q[0];
sx q[0];
rz(-0.76822922) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0529222) q[2];
sx q[2];
rz(-1.1793943) q[2];
sx q[2];
rz(0.97757593) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1369689) q[1];
sx q[1];
rz(-1.3490632) q[1];
sx q[1];
rz(-1.6349413) q[1];
x q[2];
rz(-0.20235297) q[3];
sx q[3];
rz(-0.36365299) q[3];
sx q[3];
rz(2.3285151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5006717) q[2];
sx q[2];
rz(-1.6225953) q[2];
sx q[2];
rz(2.9013157) q[2];
rz(-0.37929532) q[3];
sx q[3];
rz(-1.0255739) q[3];
sx q[3];
rz(1.7933638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76294476) q[0];
sx q[0];
rz(-0.33319107) q[0];
sx q[0];
rz(2.8906524) q[0];
rz(0.25302408) q[1];
sx q[1];
rz(-1.3809985) q[1];
sx q[1];
rz(-2.7889263) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0371373) q[0];
sx q[0];
rz(-1.4482575) q[0];
sx q[0];
rz(-0.11549581) q[0];
rz(-1.4264832) q[2];
sx q[2];
rz(-1.2735954) q[2];
sx q[2];
rz(-0.65629634) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.36665146) q[1];
sx q[1];
rz(-2.1244441) q[1];
sx q[1];
rz(1.8371546) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3064763) q[3];
sx q[3];
rz(-1.9559238) q[3];
sx q[3];
rz(1.1564099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.76715604) q[2];
sx q[2];
rz(-2.0220951) q[2];
sx q[2];
rz(-0.68332589) q[2];
rz(-1.654489) q[3];
sx q[3];
rz(-0.43928248) q[3];
sx q[3];
rz(1.1504415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7651354) q[0];
sx q[0];
rz(-2.6177804) q[0];
sx q[0];
rz(1.3002243) q[0];
rz(2.4328649) q[1];
sx q[1];
rz(-2.6649902) q[1];
sx q[1];
rz(2.2050819) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47637981) q[0];
sx q[0];
rz(-1.8139651) q[0];
sx q[0];
rz(2.5961844) q[0];
rz(0.18239393) q[2];
sx q[2];
rz(-1.8511904) q[2];
sx q[2];
rz(2.3390351) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9771802) q[1];
sx q[1];
rz(-0.74294786) q[1];
sx q[1];
rz(1.8650024) q[1];
rz(-pi) q[2];
rz(1.1141206) q[3];
sx q[3];
rz(-0.63318397) q[3];
sx q[3];
rz(-3.0999108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.54291723) q[2];
sx q[2];
rz(-2.8024555) q[2];
sx q[2];
rz(-0.014952095) q[2];
rz(2.9144918) q[3];
sx q[3];
rz(-0.98277074) q[3];
sx q[3];
rz(-1.8792413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99075714) q[0];
sx q[0];
rz(-0.80934722) q[0];
sx q[0];
rz(1.9833175) q[0];
rz(1.5630209) q[1];
sx q[1];
rz(-2.38588) q[1];
sx q[1];
rz(-0.26185782) q[1];
rz(-1.0583393) q[2];
sx q[2];
rz(-2.8556311) q[2];
sx q[2];
rz(-0.43559504) q[2];
rz(-2.6111504) q[3];
sx q[3];
rz(-1.3417288) q[3];
sx q[3];
rz(-0.19977278) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];