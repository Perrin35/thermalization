OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6948833) q[0];
sx q[0];
rz(-1.0816242) q[0];
sx q[0];
rz(2.4021436) q[0];
rz(2.3820355) q[1];
sx q[1];
rz(-1.3146725) q[1];
sx q[1];
rz(1.4256328) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65042881) q[0];
sx q[0];
rz(-0.61507498) q[0];
sx q[0];
rz(-2.0838724) q[0];
rz(-pi) q[1];
rz(-0.47503586) q[2];
sx q[2];
rz(-0.54752195) q[2];
sx q[2];
rz(0.56625596) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5162075) q[1];
sx q[1];
rz(-0.16730669) q[1];
sx q[1];
rz(1.0870666) q[1];
rz(-2.2896122) q[3];
sx q[3];
rz(-1.6036828) q[3];
sx q[3];
rz(-2.6926958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2322959) q[2];
sx q[2];
rz(-2.2984419) q[2];
sx q[2];
rz(0.24093957) q[2];
rz(-0.029189261) q[3];
sx q[3];
rz(-1.8029282) q[3];
sx q[3];
rz(1.9624814) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9772684) q[0];
sx q[0];
rz(-1.9376396) q[0];
sx q[0];
rz(0.48467317) q[0];
rz(-2.4618705) q[1];
sx q[1];
rz(-1.8599963) q[1];
sx q[1];
rz(-1.9814804) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.401448) q[0];
sx q[0];
rz(-2.8370246) q[0];
sx q[0];
rz(0.061621678) q[0];
rz(-2.2641029) q[2];
sx q[2];
rz(-2.1602767) q[2];
sx q[2];
rz(-1.9716919) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1966432) q[1];
sx q[1];
rz(-1.3232051) q[1];
sx q[1];
rz(0.61243122) q[1];
x q[2];
rz(-0.5228225) q[3];
sx q[3];
rz(-1.8588052) q[3];
sx q[3];
rz(3.0436181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.29740563) q[2];
sx q[2];
rz(-0.30218267) q[2];
sx q[2];
rz(-2.6837132) q[2];
rz(-1.1229905) q[3];
sx q[3];
rz(-0.99959683) q[3];
sx q[3];
rz(-0.56639731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8726525) q[0];
sx q[0];
rz(-0.70916969) q[0];
sx q[0];
rz(-3.1100682) q[0];
rz(-0.28745502) q[1];
sx q[1];
rz(-0.87702409) q[1];
sx q[1];
rz(-1.9256176) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7087962) q[0];
sx q[0];
rz(-2.2026718) q[0];
sx q[0];
rz(-2.2461476) q[0];
rz(2.7599665) q[2];
sx q[2];
rz(-0.51593057) q[2];
sx q[2];
rz(-2.6037773) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3423914) q[1];
sx q[1];
rz(-1.8441797) q[1];
sx q[1];
rz(-2.8552516) q[1];
x q[2];
rz(0.6612571) q[3];
sx q[3];
rz(-1.0399858) q[3];
sx q[3];
rz(-0.18839934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1027801) q[2];
sx q[2];
rz(-1.5422042) q[2];
sx q[2];
rz(-2.3056324) q[2];
rz(-1.7838259) q[3];
sx q[3];
rz(-2.416555) q[3];
sx q[3];
rz(2.3939705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.2826071) q[0];
sx q[0];
rz(-1.7361807) q[0];
sx q[0];
rz(-0.080168515) q[0];
rz(-2.0591586) q[1];
sx q[1];
rz(-2.9083462) q[1];
sx q[1];
rz(-1.3287883) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69705582) q[0];
sx q[0];
rz(-1.8955243) q[0];
sx q[0];
rz(2.6022807) q[0];
x q[1];
rz(2.8187739) q[2];
sx q[2];
rz(-1.4070961) q[2];
sx q[2];
rz(-1.1522918) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7419974) q[1];
sx q[1];
rz(-0.4542225) q[1];
sx q[1];
rz(-1.2356367) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8058047) q[3];
sx q[3];
rz(-1.3779252) q[3];
sx q[3];
rz(-1.4953116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7666011) q[2];
sx q[2];
rz(-0.98321715) q[2];
sx q[2];
rz(0.90744606) q[2];
rz(2.4713016) q[3];
sx q[3];
rz(-0.82796103) q[3];
sx q[3];
rz(1.7214187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9012673) q[0];
sx q[0];
rz(-2.2690052) q[0];
sx q[0];
rz(2.7101044) q[0];
rz(-2.0101428) q[1];
sx q[1];
rz(-1.9624458) q[1];
sx q[1];
rz(-2.7010837) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2402396) q[0];
sx q[0];
rz(-2.6640737) q[0];
sx q[0];
rz(1.3130472) q[0];
rz(1.126664) q[2];
sx q[2];
rz(-0.45373771) q[2];
sx q[2];
rz(-0.96378381) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9677093) q[1];
sx q[1];
rz(-1.5313799) q[1];
sx q[1];
rz(2.8110678) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.308485) q[3];
sx q[3];
rz(-0.78189497) q[3];
sx q[3];
rz(1.9946919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.12209192) q[2];
sx q[2];
rz(-2.2152405) q[2];
sx q[2];
rz(2.7276373) q[2];
rz(1.7533938) q[3];
sx q[3];
rz(-1.5940758) q[3];
sx q[3];
rz(-0.12510124) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6269161) q[0];
sx q[0];
rz(-1.4056982) q[0];
sx q[0];
rz(2.2077014) q[0];
rz(-2.4712708) q[1];
sx q[1];
rz(-1.8972081) q[1];
sx q[1];
rz(-1.8062887) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40437296) q[0];
sx q[0];
rz(-1.6707194) q[0];
sx q[0];
rz(0.96802442) q[0];
rz(-pi) q[1];
rz(2.9839758) q[2];
sx q[2];
rz(-2.981957) q[2];
sx q[2];
rz(-0.0052513382) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0508479) q[1];
sx q[1];
rz(-0.11493348) q[1];
sx q[1];
rz(1.7687709) q[1];
rz(-pi) q[2];
rz(1.5794831) q[3];
sx q[3];
rz(-0.93378769) q[3];
sx q[3];
rz(2.6125233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.89207092) q[2];
sx q[2];
rz(-2.8688909) q[2];
sx q[2];
rz(1.2707233) q[2];
rz(-1.3696085) q[3];
sx q[3];
rz(-1.9000051) q[3];
sx q[3];
rz(2.6197267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9524566) q[0];
sx q[0];
rz(-2.554775) q[0];
sx q[0];
rz(2.9879046) q[0];
rz(0.20248374) q[1];
sx q[1];
rz(-1.9053562) q[1];
sx q[1];
rz(-2.1930146) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21009203) q[0];
sx q[0];
rz(-2.4903959) q[0];
sx q[0];
rz(-0.33588846) q[0];
rz(1.1467298) q[2];
sx q[2];
rz(-0.68409195) q[2];
sx q[2];
rz(2.7685194) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.9002418) q[1];
sx q[1];
rz(-2.7903922) q[1];
sx q[1];
rz(-2.3022149) q[1];
rz(-pi) q[2];
rz(2.2205686) q[3];
sx q[3];
rz(-1.463784) q[3];
sx q[3];
rz(1.8821074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1367246) q[2];
sx q[2];
rz(-2.0437045) q[2];
sx q[2];
rz(3.1401805) q[2];
rz(0.062156113) q[3];
sx q[3];
rz(-1.4096189) q[3];
sx q[3];
rz(-2.1803161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.38178) q[0];
sx q[0];
rz(-0.75730046) q[0];
sx q[0];
rz(1.9267474) q[0];
rz(-2.3836721) q[1];
sx q[1];
rz(-0.52752033) q[1];
sx q[1];
rz(-2.5835999) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76282952) q[0];
sx q[0];
rz(-2.7192594) q[0];
sx q[0];
rz(-0.93393737) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9079014) q[2];
sx q[2];
rz(-2.5776641) q[2];
sx q[2];
rz(0.83966161) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.87822014) q[1];
sx q[1];
rz(-2.0505295) q[1];
sx q[1];
rz(-0.31444506) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4766944) q[3];
sx q[3];
rz(-2.8633139) q[3];
sx q[3];
rz(-0.32839963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5672292) q[2];
sx q[2];
rz(-1.1939253) q[2];
sx q[2];
rz(0.26926678) q[2];
rz(1.9783798) q[3];
sx q[3];
rz(-2.5389157) q[3];
sx q[3];
rz(2.9009624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4686541) q[0];
sx q[0];
rz(-1.0373632) q[0];
sx q[0];
rz(-2.0682251) q[0];
rz(-2.0138373) q[1];
sx q[1];
rz(-0.22722166) q[1];
sx q[1];
rz(0.55666298) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23281413) q[0];
sx q[0];
rz(-1.4289843) q[0];
sx q[0];
rz(1.892426) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.9999335) q[2];
sx q[2];
rz(-2.0682671) q[2];
sx q[2];
rz(-1.9299049) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.6491739) q[1];
sx q[1];
rz(-2.2606876) q[1];
sx q[1];
rz(2.0418732) q[1];
x q[2];
rz(-1.6708371) q[3];
sx q[3];
rz(-1.8042068) q[3];
sx q[3];
rz(-1.4518713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2399981) q[2];
sx q[2];
rz(-2.0549213) q[2];
sx q[2];
rz(1.979801) q[2];
rz(-0.53317541) q[3];
sx q[3];
rz(-2.6170001) q[3];
sx q[3];
rz(-2.9840792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6530782) q[0];
sx q[0];
rz(-2.1670659) q[0];
sx q[0];
rz(-2.5467806) q[0];
rz(0.57890233) q[1];
sx q[1];
rz(-1.6659104) q[1];
sx q[1];
rz(2.4822809) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1950451) q[0];
sx q[0];
rz(-2.5016245) q[0];
sx q[0];
rz(-2.9165687) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60229199) q[2];
sx q[2];
rz(-2.4458439) q[2];
sx q[2];
rz(2.2056923) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0486006) q[1];
sx q[1];
rz(-1.326418) q[1];
sx q[1];
rz(2.9952769) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3385184) q[3];
sx q[3];
rz(-2.7702906) q[3];
sx q[3];
rz(-0.40483958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.79002964) q[2];
sx q[2];
rz(-1.943925) q[2];
sx q[2];
rz(0.053038049) q[2];
rz(-1.3430345) q[3];
sx q[3];
rz(-1.2767295) q[3];
sx q[3];
rz(-0.0742577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2863083) q[0];
sx q[0];
rz(-1.7779779) q[0];
sx q[0];
rz(1.5386982) q[0];
rz(-0.098943624) q[1];
sx q[1];
rz(-1.7715441) q[1];
sx q[1];
rz(-3.0861707) q[1];
rz(3.121539) q[2];
sx q[2];
rz(-2.7360624) q[2];
sx q[2];
rz(-0.61780203) q[2];
rz(0.17114279) q[3];
sx q[3];
rz(-2.3928693) q[3];
sx q[3];
rz(0.00015043845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
