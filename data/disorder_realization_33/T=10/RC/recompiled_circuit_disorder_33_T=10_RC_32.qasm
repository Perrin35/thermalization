OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6361976) q[0];
sx q[0];
rz(6.0072748) q[0];
sx q[0];
rz(10.732565) q[0];
rz(-2.0055327) q[1];
sx q[1];
rz(4.0772822) q[1];
sx q[1];
rz(4.7128591) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21068621) q[0];
sx q[0];
rz(-1.9084198) q[0];
sx q[0];
rz(-2.7772285) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2575126) q[2];
sx q[2];
rz(-0.93187098) q[2];
sx q[2];
rz(2.0761348) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1085514) q[1];
sx q[1];
rz(-1.812495) q[1];
sx q[1];
rz(-0.34717314) q[1];
rz(-1.0899815) q[3];
sx q[3];
rz(-0.40502031) q[3];
sx q[3];
rz(2.4570176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.87542614) q[2];
sx q[2];
rz(-2.8484919) q[2];
sx q[2];
rz(2.0092633) q[2];
rz(-1.6752361) q[3];
sx q[3];
rz(-1.8050067) q[3];
sx q[3];
rz(-2.1291389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19673008) q[0];
sx q[0];
rz(-0.20962993) q[0];
sx q[0];
rz(-2.9557513) q[0];
rz(-0.56022412) q[1];
sx q[1];
rz(-1.2954243) q[1];
sx q[1];
rz(0.21683189) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5179948) q[0];
sx q[0];
rz(-1.2835842) q[0];
sx q[0];
rz(2.2395796) q[0];
rz(-pi) q[1];
rz(-1.0160604) q[2];
sx q[2];
rz(-1.9811355) q[2];
sx q[2];
rz(-2.1330657) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5237907) q[1];
sx q[1];
rz(-2.4317867) q[1];
sx q[1];
rz(1.2426504) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84083765) q[3];
sx q[3];
rz(-2.3008122) q[3];
sx q[3];
rz(0.18883146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8314787) q[2];
sx q[2];
rz(-0.82565132) q[2];
sx q[2];
rz(1.8537834) q[2];
rz(0.76256049) q[3];
sx q[3];
rz(-1.1688787) q[3];
sx q[3];
rz(0.30502239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6771616) q[0];
sx q[0];
rz(-0.34496775) q[0];
sx q[0];
rz(2.537354) q[0];
rz(-1.3263946) q[1];
sx q[1];
rz(-1.3605958) q[1];
sx q[1];
rz(2.2089829) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5843825) q[0];
sx q[0];
rz(-0.32214468) q[0];
sx q[0];
rz(-1.4003217) q[0];
rz(-pi) q[1];
rz(2.2432125) q[2];
sx q[2];
rz(-0.2873688) q[2];
sx q[2];
rz(-0.76295602) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4438666) q[1];
sx q[1];
rz(-1.4523456) q[1];
sx q[1];
rz(-2.5812134) q[1];
rz(2.695735) q[3];
sx q[3];
rz(-2.0816396) q[3];
sx q[3];
rz(-1.0872935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9937667) q[2];
sx q[2];
rz(-1.0819165) q[2];
sx q[2];
rz(-1.0926584) q[2];
rz(0.5422194) q[3];
sx q[3];
rz(-2.0565624) q[3];
sx q[3];
rz(-0.96737635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3595235) q[0];
sx q[0];
rz(-3.0451267) q[0];
sx q[0];
rz(2.6413667) q[0];
rz(0.80530986) q[1];
sx q[1];
rz(-1.1601245) q[1];
sx q[1];
rz(1.6436228) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26101199) q[0];
sx q[0];
rz(-1.8646761) q[0];
sx q[0];
rz(1.0512933) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.285032) q[2];
sx q[2];
rz(-0.19837241) q[2];
sx q[2];
rz(-2.6464268) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5768891) q[1];
sx q[1];
rz(-2.7895045) q[1];
sx q[1];
rz(2.4114154) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3718932) q[3];
sx q[3];
rz(-0.61004988) q[3];
sx q[3];
rz(-1.9609914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3952289) q[2];
sx q[2];
rz(-0.56240288) q[2];
sx q[2];
rz(0.70181075) q[2];
rz(-0.83135215) q[3];
sx q[3];
rz(-2.1777007) q[3];
sx q[3];
rz(0.62197661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9005301) q[0];
sx q[0];
rz(-0.59589544) q[0];
sx q[0];
rz(0.81533122) q[0];
rz(1.5218081) q[1];
sx q[1];
rz(-0.83414572) q[1];
sx q[1];
rz(-2.0933847) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5532903) q[0];
sx q[0];
rz(-0.36670812) q[0];
sx q[0];
rz(-2.328863) q[0];
rz(-pi) q[1];
rz(1.5830718) q[2];
sx q[2];
rz(-0.91931146) q[2];
sx q[2];
rz(-0.94142454) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.198846) q[1];
sx q[1];
rz(-1.2578576) q[1];
sx q[1];
rz(-1.320977) q[1];
x q[2];
rz(2.7311677) q[3];
sx q[3];
rz(-2.1459208) q[3];
sx q[3];
rz(2.9068973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6158225) q[2];
sx q[2];
rz(-2.5746391) q[2];
sx q[2];
rz(1.099951) q[2];
rz(-0.82529092) q[3];
sx q[3];
rz(-1.0422948) q[3];
sx q[3];
rz(-0.88551372) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0681756) q[0];
sx q[0];
rz(-0.59403479) q[0];
sx q[0];
rz(0.90240479) q[0];
rz(2.1249318) q[1];
sx q[1];
rz(-2.0817751) q[1];
sx q[1];
rz(-0.12983233) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5141402) q[0];
sx q[0];
rz(-2.1462626) q[0];
sx q[0];
rz(2.0155725) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1075222) q[2];
sx q[2];
rz(-0.64646361) q[2];
sx q[2];
rz(-1.5922286) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8772014) q[1];
sx q[1];
rz(-1.0069205) q[1];
sx q[1];
rz(-2.0045723) q[1];
rz(1.3744266) q[3];
sx q[3];
rz(-1.0305627) q[3];
sx q[3];
rz(-1.2353209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3123902) q[2];
sx q[2];
rz(-0.94909334) q[2];
sx q[2];
rz(-0.20425805) q[2];
rz(1.9355109) q[3];
sx q[3];
rz(-1.5217425) q[3];
sx q[3];
rz(0.23541418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7234574) q[0];
sx q[0];
rz(-1.3292987) q[0];
sx q[0];
rz(1.6947421) q[0];
rz(1.2591259) q[1];
sx q[1];
rz(-2.1513758) q[1];
sx q[1];
rz(2.4553305) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2080363) q[0];
sx q[0];
rz(-2.4718923) q[0];
sx q[0];
rz(-0.74525381) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5574066) q[2];
sx q[2];
rz(-0.91911941) q[2];
sx q[2];
rz(-0.78648957) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2346748) q[1];
sx q[1];
rz(-1.385681) q[1];
sx q[1];
rz(-3.0659552) q[1];
rz(-pi) q[2];
rz(2.5043082) q[3];
sx q[3];
rz(-2.6316959) q[3];
sx q[3];
rz(-0.84264681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.69616047) q[2];
sx q[2];
rz(-1.3779209) q[2];
sx q[2];
rz(-3.1398204) q[2];
rz(-0.56162515) q[3];
sx q[3];
rz(-2.2300945) q[3];
sx q[3];
rz(1.6368438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6034265) q[0];
sx q[0];
rz(-2.4551233) q[0];
sx q[0];
rz(-1.6954533) q[0];
rz(2.360545) q[1];
sx q[1];
rz(-1.3054409) q[1];
sx q[1];
rz(1.6400281) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1743463) q[0];
sx q[0];
rz(-1.5513199) q[0];
sx q[0];
rz(2.0635701) q[0];
rz(-1.5090452) q[2];
sx q[2];
rz(-1.5973063) q[2];
sx q[2];
rz(2.3960631) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9403119) q[1];
sx q[1];
rz(-0.33946013) q[1];
sx q[1];
rz(-2.8708354) q[1];
x q[2];
rz(-0.99174188) q[3];
sx q[3];
rz(-1.0240882) q[3];
sx q[3];
rz(0.21608298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.35187307) q[2];
sx q[2];
rz(-1.3616273) q[2];
sx q[2];
rz(1.3191351) q[2];
rz(1.9296648) q[3];
sx q[3];
rz(-1.8550248) q[3];
sx q[3];
rz(-2.8222728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8050352) q[0];
sx q[0];
rz(-0.55258495) q[0];
sx q[0];
rz(-1.9375027) q[0];
rz(2.7583292) q[1];
sx q[1];
rz(-0.52572322) q[1];
sx q[1];
rz(-0.35167545) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4708913) q[0];
sx q[0];
rz(-1.1180709) q[0];
sx q[0];
rz(0.41505138) q[0];
rz(-pi) q[1];
rz(1.2665389) q[2];
sx q[2];
rz(-0.90196246) q[2];
sx q[2];
rz(2.5077016) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5896776) q[1];
sx q[1];
rz(-0.81112408) q[1];
sx q[1];
rz(-0.07304904) q[1];
rz(-0.10961253) q[3];
sx q[3];
rz(-1.5188367) q[3];
sx q[3];
rz(0.51652858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.3433156) q[2];
sx q[2];
rz(-1.1077935) q[2];
sx q[2];
rz(1.2822255) q[2];
rz(-1.6451689) q[3];
sx q[3];
rz(-1.5346425) q[3];
sx q[3];
rz(-1.055868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6431817) q[0];
sx q[0];
rz(-1.2675985) q[0];
sx q[0];
rz(2.9472651) q[0];
rz(-1.0378029) q[1];
sx q[1];
rz(-0.56832814) q[1];
sx q[1];
rz(1.0338354) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72458306) q[0];
sx q[0];
rz(-1.7010744) q[0];
sx q[0];
rz(-2.2307322) q[0];
rz(0.86976544) q[2];
sx q[2];
rz(-1.1681721) q[2];
sx q[2];
rz(-2.082777) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0018113) q[1];
sx q[1];
rz(-2.1404631) q[1];
sx q[1];
rz(0.12057481) q[1];
rz(2.8909573) q[3];
sx q[3];
rz(-0.72892979) q[3];
sx q[3];
rz(0.38754101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0620492) q[2];
sx q[2];
rz(-2.1958308) q[2];
sx q[2];
rz(-2.5058084) q[2];
rz(-2.87129) q[3];
sx q[3];
rz(-0.79939866) q[3];
sx q[3];
rz(-1.5283782) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4476267) q[0];
sx q[0];
rz(-1.8287369) q[0];
sx q[0];
rz(1.0736314) q[0];
rz(1.4355961) q[1];
sx q[1];
rz(-1.5789079) q[1];
sx q[1];
rz(0.78067738) q[1];
rz(1.5244665) q[2];
sx q[2];
rz(-0.6033069) q[2];
sx q[2];
rz(2.6848007) q[2];
rz(-0.070449645) q[3];
sx q[3];
rz(-1.2415213) q[3];
sx q[3];
rz(0.51125676) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];