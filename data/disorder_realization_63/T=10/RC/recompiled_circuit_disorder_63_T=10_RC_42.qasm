OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1749984) q[0];
sx q[0];
rz(-0.35342616) q[0];
sx q[0];
rz(1.0647635) q[0];
rz(-2.3454173) q[1];
sx q[1];
rz(-1.2086955) q[1];
sx q[1];
rz(2.605521) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.160643) q[0];
sx q[0];
rz(-1.0796709) q[0];
sx q[0];
rz(2.8777697) q[0];
rz(1.6494765) q[2];
sx q[2];
rz(-2.4298551) q[2];
sx q[2];
rz(-0.56456883) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.55878996) q[1];
sx q[1];
rz(-1.3830118) q[1];
sx q[1];
rz(2.5488528) q[1];
x q[2];
rz(0.53341289) q[3];
sx q[3];
rz(-2.3331113) q[3];
sx q[3];
rz(-1.7318788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7321695) q[2];
sx q[2];
rz(-0.44115856) q[2];
sx q[2];
rz(-2.1146963) q[2];
rz(-2.8895767) q[3];
sx q[3];
rz(-1.9988632) q[3];
sx q[3];
rz(-0.86565971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.058763) q[0];
sx q[0];
rz(-2.7936462) q[0];
sx q[0];
rz(-1.7513562) q[0];
rz(-0.83579666) q[1];
sx q[1];
rz(-0.73671571) q[1];
sx q[1];
rz(-2.4332411) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10193292) q[0];
sx q[0];
rz(-1.9358557) q[0];
sx q[0];
rz(1.2714766) q[0];
rz(-pi) q[1];
rz(1.3588261) q[2];
sx q[2];
rz(-2.5296629) q[2];
sx q[2];
rz(1.7906534) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2409776) q[1];
sx q[1];
rz(-0.88417378) q[1];
sx q[1];
rz(0.91614206) q[1];
rz(-pi) q[2];
rz(2.0833011) q[3];
sx q[3];
rz(-1.1304434) q[3];
sx q[3];
rz(0.022170598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11671242) q[2];
sx q[2];
rz(-0.79170266) q[2];
sx q[2];
rz(-0.29176816) q[2];
rz(-0.10270384) q[3];
sx q[3];
rz(-1.4029968) q[3];
sx q[3];
rz(-1.5244938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8125238) q[0];
sx q[0];
rz(-0.085443184) q[0];
sx q[0];
rz(-2.8318751) q[0];
rz(-1.6614871) q[1];
sx q[1];
rz(-1.8146351) q[1];
sx q[1];
rz(0.57166878) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088061995) q[0];
sx q[0];
rz(-0.54170875) q[0];
sx q[0];
rz(-0.17266973) q[0];
rz(-2.2505635) q[2];
sx q[2];
rz(-2.6030412) q[2];
sx q[2];
rz(2.0856759) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3571615) q[1];
sx q[1];
rz(-0.845134) q[1];
sx q[1];
rz(-1.3933338) q[1];
rz(-pi) q[2];
x q[2];
rz(0.27922697) q[3];
sx q[3];
rz(-2.673827) q[3];
sx q[3];
rz(-1.8750909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9812575) q[2];
sx q[2];
rz(-2.2312639) q[2];
sx q[2];
rz(0.70880115) q[2];
rz(-2.839084) q[3];
sx q[3];
rz(-1.442028) q[3];
sx q[3];
rz(-1.9992874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70808327) q[0];
sx q[0];
rz(-2.0286562) q[0];
sx q[0];
rz(2.3420912) q[0];
rz(-3.0917621) q[1];
sx q[1];
rz(-2.2466876) q[1];
sx q[1];
rz(2.9611011) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9468964) q[0];
sx q[0];
rz(-1.4799656) q[0];
sx q[0];
rz(-0.38740654) q[0];
rz(-pi) q[1];
rz(1.7454342) q[2];
sx q[2];
rz(-1.2067814) q[2];
sx q[2];
rz(0.82496914) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.71112878) q[1];
sx q[1];
rz(-2.8006878) q[1];
sx q[1];
rz(0.35699637) q[1];
x q[2];
rz(1.5933851) q[3];
sx q[3];
rz(-0.49593192) q[3];
sx q[3];
rz(0.13726928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0522456) q[2];
sx q[2];
rz(-1.1648488) q[2];
sx q[2];
rz(1.583741) q[2];
rz(-1.3752939) q[3];
sx q[3];
rz(-1.3170653) q[3];
sx q[3];
rz(2.2089675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.76064008) q[0];
sx q[0];
rz(-0.33712688) q[0];
sx q[0];
rz(-0.98651648) q[0];
rz(-1.1622693) q[1];
sx q[1];
rz(-1.9245851) q[1];
sx q[1];
rz(-2.8900237) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0753206) q[0];
sx q[0];
rz(-2.2499654) q[0];
sx q[0];
rz(-2.9156296) q[0];
x q[1];
rz(-1.4839843) q[2];
sx q[2];
rz(-1.112339) q[2];
sx q[2];
rz(0.77755962) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1269762) q[1];
sx q[1];
rz(-1.7661621) q[1];
sx q[1];
rz(-1.5429392) q[1];
rz(-pi) q[2];
rz(2.1523347) q[3];
sx q[3];
rz(-0.79205081) q[3];
sx q[3];
rz(-2.0056412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.183737) q[2];
sx q[2];
rz(-0.541406) q[2];
sx q[2];
rz(-1.0305369) q[2];
rz(-1.41097) q[3];
sx q[3];
rz(-2.1822699) q[3];
sx q[3];
rz(-2.5103536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75491607) q[0];
sx q[0];
rz(-0.47575352) q[0];
sx q[0];
rz(-0.85012287) q[0];
rz(-1.9592346) q[1];
sx q[1];
rz(-1.9071002) q[1];
sx q[1];
rz(-2.7485671) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.128199) q[0];
sx q[0];
rz(-1.734874) q[0];
sx q[0];
rz(-1.8830995) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8259465) q[2];
sx q[2];
rz(-2.3806551) q[2];
sx q[2];
rz(1.4756502) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.46950754) q[1];
sx q[1];
rz(-2.2219594) q[1];
sx q[1];
rz(0.0068411946) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29720184) q[3];
sx q[3];
rz(-1.3085877) q[3];
sx q[3];
rz(-2.8924243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7810016) q[2];
sx q[2];
rz(-2.0157308) q[2];
sx q[2];
rz(-2.1684516) q[2];
rz(-0.94758236) q[3];
sx q[3];
rz(-1.1004473) q[3];
sx q[3];
rz(-1.9472286) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6958375) q[0];
sx q[0];
rz(-2.7212302) q[0];
sx q[0];
rz(1.6181035) q[0];
rz(2.7667926) q[1];
sx q[1];
rz(-1.8811036) q[1];
sx q[1];
rz(-3.135625) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8864266) q[0];
sx q[0];
rz(-0.92029858) q[0];
sx q[0];
rz(-0.35264539) q[0];
rz(-pi) q[1];
rz(2.3890424) q[2];
sx q[2];
rz(-1.9117022) q[2];
sx q[2];
rz(-2.7698851) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55885) q[1];
sx q[1];
rz(-0.55816459) q[1];
sx q[1];
rz(-1.6911669) q[1];
rz(-pi) q[2];
rz(1.5238477) q[3];
sx q[3];
rz(-1.1991683) q[3];
sx q[3];
rz(2.9058128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.57198793) q[2];
sx q[2];
rz(-2.5964952) q[2];
sx q[2];
rz(0.19006426) q[2];
rz(-0.41641411) q[3];
sx q[3];
rz(-2.1834686) q[3];
sx q[3];
rz(0.73474187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1104601) q[0];
sx q[0];
rz(-1.0460331) q[0];
sx q[0];
rz(-2.6053612) q[0];
rz(0.93332943) q[1];
sx q[1];
rz(-1.5678762) q[1];
sx q[1];
rz(-0.94820625) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95420102) q[0];
sx q[0];
rz(-1.3905977) q[0];
sx q[0];
rz(-2.3733634) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2657329) q[2];
sx q[2];
rz(-0.63820733) q[2];
sx q[2];
rz(-1.1832331) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1369689) q[1];
sx q[1];
rz(-1.3490632) q[1];
sx q[1];
rz(1.5066513) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6471345) q[3];
sx q[3];
rz(-1.2148972) q[3];
sx q[3];
rz(2.5446041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.640921) q[2];
sx q[2];
rz(-1.6225953) q[2];
sx q[2];
rz(2.9013157) q[2];
rz(-0.37929532) q[3];
sx q[3];
rz(-1.0255739) q[3];
sx q[3];
rz(-1.3482288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76294476) q[0];
sx q[0];
rz(-2.8084016) q[0];
sx q[0];
rz(-0.25094029) q[0];
rz(0.25302408) q[1];
sx q[1];
rz(-1.7605942) q[1];
sx q[1];
rz(-0.35266638) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34516016) q[0];
sx q[0];
rz(-2.9734018) q[0];
sx q[0];
rz(2.3229984) q[0];
x q[1];
rz(2.7025931) q[2];
sx q[2];
rz(-2.8121431) q[2];
sx q[2];
rz(1.1169369) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0295804) q[1];
sx q[1];
rz(-2.5332846) q[1];
sx q[1];
rz(0.40257247) q[1];
rz(-pi) q[2];
x q[2];
rz(2.114186) q[3];
sx q[3];
rz(-0.81334844) q[3];
sx q[3];
rz(3.120595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3744366) q[2];
sx q[2];
rz(-2.0220951) q[2];
sx q[2];
rz(-2.4582668) q[2];
rz(-1.654489) q[3];
sx q[3];
rz(-2.7023102) q[3];
sx q[3];
rz(-1.1504415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3764573) q[0];
sx q[0];
rz(-0.52381223) q[0];
sx q[0];
rz(1.3002243) q[0];
rz(-0.70872778) q[1];
sx q[1];
rz(-2.6649902) q[1];
sx q[1];
rz(2.2050819) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4721601) q[0];
sx q[0];
rz(-0.59211187) q[0];
sx q[0];
rz(0.44606146) q[0];
x q[1];
rz(0.18239393) q[2];
sx q[2];
rz(-1.2904022) q[2];
sx q[2];
rz(0.80255752) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9156956) q[1];
sx q[1];
rz(-0.86663336) q[1];
sx q[1];
rz(-0.26030635) q[1];
rz(-pi) q[2];
rz(0.98827045) q[3];
sx q[3];
rz(-1.8347782) q[3];
sx q[3];
rz(1.1519983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5986754) q[2];
sx q[2];
rz(-0.3391372) q[2];
sx q[2];
rz(-3.1266406) q[2];
rz(0.22710083) q[3];
sx q[3];
rz(-0.98277074) q[3];
sx q[3];
rz(-1.2623513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1508355) q[0];
sx q[0];
rz(-0.80934722) q[0];
sx q[0];
rz(1.9833175) q[0];
rz(1.5630209) q[1];
sx q[1];
rz(-2.38588) q[1];
sx q[1];
rz(-0.26185782) q[1];
rz(2.0832534) q[2];
sx q[2];
rz(-2.8556311) q[2];
sx q[2];
rz(-0.43559504) q[2];
rz(2.7097377) q[3];
sx q[3];
rz(-2.5681744) q[3];
sx q[3];
rz(-1.4011866) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
