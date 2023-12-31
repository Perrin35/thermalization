OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.66184008) q[0];
sx q[0];
rz(2.2979484) q[0];
sx q[0];
rz(9.2568682) q[0];
rz(1.1711988) q[1];
sx q[1];
rz(3.436915) q[1];
sx q[1];
rz(9.480939) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4958772) q[0];
sx q[0];
rz(-2.0464532) q[0];
sx q[0];
rz(-0.23947421) q[0];
rz(2.216823) q[2];
sx q[2];
rz(-1.3999108) q[2];
sx q[2];
rz(-0.31121635) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1311156) q[1];
sx q[1];
rz(-1.8382204) q[1];
sx q[1];
rz(-2.1285776) q[1];
x q[2];
rz(-0.25986259) q[3];
sx q[3];
rz(-1.5279603) q[3];
sx q[3];
rz(2.6916137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7636259) q[2];
sx q[2];
rz(-2.8597735) q[2];
sx q[2];
rz(-2.7089233) q[2];
rz(-1.9487322) q[3];
sx q[3];
rz(-1.2377219) q[3];
sx q[3];
rz(-2.7584934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52779657) q[0];
sx q[0];
rz(-0.48848099) q[0];
sx q[0];
rz(-1.8288076) q[0];
rz(2.9361172) q[1];
sx q[1];
rz(-2.165129) q[1];
sx q[1];
rz(-1.1516494) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0002999) q[0];
sx q[0];
rz(-2.2465758) q[0];
sx q[0];
rz(-2.7594901) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0863016) q[2];
sx q[2];
rz(-2.0995579) q[2];
sx q[2];
rz(0.27258401) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.18560219) q[1];
sx q[1];
rz(-2.5163109) q[1];
sx q[1];
rz(2.37466) q[1];
rz(-pi) q[2];
x q[2];
rz(2.504911) q[3];
sx q[3];
rz(-2.5425306) q[3];
sx q[3];
rz(0.77081313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1318704) q[2];
sx q[2];
rz(-1.4902318) q[2];
sx q[2];
rz(0.22182626) q[2];
rz(-0.37718537) q[3];
sx q[3];
rz(-2.714034) q[3];
sx q[3];
rz(-2.3691573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31056988) q[0];
sx q[0];
rz(-0.092386827) q[0];
sx q[0];
rz(-3.1047399) q[0];
rz(-0.82551461) q[1];
sx q[1];
rz(-1.8258391) q[1];
sx q[1];
rz(-0.056578606) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7650334) q[0];
sx q[0];
rz(-2.3925836) q[0];
sx q[0];
rz(2.2545933) q[0];
x q[1];
rz(2.668004) q[2];
sx q[2];
rz(-2.8550672) q[2];
sx q[2];
rz(0.63602704) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6371582) q[1];
sx q[1];
rz(-1.8706026) q[1];
sx q[1];
rz(-1.7791041) q[1];
rz(-pi) q[2];
x q[2];
rz(2.644396) q[3];
sx q[3];
rz(-2.433784) q[3];
sx q[3];
rz(1.7266112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.27292192) q[2];
sx q[2];
rz(-1.6438831) q[2];
sx q[2];
rz(-0.92612129) q[2];
rz(2.5849294) q[3];
sx q[3];
rz(-0.29354468) q[3];
sx q[3];
rz(-2.0986957) q[3];
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
rz(-pi) q[0];
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
rz(0.91519231) q[0];
sx q[0];
rz(-2.4214348) q[0];
sx q[0];
rz(0.74209374) q[0];
rz(1.1391976) q[1];
sx q[1];
rz(-2.6622055) q[1];
sx q[1];
rz(-0.46359584) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0915506) q[0];
sx q[0];
rz(-0.76489641) q[0];
sx q[0];
rz(2.0043623) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6323339) q[2];
sx q[2];
rz(-1.9029641) q[2];
sx q[2];
rz(-2.5663944) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.61958085) q[1];
sx q[1];
rz(-1.7420235) q[1];
sx q[1];
rz(-0.79230688) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2673244) q[3];
sx q[3];
rz(-2.4618751) q[3];
sx q[3];
rz(-0.10248871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3670369) q[2];
sx q[2];
rz(-2.98428) q[2];
sx q[2];
rz(3.0920933) q[2];
rz(3.0130623) q[3];
sx q[3];
rz(-1.5481719) q[3];
sx q[3];
rz(-3.0226504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0054935) q[0];
sx q[0];
rz(-2.7042784) q[0];
sx q[0];
rz(-0.29770011) q[0];
rz(-2.659335) q[1];
sx q[1];
rz(-2.3869956) q[1];
sx q[1];
rz(-0.94435) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1221065) q[0];
sx q[0];
rz(-1.6158316) q[0];
sx q[0];
rz(-1.7800063) q[0];
rz(-pi) q[1];
rz(0.77563939) q[2];
sx q[2];
rz(-1.8619985) q[2];
sx q[2];
rz(1.2736543) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.797232) q[1];
sx q[1];
rz(-1.5096944) q[1];
sx q[1];
rz(-0.0021211591) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7086908) q[3];
sx q[3];
rz(-1.3372984) q[3];
sx q[3];
rz(1.7683065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8828316) q[2];
sx q[2];
rz(-2.0488887) q[2];
sx q[2];
rz(0.10822254) q[2];
rz(0.0023068874) q[3];
sx q[3];
rz(-1.6121515) q[3];
sx q[3];
rz(0.32430696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7047983) q[0];
sx q[0];
rz(-0.44610766) q[0];
sx q[0];
rz(0.58445245) q[0];
rz(2.2553518) q[1];
sx q[1];
rz(-0.61683547) q[1];
sx q[1];
rz(0.054919682) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14558218) q[0];
sx q[0];
rz(-1.8390363) q[0];
sx q[0];
rz(0.085573816) q[0];
x q[1];
rz(1.6266277) q[2];
sx q[2];
rz(-2.8179114) q[2];
sx q[2];
rz(-1.2602381) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0262895) q[1];
sx q[1];
rz(-0.68740986) q[1];
sx q[1];
rz(-2.5251212) q[1];
x q[2];
rz(-0.98723282) q[3];
sx q[3];
rz(-1.1757441) q[3];
sx q[3];
rz(2.6810255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.75446689) q[2];
sx q[2];
rz(-0.12067623) q[2];
sx q[2];
rz(2.1248655) q[2];
rz(2.5975442) q[3];
sx q[3];
rz(-0.35990158) q[3];
sx q[3];
rz(-1.8090766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.6761557) q[0];
sx q[0];
rz(-2.1570719) q[0];
sx q[0];
rz(2.8570535) q[0];
rz(-0.94447213) q[1];
sx q[1];
rz(-1.1962793) q[1];
sx q[1];
rz(0.91032666) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1658265) q[0];
sx q[0];
rz(-1.5162139) q[0];
sx q[0];
rz(1.5058917) q[0];
rz(0.89332135) q[2];
sx q[2];
rz(-2.6258694) q[2];
sx q[2];
rz(-2.233778) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4714204) q[1];
sx q[1];
rz(-0.52392611) q[1];
sx q[1];
rz(-2.7736204) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5038239) q[3];
sx q[3];
rz(-2.31156) q[3];
sx q[3];
rz(-0.72731599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7515144) q[2];
sx q[2];
rz(-3.0780767) q[2];
sx q[2];
rz(-2.2195623) q[2];
rz(-2.5743124) q[3];
sx q[3];
rz(-1.4138979) q[3];
sx q[3];
rz(-1.0197619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89408016) q[0];
sx q[0];
rz(-0.67665726) q[0];
sx q[0];
rz(-0.12938736) q[0];
rz(2.5091876) q[1];
sx q[1];
rz(-1.0267195) q[1];
sx q[1];
rz(0.30050373) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6209517) q[0];
sx q[0];
rz(-1.2004939) q[0];
sx q[0];
rz(0.83604367) q[0];
rz(0.95327611) q[2];
sx q[2];
rz(-2.2308908) q[2];
sx q[2];
rz(2.2613139) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0771675) q[1];
sx q[1];
rz(-0.38591138) q[1];
sx q[1];
rz(-2.6618631) q[1];
rz(0.81340202) q[3];
sx q[3];
rz(-1.3374995) q[3];
sx q[3];
rz(-2.5113311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5552716) q[2];
sx q[2];
rz(-2.1466612) q[2];
sx q[2];
rz(-2.3596181) q[2];
rz(2.590495) q[3];
sx q[3];
rz(-1.7588153) q[3];
sx q[3];
rz(0.15792318) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5683811) q[0];
sx q[0];
rz(-1.1567572) q[0];
sx q[0];
rz(3.0138299) q[0];
rz(2.5993775) q[1];
sx q[1];
rz(-0.95710373) q[1];
sx q[1];
rz(0.75884563) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6700867) q[0];
sx q[0];
rz(-0.42450464) q[0];
sx q[0];
rz(-1.5295117) q[0];
rz(-pi) q[1];
rz(-1.7636289) q[2];
sx q[2];
rz(-0.97059965) q[2];
sx q[2];
rz(-0.45229518) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.54770494) q[1];
sx q[1];
rz(-0.87915671) q[1];
sx q[1];
rz(1.7734852) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2364053) q[3];
sx q[3];
rz(-1.5593411) q[3];
sx q[3];
rz(0.36597914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0163429) q[2];
sx q[2];
rz(-1.7813851) q[2];
sx q[2];
rz(2.6515567) q[2];
rz(-1.4222493) q[3];
sx q[3];
rz(-1.2079206) q[3];
sx q[3];
rz(1.0629883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7816417) q[0];
sx q[0];
rz(-2.521823) q[0];
sx q[0];
rz(-3.066257) q[0];
rz(-2.244859) q[1];
sx q[1];
rz(-1.9672085) q[1];
sx q[1];
rz(0.60992253) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1760575) q[0];
sx q[0];
rz(-0.81747222) q[0];
sx q[0];
rz(0.39359351) q[0];
x q[1];
rz(0.62060771) q[2];
sx q[2];
rz(-2.6891516) q[2];
sx q[2];
rz(-2.0394182) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8049106) q[1];
sx q[1];
rz(-1.8098127) q[1];
sx q[1];
rz(2.3906624) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3992357) q[3];
sx q[3];
rz(-1.9543813) q[3];
sx q[3];
rz(2.50768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.909409) q[2];
sx q[2];
rz(-0.82513088) q[2];
sx q[2];
rz(-0.71371901) q[2];
rz(0.37832007) q[3];
sx q[3];
rz(-0.49351966) q[3];
sx q[3];
rz(2.2617214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.338035) q[0];
sx q[0];
rz(-1.9914347) q[0];
sx q[0];
rz(1.5557355) q[0];
rz(-2.4907885) q[1];
sx q[1];
rz(-1.6497859) q[1];
sx q[1];
rz(-0.12129687) q[1];
rz(-1.7667608) q[2];
sx q[2];
rz(-1.6008196) q[2];
sx q[2];
rz(0.65884789) q[2];
rz(0.25272947) q[3];
sx q[3];
rz(-1.1128294) q[3];
sx q[3];
rz(2.2313234) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
