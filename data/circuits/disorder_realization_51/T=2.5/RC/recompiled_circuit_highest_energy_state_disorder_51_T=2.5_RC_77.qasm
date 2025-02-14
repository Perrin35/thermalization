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
rz(2.3866374) q[0];
sx q[0];
rz(-0.75140262) q[0];
sx q[0];
rz(-1.0487392) q[0];
rz(-2.7222848) q[1];
sx q[1];
rz(4.5276548) q[1];
sx q[1];
rz(6.164896) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48003475) q[0];
sx q[0];
rz(-1.841651) q[0];
sx q[0];
rz(0.99390219) q[0];
rz(-pi) q[1];
rz(-0.050968214) q[2];
sx q[2];
rz(-1.5545003) q[2];
sx q[2];
rz(2.7560134) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1227545) q[1];
sx q[1];
rz(-1.5411106) q[1];
sx q[1];
rz(-1.5537986) q[1];
rz(1.6476326) q[3];
sx q[3];
rz(-1.8733403) q[3];
sx q[3];
rz(-1.8267814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.92363858) q[2];
sx q[2];
rz(-0.0035692735) q[2];
sx q[2];
rz(-2.9812319) q[2];
rz(1.977836) q[3];
sx q[3];
rz(-2.0573503) q[3];
sx q[3];
rz(0.81419301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1639975) q[0];
sx q[0];
rz(-1.4871335) q[0];
sx q[0];
rz(-2.1556222) q[0];
rz(-1.5522955) q[1];
sx q[1];
rz(-2.8542216) q[1];
sx q[1];
rz(1.5812965) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2575657) q[0];
sx q[0];
rz(-2.0456893) q[0];
sx q[0];
rz(-2.7090461) q[0];
rz(3.1104187) q[2];
sx q[2];
rz(-0.59470648) q[2];
sx q[2];
rz(0.020356914) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6735005) q[1];
sx q[1];
rz(-1.4080268) q[1];
sx q[1];
rz(-2.7767608) q[1];
x q[2];
rz(-0.12388568) q[3];
sx q[3];
rz(-2.200211) q[3];
sx q[3];
rz(0.38485011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5220149) q[2];
sx q[2];
rz(-1.3238944) q[2];
sx q[2];
rz(-2.1066693) q[2];
rz(-0.50152913) q[3];
sx q[3];
rz(-0.087567121) q[3];
sx q[3];
rz(0.97626221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72306776) q[0];
sx q[0];
rz(-1.9880966) q[0];
sx q[0];
rz(-1.204741) q[0];
rz(-1.0459666) q[1];
sx q[1];
rz(-0.090066411) q[1];
sx q[1];
rz(-0.14030309) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3922962) q[0];
sx q[0];
rz(-0.97784144) q[0];
sx q[0];
rz(2.8234574) q[0];
rz(-0.42201747) q[2];
sx q[2];
rz(-0.8895424) q[2];
sx q[2];
rz(1.8651419) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2724803) q[1];
sx q[1];
rz(-1.7971874) q[1];
sx q[1];
rz(2.5416944) q[1];
rz(-0.14032538) q[3];
sx q[3];
rz(-1.9092536) q[3];
sx q[3];
rz(0.1986157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3673765) q[2];
sx q[2];
rz(-2.1420631) q[2];
sx q[2];
rz(2.9191169) q[2];
rz(3.0761278) q[3];
sx q[3];
rz(-1.2832578) q[3];
sx q[3];
rz(0.84622598) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9533933) q[0];
sx q[0];
rz(-0.024217483) q[0];
sx q[0];
rz(-0.54704332) q[0];
rz(0.25829265) q[1];
sx q[1];
rz(-0.02198418) q[1];
sx q[1];
rz(-2.8000854) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19778457) q[0];
sx q[0];
rz(-0.0066692624) q[0];
sx q[0];
rz(-1.7775373) q[0];
rz(-pi) q[1];
rz(1.1522305) q[2];
sx q[2];
rz(-1.2593049) q[2];
sx q[2];
rz(-2.3626987) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8402183) q[1];
sx q[1];
rz(-2.3010588) q[1];
sx q[1];
rz(-0.45481843) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3157201) q[3];
sx q[3];
rz(-2.4950728) q[3];
sx q[3];
rz(2.7310128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7121938) q[2];
sx q[2];
rz(-1.8455576) q[2];
sx q[2];
rz(-2.1923501) q[2];
rz(2.3927169) q[3];
sx q[3];
rz(-1.8604934) q[3];
sx q[3];
rz(-3.0794028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5267938) q[0];
sx q[0];
rz(-3.1068046) q[0];
sx q[0];
rz(1.5507966) q[0];
rz(1.3560449) q[1];
sx q[1];
rz(-0.0043914774) q[1];
sx q[1];
rz(0.063025085) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18308355) q[0];
sx q[0];
rz(-1.5048072) q[0];
sx q[0];
rz(-1.5404758) q[0];
rz(-pi) q[1];
x q[1];
rz(0.56383987) q[2];
sx q[2];
rz(-0.83602521) q[2];
sx q[2];
rz(-0.47845248) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1511475) q[1];
sx q[1];
rz(-1.778981) q[1];
sx q[1];
rz(3.1143536) q[1];
rz(-pi) q[2];
rz(0.47565094) q[3];
sx q[3];
rz(-0.84053549) q[3];
sx q[3];
rz(0.85640872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4847792) q[2];
sx q[2];
rz(-1.8317089) q[2];
sx q[2];
rz(-2.4978034) q[2];
rz(-2.2667609) q[3];
sx q[3];
rz(-0.30431408) q[3];
sx q[3];
rz(2.3410102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0463882) q[0];
sx q[0];
rz(-0.055618532) q[0];
sx q[0];
rz(-0.355542) q[0];
rz(2.9429759) q[1];
sx q[1];
rz(-0.0067409975) q[1];
sx q[1];
rz(-2.9933062) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4041408) q[0];
sx q[0];
rz(-0.1830398) q[0];
sx q[0];
rz(-3.1279081) q[0];
rz(-pi) q[1];
rz(-2.9695187) q[2];
sx q[2];
rz(-2.7248757) q[2];
sx q[2];
rz(-0.91815776) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4365789) q[1];
sx q[1];
rz(-1.3376029) q[1];
sx q[1];
rz(1.7089273) q[1];
x q[2];
rz(1.8174174) q[3];
sx q[3];
rz(-2.518836) q[3];
sx q[3];
rz(-2.0758219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4699012) q[2];
sx q[2];
rz(-2.9000059) q[2];
sx q[2];
rz(-3.0174603) q[2];
rz(0.5747253) q[3];
sx q[3];
rz(-0.14437965) q[3];
sx q[3];
rz(-0.1709443) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8827051) q[0];
sx q[0];
rz(-3.0165065) q[0];
sx q[0];
rz(0.74147725) q[0];
rz(2.8575836) q[1];
sx q[1];
rz(-3.1378855) q[1];
sx q[1];
rz(0.31518087) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7595939) q[0];
sx q[0];
rz(-1.505251) q[0];
sx q[0];
rz(-3.1146953) q[0];
rz(-pi) q[1];
rz(0.20756794) q[2];
sx q[2];
rz(-1.1313442) q[2];
sx q[2];
rz(-2.07711) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.3005581) q[1];
sx q[1];
rz(-1.729064) q[1];
sx q[1];
rz(-0.59945089) q[1];
rz(-pi) q[2];
rz(-0.02353286) q[3];
sx q[3];
rz(-1.7943903) q[3];
sx q[3];
rz(1.6361039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2726941) q[2];
sx q[2];
rz(-2.0468476) q[2];
sx q[2];
rz(2.4197253) q[2];
rz(-2.7591211) q[3];
sx q[3];
rz(-1.1451274) q[3];
sx q[3];
rz(-1.0684048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5751936) q[0];
sx q[0];
rz(-0.02481758) q[0];
sx q[0];
rz(1.5665293) q[0];
rz(-2.9381835) q[1];
sx q[1];
rz(-1.2982439) q[1];
sx q[1];
rz(0.64483109) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1328076) q[0];
sx q[0];
rz(-1.5624579) q[0];
sx q[0];
rz(-1.7783283) q[0];
rz(-pi) q[1];
rz(-1.5956085) q[2];
sx q[2];
rz(-2.1729219) q[2];
sx q[2];
rz(0.40222049) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3159807) q[1];
sx q[1];
rz(-0.54345268) q[1];
sx q[1];
rz(1.4280591) q[1];
rz(-pi) q[2];
rz(1.1172148) q[3];
sx q[3];
rz(-2.0189813) q[3];
sx q[3];
rz(-0.77806015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5975534) q[2];
sx q[2];
rz(-0.35352239) q[2];
sx q[2];
rz(-2.7618347) q[2];
rz(1.0455658) q[3];
sx q[3];
rz(-1.2328204) q[3];
sx q[3];
rz(1.9426965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3658635) q[0];
sx q[0];
rz(-3.1079223) q[0];
sx q[0];
rz(-1.7849543) q[0];
rz(2.7011073) q[1];
sx q[1];
rz(-2.0511274) q[1];
sx q[1];
rz(0.7007362) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8829699) q[0];
sx q[0];
rz(-2.077266) q[0];
sx q[0];
rz(2.3372786) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9435025) q[2];
sx q[2];
rz(-2.2513394) q[2];
sx q[2];
rz(2.1777505) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8776013) q[1];
sx q[1];
rz(-0.7670247) q[1];
sx q[1];
rz(2.8061538) q[1];
rz(0.0088720284) q[3];
sx q[3];
rz(-2.7089467) q[3];
sx q[3];
rz(0.052730058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7880154) q[2];
sx q[2];
rz(-2.7699296) q[2];
sx q[2];
rz(-1.8549982) q[2];
rz(2.6364117) q[3];
sx q[3];
rz(-0.44613871) q[3];
sx q[3];
rz(1.9193468) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6219567) q[0];
sx q[0];
rz(-3.0918047) q[0];
sx q[0];
rz(1.5995481) q[0];
rz(-2.385335) q[1];
sx q[1];
rz(-0.007096346) q[1];
sx q[1];
rz(-0.33682987) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3339506) q[0];
sx q[0];
rz(-1.5686638) q[0];
sx q[0];
rz(0.3176519) q[0];
x q[1];
rz(2.1977399) q[2];
sx q[2];
rz(-0.34239951) q[2];
sx q[2];
rz(-1.8998208) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5440798) q[1];
sx q[1];
rz(-1.656053) q[1];
sx q[1];
rz(1.6823177) q[1];
rz(-2.7261717) q[3];
sx q[3];
rz(-1.0491228) q[3];
sx q[3];
rz(0.085062438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.51356641) q[2];
sx q[2];
rz(-2.1921373) q[2];
sx q[2];
rz(-0.27964082) q[2];
rz(2.5102992) q[3];
sx q[3];
rz(-0.93835962) q[3];
sx q[3];
rz(2.5173371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30017988) q[0];
sx q[0];
rz(-1.5501839) q[0];
sx q[0];
rz(-1.3612904) q[0];
rz(-2.367876) q[1];
sx q[1];
rz(-0.63540375) q[1];
sx q[1];
rz(0.2159963) q[1];
rz(2.5071267) q[2];
sx q[2];
rz(-0.96426156) q[2];
sx q[2];
rz(-2.5360863) q[2];
rz(-0.37626304) q[3];
sx q[3];
rz(-1.5496764) q[3];
sx q[3];
rz(-1.5839034) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
