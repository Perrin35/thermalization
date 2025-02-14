OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2090981) q[0];
sx q[0];
rz(0.068109186) q[0];
sx q[0];
rz(13.343233) q[0];
rz(-4.4486899) q[1];
sx q[1];
rz(-2.1721462) q[1];
sx q[1];
rz(8.1028508) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37047077) q[0];
sx q[0];
rz(-1.7092472) q[0];
sx q[0];
rz(-0.8531424) q[0];
x q[1];
rz(2.5547682) q[2];
sx q[2];
rz(-0.57673645) q[2];
sx q[2];
rz(-1.5453892) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2474587) q[1];
sx q[1];
rz(-1.1620635) q[1];
sx q[1];
rz(0.031886727) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8613937) q[3];
sx q[3];
rz(-1.3826256) q[3];
sx q[3];
rz(3.1332603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.063735828) q[2];
sx q[2];
rz(-1.7813762) q[2];
sx q[2];
rz(-0.35749164) q[2];
rz(0.56719559) q[3];
sx q[3];
rz(-1.7287858) q[3];
sx q[3];
rz(0.43958694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95328632) q[0];
sx q[0];
rz(-0.28306857) q[0];
sx q[0];
rz(-1.333492) q[0];
rz(2.7045344) q[1];
sx q[1];
rz(-0.60085618) q[1];
sx q[1];
rz(0.83998799) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9962358) q[0];
sx q[0];
rz(-2.0780434) q[0];
sx q[0];
rz(2.9661687) q[0];
rz(-pi) q[1];
rz(2.5469668) q[2];
sx q[2];
rz(-1.42808) q[2];
sx q[2];
rz(-0.12446257) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.67761356) q[1];
sx q[1];
rz(-2.2871454) q[1];
sx q[1];
rz(0.78732995) q[1];
rz(3.1026353) q[3];
sx q[3];
rz(-1.0213199) q[3];
sx q[3];
rz(2.6770858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7175032) q[2];
sx q[2];
rz(-1.6987897) q[2];
sx q[2];
rz(-1.7355512) q[2];
rz(2.9794335) q[3];
sx q[3];
rz(-1.5608965) q[3];
sx q[3];
rz(0.46419188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9185987) q[0];
sx q[0];
rz(-0.080568947) q[0];
sx q[0];
rz(-2.719847) q[0];
rz(2.2987507) q[1];
sx q[1];
rz(-1.3858567) q[1];
sx q[1];
rz(0.27580321) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0378135) q[0];
sx q[0];
rz(-0.40291726) q[0];
sx q[0];
rz(-0.62312868) q[0];
x q[1];
rz(-0.47068542) q[2];
sx q[2];
rz(-1.3725201) q[2];
sx q[2];
rz(0.40775611) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2308064) q[1];
sx q[1];
rz(-1.2503997) q[1];
sx q[1];
rz(-1.5159392) q[1];
x q[2];
rz(2.2725676) q[3];
sx q[3];
rz(-1.4956253) q[3];
sx q[3];
rz(-1.6628671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3383823) q[2];
sx q[2];
rz(-0.79757491) q[2];
sx q[2];
rz(-2.4513643) q[2];
rz(0.35987443) q[3];
sx q[3];
rz(-1.7706324) q[3];
sx q[3];
rz(-0.38351044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91367078) q[0];
sx q[0];
rz(-1.0424732) q[0];
sx q[0];
rz(0.87164718) q[0];
rz(-2.7323885) q[1];
sx q[1];
rz(-2.6458461) q[1];
sx q[1];
rz(-0.31235487) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46074122) q[0];
sx q[0];
rz(-0.4609403) q[0];
sx q[0];
rz(0.26682202) q[0];
rz(2.9351685) q[2];
sx q[2];
rz(-1.4731493) q[2];
sx q[2];
rz(-0.7266149) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4620381) q[1];
sx q[1];
rz(-1.3674539) q[1];
sx q[1];
rz(2.773079) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7658224) q[3];
sx q[3];
rz(-2.2100976) q[3];
sx q[3];
rz(-2.2488058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.13545869) q[2];
sx q[2];
rz(-1.6929071) q[2];
sx q[2];
rz(-0.82497605) q[2];
rz(1.4247591) q[3];
sx q[3];
rz(-2.8202839) q[3];
sx q[3];
rz(0.43535522) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.41998) q[0];
sx q[0];
rz(-1.5897607) q[0];
sx q[0];
rz(2.2671674) q[0];
rz(-2.8127316) q[1];
sx q[1];
rz(-1.7560274) q[1];
sx q[1];
rz(1.0985451) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6499574) q[0];
sx q[0];
rz(-1.1960317) q[0];
sx q[0];
rz(-1.4521506) q[0];
rz(-pi) q[1];
x q[1];
rz(0.036089049) q[2];
sx q[2];
rz(-0.83613013) q[2];
sx q[2];
rz(-0.058908894) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.23197996) q[1];
sx q[1];
rz(-1.2709588) q[1];
sx q[1];
rz(-2.0832056) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9807113) q[3];
sx q[3];
rz(-1.8030589) q[3];
sx q[3];
rz(1.6764318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.99981368) q[2];
sx q[2];
rz(-1.1835316) q[2];
sx q[2];
rz(-2.4675274) q[2];
rz(1.6541803) q[3];
sx q[3];
rz(-1.4451278) q[3];
sx q[3];
rz(-0.28356799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11033002) q[0];
sx q[0];
rz(-1.0227579) q[0];
sx q[0];
rz(-1.3856101) q[0];
rz(-2.022187) q[1];
sx q[1];
rz(-1.6740572) q[1];
sx q[1];
rz(1.7562235) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.114217) q[0];
sx q[0];
rz(-1.8226591) q[0];
sx q[0];
rz(-2.2172539) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9846956) q[2];
sx q[2];
rz(-0.50373915) q[2];
sx q[2];
rz(-2.4258326) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0138356) q[1];
sx q[1];
rz(-2.7633939) q[1];
sx q[1];
rz(-2.5950123) q[1];
x q[2];
rz(-1.6072261) q[3];
sx q[3];
rz(-1.4163989) q[3];
sx q[3];
rz(1.2702219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0783656) q[2];
sx q[2];
rz(-2.207022) q[2];
sx q[2];
rz(-0.37609491) q[2];
rz(-1.1345351) q[3];
sx q[3];
rz(-0.36342707) q[3];
sx q[3];
rz(-0.80662066) q[3];
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
rz(3.1217839) q[0];
sx q[0];
rz(-2.5406295) q[0];
sx q[0];
rz(-3.0390749) q[0];
rz(-0.58492297) q[1];
sx q[1];
rz(-0.94766098) q[1];
sx q[1];
rz(-1.9580511) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85921209) q[0];
sx q[0];
rz(-1.6995653) q[0];
sx q[0];
rz(-0.14315258) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7223515) q[2];
sx q[2];
rz(-2.4346874) q[2];
sx q[2];
rz(2.8173994) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9283824) q[1];
sx q[1];
rz(-2.2021535) q[1];
sx q[1];
rz(2.8246095) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6449219) q[3];
sx q[3];
rz(-2.6210945) q[3];
sx q[3];
rz(2.6561873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.79503235) q[2];
sx q[2];
rz(-0.9298032) q[2];
sx q[2];
rz(-2.9417876) q[2];
rz(1.0605158) q[3];
sx q[3];
rz(-0.54148713) q[3];
sx q[3];
rz(-3.0919891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26315966) q[0];
sx q[0];
rz(-1.4819772) q[0];
sx q[0];
rz(2.8286381) q[0];
rz(-0.9494268) q[1];
sx q[1];
rz(-1.7890309) q[1];
sx q[1];
rz(1.4535646) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2672103) q[0];
sx q[0];
rz(-1.4802762) q[0];
sx q[0];
rz(0.1874013) q[0];
rz(-pi) q[1];
rz(1.4916999) q[2];
sx q[2];
rz(-2.1888615) q[2];
sx q[2];
rz(1.5072418) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1968699) q[1];
sx q[1];
rz(-2.2579282) q[1];
sx q[1];
rz(-1.8623167) q[1];
rz(-2.6864004) q[3];
sx q[3];
rz(-1.2641915) q[3];
sx q[3];
rz(2.2835177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.35385418) q[2];
sx q[2];
rz(-2.165803) q[2];
sx q[2];
rz(-1.3346416) q[2];
rz(1.3245964) q[3];
sx q[3];
rz(-0.66671222) q[3];
sx q[3];
rz(1.9273531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1052833) q[0];
sx q[0];
rz(-2.5258625) q[0];
sx q[0];
rz(-0.69806725) q[0];
rz(-2.2659194) q[1];
sx q[1];
rz(-2.0798637) q[1];
sx q[1];
rz(0.36044136) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9746285) q[0];
sx q[0];
rz(-0.95050838) q[0];
sx q[0];
rz(-1.0247158) q[0];
rz(-pi) q[1];
rz(-1.6351885) q[2];
sx q[2];
rz(-1.1887822) q[2];
sx q[2];
rz(-2.7329993) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0086509) q[1];
sx q[1];
rz(-1.9509778) q[1];
sx q[1];
rz(-0.93312414) q[1];
rz(-pi) q[2];
rz(-2.5927173) q[3];
sx q[3];
rz(-1.4272808) q[3];
sx q[3];
rz(1.9560069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2299819) q[2];
sx q[2];
rz(-1.7155827) q[2];
sx q[2];
rz(2.3243375) q[2];
rz(-2.3115555) q[3];
sx q[3];
rz(-0.54098141) q[3];
sx q[3];
rz(-0.219492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63477409) q[0];
sx q[0];
rz(-2.2932678) q[0];
sx q[0];
rz(-0.032055227) q[0];
rz(-1.3196779) q[1];
sx q[1];
rz(-1.7664884) q[1];
sx q[1];
rz(2.4874036) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35857707) q[0];
sx q[0];
rz(-0.23046432) q[0];
sx q[0];
rz(-0.5490659) q[0];
rz(2.0243353) q[2];
sx q[2];
rz(-0.83570601) q[2];
sx q[2];
rz(-2.105956) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5947551) q[1];
sx q[1];
rz(-2.0528626) q[1];
sx q[1];
rz(-1.0422816) q[1];
x q[2];
rz(-0.93864949) q[3];
sx q[3];
rz(-1.8018) q[3];
sx q[3];
rz(-1.7836026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1303945) q[2];
sx q[2];
rz(-1.0426714) q[2];
sx q[2];
rz(-1.0065669) q[2];
rz(0.69342962) q[3];
sx q[3];
rz(-1.3207685) q[3];
sx q[3];
rz(-1.8600195) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4257767) q[0];
sx q[0];
rz(-2.4632813) q[0];
sx q[0];
rz(3.0671469) q[0];
rz(-0.88059942) q[1];
sx q[1];
rz(-1.1136628) q[1];
sx q[1];
rz(-2.640092) q[1];
rz(-0.87133741) q[2];
sx q[2];
rz(-1.9819145) q[2];
sx q[2];
rz(-2.869538) q[2];
rz(1.1893336) q[3];
sx q[3];
rz(-0.91100024) q[3];
sx q[3];
rz(-0.80445214) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
