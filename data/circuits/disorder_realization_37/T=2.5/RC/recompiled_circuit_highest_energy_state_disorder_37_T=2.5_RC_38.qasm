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
rz(-1.9395186) q[0];
sx q[0];
rz(2.095686) q[0];
sx q[0];
rz(8.6742824) q[0];
rz(0.22076386) q[1];
sx q[1];
rz(-2.0808487) q[1];
sx q[1];
rz(1.9579252) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96776828) q[0];
sx q[0];
rz(-2.2319751) q[0];
sx q[0];
rz(1.038616) q[0];
x q[1];
rz(-0.12678876) q[2];
sx q[2];
rz(-1.7361963) q[2];
sx q[2];
rz(-0.70387712) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.88088502) q[1];
sx q[1];
rz(-1.5591677) q[1];
sx q[1];
rz(-0.67990085) q[1];
rz(-pi) q[2];
rz(1.2624083) q[3];
sx q[3];
rz(-1.3779791) q[3];
sx q[3];
rz(0.28379019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3262647) q[2];
sx q[2];
rz(-1.999819) q[2];
sx q[2];
rz(-3.0495138) q[2];
rz(2.0103256) q[3];
sx q[3];
rz(-2.0720033) q[3];
sx q[3];
rz(-0.85589516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49694127) q[0];
sx q[0];
rz(-2.6583789) q[0];
sx q[0];
rz(0.53571969) q[0];
rz(-3.1247395) q[1];
sx q[1];
rz(-2.2622175) q[1];
sx q[1];
rz(-3.1244315) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099875256) q[0];
sx q[0];
rz(-1.2493396) q[0];
sx q[0];
rz(-1.0519371) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.13249915) q[2];
sx q[2];
rz(-2.6713058) q[2];
sx q[2];
rz(0.4343701) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29083272) q[1];
sx q[1];
rz(-2.8597745) q[1];
sx q[1];
rz(2.5435104) q[1];
rz(-pi) q[2];
rz(3.0073418) q[3];
sx q[3];
rz(-1.6608616) q[3];
sx q[3];
rz(1.4081189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0993593) q[2];
sx q[2];
rz(-1.4734273) q[2];
sx q[2];
rz(-2.8933706) q[2];
rz(0.92784268) q[3];
sx q[3];
rz(-0.78301269) q[3];
sx q[3];
rz(-0.46970126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27702734) q[0];
sx q[0];
rz(-1.1772573) q[0];
sx q[0];
rz(0.16265854) q[0];
rz(-2.3795369) q[1];
sx q[1];
rz(-2.9167852) q[1];
sx q[1];
rz(2.3077097) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57788173) q[0];
sx q[0];
rz(-1.0037918) q[0];
sx q[0];
rz(1.4253229) q[0];
rz(0.32270785) q[2];
sx q[2];
rz(-1.4526748) q[2];
sx q[2];
rz(-0.28216991) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2411338) q[1];
sx q[1];
rz(-0.83158439) q[1];
sx q[1];
rz(0.99947387) q[1];
rz(-0.44077222) q[3];
sx q[3];
rz(-1.5197872) q[3];
sx q[3];
rz(2.9875502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.062408058) q[2];
sx q[2];
rz(-2.0024313) q[2];
sx q[2];
rz(2.8774234) q[2];
rz(2.0832113) q[3];
sx q[3];
rz(-0.94401413) q[3];
sx q[3];
rz(3.1112166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77085483) q[0];
sx q[0];
rz(-2.3069032) q[0];
sx q[0];
rz(-1.4282861) q[0];
rz(-0.76438534) q[1];
sx q[1];
rz(-1.9995707) q[1];
sx q[1];
rz(1.8199325) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7390838) q[0];
sx q[0];
rz(-0.54259491) q[0];
sx q[0];
rz(0.71904166) q[0];
rz(-pi) q[1];
rz(-1.7197505) q[2];
sx q[2];
rz(-2.8064054) q[2];
sx q[2];
rz(-2.6870948) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.053406) q[1];
sx q[1];
rz(-1.2866396) q[1];
sx q[1];
rz(-2.8233145) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2854487) q[3];
sx q[3];
rz(-1.2309181) q[3];
sx q[3];
rz(1.651317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4798212) q[2];
sx q[2];
rz(-1.4134553) q[2];
sx q[2];
rz(-2.856355) q[2];
rz(-1.8789004) q[3];
sx q[3];
rz(-1.9243536) q[3];
sx q[3];
rz(1.4234022) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30094639) q[0];
sx q[0];
rz(-2.3941289) q[0];
sx q[0];
rz(-0.88229156) q[0];
rz(2.4709985) q[1];
sx q[1];
rz(-2.0457485) q[1];
sx q[1];
rz(-0.98463279) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1272972) q[0];
sx q[0];
rz(-1.5080351) q[0];
sx q[0];
rz(1.5393637) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1124723) q[2];
sx q[2];
rz(-1.6075252) q[2];
sx q[2];
rz(-2.6482895) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4337816) q[1];
sx q[1];
rz(-2.1137775) q[1];
sx q[1];
rz(2.19853) q[1];
x q[2];
rz(0.71159848) q[3];
sx q[3];
rz(-0.67005537) q[3];
sx q[3];
rz(1.3942476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.275445) q[2];
sx q[2];
rz(-2.6068942) q[2];
sx q[2];
rz(2.5246485) q[2];
rz(0.75211891) q[3];
sx q[3];
rz(-1.3277466) q[3];
sx q[3];
rz(0.45929685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1521456) q[0];
sx q[0];
rz(-2.4247657) q[0];
sx q[0];
rz(-0.77051198) q[0];
rz(0.72692263) q[1];
sx q[1];
rz(-0.22218552) q[1];
sx q[1];
rz(-2.4368584) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4024538) q[0];
sx q[0];
rz(-1.0346892) q[0];
sx q[0];
rz(0.86748634) q[0];
rz(0.36811916) q[2];
sx q[2];
rz(-1.0794753) q[2];
sx q[2];
rz(1.4368601) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.875346) q[1];
sx q[1];
rz(-0.35142144) q[1];
sx q[1];
rz(-0.09697341) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1416311) q[3];
sx q[3];
rz(-1.785716) q[3];
sx q[3];
rz(-2.8254912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7695693) q[2];
sx q[2];
rz(-1.4906887) q[2];
sx q[2];
rz(-2.7083569) q[2];
rz(-2.9406934) q[3];
sx q[3];
rz(-1.9379987) q[3];
sx q[3];
rz(-2.6065629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.83605003) q[0];
sx q[0];
rz(-1.9925646) q[0];
sx q[0];
rz(-0.38247633) q[0];
rz(1.0410694) q[1];
sx q[1];
rz(-1.46547) q[1];
sx q[1];
rz(2.6952851) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9093626) q[0];
sx q[0];
rz(-1.5188471) q[0];
sx q[0];
rz(-2.8828017) q[0];
rz(-pi) q[1];
rz(1.3417753) q[2];
sx q[2];
rz(-1.9441368) q[2];
sx q[2];
rz(0.41002204) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1416229) q[1];
sx q[1];
rz(-1.9239307) q[1];
sx q[1];
rz(0.94393877) q[1];
rz(-pi) q[2];
rz(-0.71766149) q[3];
sx q[3];
rz(-0.98472408) q[3];
sx q[3];
rz(-1.7505898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3216766) q[2];
sx q[2];
rz(-2.0240462) q[2];
sx q[2];
rz(-0.72810158) q[2];
rz(2.8150832) q[3];
sx q[3];
rz(-1.6312586) q[3];
sx q[3];
rz(-2.2842893) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2364872) q[0];
sx q[0];
rz(-0.40920722) q[0];
sx q[0];
rz(-1.8439199) q[0];
rz(-2.1406651) q[1];
sx q[1];
rz(-2.3960254) q[1];
sx q[1];
rz(-2.7330858) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4057874) q[0];
sx q[0];
rz(-1.5716432) q[0];
sx q[0];
rz(1.5566096) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0922008) q[2];
sx q[2];
rz(-1.166162) q[2];
sx q[2];
rz(1.7737845) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0169819) q[1];
sx q[1];
rz(-1.4370185) q[1];
sx q[1];
rz(0.11511187) q[1];
x q[2];
rz(2.767603) q[3];
sx q[3];
rz(-1.9714103) q[3];
sx q[3];
rz(-3.0246322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0268176) q[2];
sx q[2];
rz(-1.5207745) q[2];
sx q[2];
rz(-2.3180023) q[2];
rz(0.94281998) q[3];
sx q[3];
rz(-1.4225682) q[3];
sx q[3];
rz(-1.6387117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5470062) q[0];
sx q[0];
rz(-0.44438812) q[0];
sx q[0];
rz(1.8894926) q[0];
rz(1.7777255) q[1];
sx q[1];
rz(-1.6412647) q[1];
sx q[1];
rz(1.3346671) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8060243) q[0];
sx q[0];
rz(-2.8360406) q[0];
sx q[0];
rz(-3.0204885) q[0];
rz(-pi) q[1];
rz(-0.23128839) q[2];
sx q[2];
rz(-2.4866001) q[2];
sx q[2];
rz(-0.97823036) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4542983) q[1];
sx q[1];
rz(-2.0210938) q[1];
sx q[1];
rz(-1.4895093) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1030508) q[3];
sx q[3];
rz(-1.333101) q[3];
sx q[3];
rz(2.2464744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.021726457) q[2];
sx q[2];
rz(-2.848337) q[2];
sx q[2];
rz(-0.29652706) q[2];
rz(-1.4793388) q[3];
sx q[3];
rz(-1.5115503) q[3];
sx q[3];
rz(-0.1571981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6327561) q[0];
sx q[0];
rz(-2.4907676) q[0];
sx q[0];
rz(0.79488361) q[0];
rz(-2.5859313) q[1];
sx q[1];
rz(-1.2763005) q[1];
sx q[1];
rz(1.4478252) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5710766) q[0];
sx q[0];
rz(-2.2865642) q[0];
sx q[0];
rz(1.060673) q[0];
x q[1];
rz(-3.0270711) q[2];
sx q[2];
rz(-2.6410069) q[2];
sx q[2];
rz(1.3040744) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70778479) q[1];
sx q[1];
rz(-2.6205739) q[1];
sx q[1];
rz(1.9791114) q[1];
rz(-0.9742188) q[3];
sx q[3];
rz(-2.0250626) q[3];
sx q[3];
rz(2.5429986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.29622233) q[2];
sx q[2];
rz(-1.5105379) q[2];
sx q[2];
rz(-1.4307865) q[2];
rz(0.75302643) q[3];
sx q[3];
rz(-0.81792653) q[3];
sx q[3];
rz(2.9770765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13567781) q[0];
sx q[0];
rz(-2.4527241) q[0];
sx q[0];
rz(1.2149568) q[0];
rz(1.6750492) q[1];
sx q[1];
rz(-2.2129682) q[1];
sx q[1];
rz(0.60563544) q[1];
rz(0.56500021) q[2];
sx q[2];
rz(-1.0641268) q[2];
sx q[2];
rz(0.63882154) q[2];
rz(-1.4004272) q[3];
sx q[3];
rz(-1.4650678) q[3];
sx q[3];
rz(2.3847945) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
