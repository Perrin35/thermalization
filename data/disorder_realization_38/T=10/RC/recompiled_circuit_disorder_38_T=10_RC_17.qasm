OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4296071) q[0];
sx q[0];
rz(-0.40661943) q[0];
sx q[0];
rz(-2.8924195) q[0];
rz(3.0781526) q[1];
sx q[1];
rz(-0.97172207) q[1];
sx q[1];
rz(2.5914153) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30480706) q[0];
sx q[0];
rz(-1.3512304) q[0];
sx q[0];
rz(3.1138793) q[0];
rz(-pi) q[1];
rz(-1.0371738) q[2];
sx q[2];
rz(-1.1871157) q[2];
sx q[2];
rz(-1.2984315) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8093811) q[1];
sx q[1];
rz(-0.84377938) q[1];
sx q[1];
rz(1.0512645) q[1];
rz(2.1816741) q[3];
sx q[3];
rz(-1.6392518) q[3];
sx q[3];
rz(-1.6144891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7258519) q[2];
sx q[2];
rz(-2.6927413) q[2];
sx q[2];
rz(0.63981167) q[2];
rz(-2.2885627) q[3];
sx q[3];
rz(-0.60522389) q[3];
sx q[3];
rz(0.38133347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(1.8752276) q[0];
sx q[0];
rz(-1.1667644) q[0];
sx q[0];
rz(-0.27045989) q[0];
rz(-2.4282783) q[1];
sx q[1];
rz(-1.0353054) q[1];
sx q[1];
rz(-1.5126022) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7249811) q[0];
sx q[0];
rz(-2.1583301) q[0];
sx q[0];
rz(0.38831098) q[0];
rz(0.28378758) q[2];
sx q[2];
rz(-1.5740526) q[2];
sx q[2];
rz(0.46797215) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4334129) q[1];
sx q[1];
rz(-0.69060329) q[1];
sx q[1];
rz(1.7811) q[1];
rz(-2.9019722) q[3];
sx q[3];
rz(-1.8488956) q[3];
sx q[3];
rz(-0.77277377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9397395) q[2];
sx q[2];
rz(-2.8994603) q[2];
sx q[2];
rz(-2.3382323) q[2];
rz(1.057829) q[3];
sx q[3];
rz(-1.4927031) q[3];
sx q[3];
rz(-0.025618205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8818883) q[0];
sx q[0];
rz(-2.0443679) q[0];
sx q[0];
rz(-0.29552466) q[0];
rz(-2.9064536) q[1];
sx q[1];
rz(-1.4328911) q[1];
sx q[1];
rz(-0.74584109) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.495372) q[0];
sx q[0];
rz(-1.4797858) q[0];
sx q[0];
rz(-1.3489086) q[0];
x q[1];
rz(2.6150871) q[2];
sx q[2];
rz(-1.4284992) q[2];
sx q[2];
rz(-1.3640179) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5513788) q[1];
sx q[1];
rz(-2.5775902) q[1];
sx q[1];
rz(-0.53932921) q[1];
x q[2];
rz(2.6607473) q[3];
sx q[3];
rz(-0.50243176) q[3];
sx q[3];
rz(-1.06711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.088034257) q[2];
sx q[2];
rz(-3.1159248) q[2];
sx q[2];
rz(0.68874613) q[2];
rz(-3.0912494) q[3];
sx q[3];
rz(-0.91209948) q[3];
sx q[3];
rz(-1.6200199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21928366) q[0];
sx q[0];
rz(-1.0204717) q[0];
sx q[0];
rz(3.0396089) q[0];
rz(3.02137) q[1];
sx q[1];
rz(-0.52821237) q[1];
sx q[1];
rz(0.27329683) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8768339) q[0];
sx q[0];
rz(-0.33565531) q[0];
sx q[0];
rz(-1.8434974) q[0];
rz(-2.4457473) q[2];
sx q[2];
rz(-0.77866422) q[2];
sx q[2];
rz(0.7960745) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.30444333) q[1];
sx q[1];
rz(-1.4238365) q[1];
sx q[1];
rz(1.8069581) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23394211) q[3];
sx q[3];
rz(-2.4787239) q[3];
sx q[3];
rz(-2.5556504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3499202) q[2];
sx q[2];
rz(-2.4035954) q[2];
sx q[2];
rz(2.6376574) q[2];
rz(-0.079581633) q[3];
sx q[3];
rz(-1.9888398) q[3];
sx q[3];
rz(-0.3751522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85916096) q[0];
sx q[0];
rz(-1.0520042) q[0];
sx q[0];
rz(3.058847) q[0];
rz(2.4619608) q[1];
sx q[1];
rz(-1.6487164) q[1];
sx q[1];
rz(-0.98714978) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25492451) q[0];
sx q[0];
rz(-0.75037557) q[0];
sx q[0];
rz(-0.95959856) q[0];
rz(-pi) q[1];
rz(-1.5613902) q[2];
sx q[2];
rz(-0.9270037) q[2];
sx q[2];
rz(0.26088342) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.6189177) q[1];
sx q[1];
rz(-2.4213311) q[1];
sx q[1];
rz(-0.48536761) q[1];
rz(3.0994814) q[3];
sx q[3];
rz(-1.9633506) q[3];
sx q[3];
rz(0.3596572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2510898) q[2];
sx q[2];
rz(-2.3513998) q[2];
sx q[2];
rz(-0.21128543) q[2];
rz(-2.7206897) q[3];
sx q[3];
rz(-2.588721) q[3];
sx q[3];
rz(-3.0781854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6565276) q[0];
sx q[0];
rz(-1.2542897) q[0];
sx q[0];
rz(2.3497537) q[0];
rz(-2.1461398) q[1];
sx q[1];
rz(-2.1738926) q[1];
sx q[1];
rz(1.8621559) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1267705) q[0];
sx q[0];
rz(-1.556067) q[0];
sx q[0];
rz(-1.2889839) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.860414) q[2];
sx q[2];
rz(-0.98266232) q[2];
sx q[2];
rz(-2.9786125) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2020859) q[1];
sx q[1];
rz(-1.0491976) q[1];
sx q[1];
rz(-0.034394666) q[1];
rz(-pi) q[2];
rz(1.389723) q[3];
sx q[3];
rz(-1.7098134) q[3];
sx q[3];
rz(2.3982323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8906158) q[2];
sx q[2];
rz(-1.0303409) q[2];
sx q[2];
rz(-2.3708169) q[2];
rz(1.4701014) q[3];
sx q[3];
rz(-0.41430587) q[3];
sx q[3];
rz(2.9582086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8188748) q[0];
sx q[0];
rz(-2.4920431) q[0];
sx q[0];
rz(-3.0859257) q[0];
rz(0.21559134) q[1];
sx q[1];
rz(-2.3781653) q[1];
sx q[1];
rz(-3.1380222) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5647033) q[0];
sx q[0];
rz(-1.2029552) q[0];
sx q[0];
rz(-2.6951615) q[0];
rz(-2.7262906) q[2];
sx q[2];
rz(-0.79681444) q[2];
sx q[2];
rz(-1.0559527) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0842972) q[1];
sx q[1];
rz(-2.1665386) q[1];
sx q[1];
rz(2.7792395) q[1];
rz(-1.9332063) q[3];
sx q[3];
rz(-0.21167314) q[3];
sx q[3];
rz(-0.025312245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.59721649) q[2];
sx q[2];
rz(-2.1901013) q[2];
sx q[2];
rz(2.2214831) q[2];
rz(0.19872935) q[3];
sx q[3];
rz(-1.9072429) q[3];
sx q[3];
rz(-0.41771093) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51628095) q[0];
sx q[0];
rz(-1.5100864) q[0];
sx q[0];
rz(1.4177119) q[0];
rz(-2.7334546) q[1];
sx q[1];
rz(-2.1022271) q[1];
sx q[1];
rz(-0.67869854) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7413901) q[0];
sx q[0];
rz(-2.7149902) q[0];
sx q[0];
rz(-2.461117) q[0];
rz(-1.0661725) q[2];
sx q[2];
rz(-2.3373211) q[2];
sx q[2];
rz(-1.7240766) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8933967) q[1];
sx q[1];
rz(-2.0940603) q[1];
sx q[1];
rz(-3.0688973) q[1];
x q[2];
rz(-0.42550605) q[3];
sx q[3];
rz(-1.4262428) q[3];
sx q[3];
rz(-1.0021462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9644908) q[2];
sx q[2];
rz(-1.0937546) q[2];
sx q[2];
rz(-2.2237681) q[2];
rz(1.142189) q[3];
sx q[3];
rz(-0.85421383) q[3];
sx q[3];
rz(-0.93723047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98638242) q[0];
sx q[0];
rz(-2.4676403) q[0];
sx q[0];
rz(-2.881799) q[0];
rz(-0.70867509) q[1];
sx q[1];
rz(-2.8627113) q[1];
sx q[1];
rz(2.6146467) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0748358) q[0];
sx q[0];
rz(-1.9219251) q[0];
sx q[0];
rz(-1.8021276) q[0];
x q[1];
rz(-1.9138463) q[2];
sx q[2];
rz(-2.4545049) q[2];
sx q[2];
rz(-2.3234141) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7582015) q[1];
sx q[1];
rz(-1.8806629) q[1];
sx q[1];
rz(-1.8499225) q[1];
x q[2];
rz(0.55836375) q[3];
sx q[3];
rz(-0.51279587) q[3];
sx q[3];
rz(0.59190291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8156585) q[2];
sx q[2];
rz(-2.2951173) q[2];
sx q[2];
rz(-2.3507067) q[2];
rz(-2.7549426) q[3];
sx q[3];
rz(-2.1270042) q[3];
sx q[3];
rz(-2.8044243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40795046) q[0];
sx q[0];
rz(-2.9691417) q[0];
sx q[0];
rz(0.98544425) q[0];
rz(0.5685637) q[1];
sx q[1];
rz(-2.0303346) q[1];
sx q[1];
rz(0.3607761) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11689582) q[0];
sx q[0];
rz(-1.7381867) q[0];
sx q[0];
rz(-0.025339729) q[0];
rz(1.3558657) q[2];
sx q[2];
rz(-1.6882009) q[2];
sx q[2];
rz(-1.4431151) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0204748) q[1];
sx q[1];
rz(-1.4282244) q[1];
sx q[1];
rz(2.3974182) q[1];
rz(-1.4233227) q[3];
sx q[3];
rz(-1.4031938) q[3];
sx q[3];
rz(-0.22351219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0976022) q[2];
sx q[2];
rz(-0.68949914) q[2];
sx q[2];
rz(-0.94341755) q[2];
rz(2.5676981) q[3];
sx q[3];
rz(-0.4807764) q[3];
sx q[3];
rz(1.3963612) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5744793) q[0];
sx q[0];
rz(-1.6945101) q[0];
sx q[0];
rz(2.2816818) q[0];
rz(1.7815331) q[1];
sx q[1];
rz(-2.139745) q[1];
sx q[1];
rz(2.3812961) q[1];
rz(-1.1311244) q[2];
sx q[2];
rz(-0.22618539) q[2];
sx q[2];
rz(-2.4629081) q[2];
rz(0.87313575) q[3];
sx q[3];
rz(-0.94948873) q[3];
sx q[3];
rz(-0.82617847) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
