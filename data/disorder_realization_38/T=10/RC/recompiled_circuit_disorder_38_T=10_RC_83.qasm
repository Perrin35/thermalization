OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71198553) q[0];
sx q[0];
rz(-2.7349732) q[0];
sx q[0];
rz(-0.24917319) q[0];
rz(3.0781526) q[1];
sx q[1];
rz(-0.97172207) q[1];
sx q[1];
rz(2.5914153) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8367856) q[0];
sx q[0];
rz(-1.7903622) q[0];
sx q[0];
rz(3.1138793) q[0];
rz(-2.2416441) q[2];
sx q[2];
rz(-2.4953825) q[2];
sx q[2];
rz(0.83713573) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5397415) q[1];
sx q[1];
rz(-1.1907693) q[1];
sx q[1];
rz(-2.3439581) q[1];
x q[2];
rz(0.95991858) q[3];
sx q[3];
rz(-1.5023408) q[3];
sx q[3];
rz(1.5271036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.41574079) q[2];
sx q[2];
rz(-0.44885138) q[2];
sx q[2];
rz(0.63981167) q[2];
rz(0.85302991) q[3];
sx q[3];
rz(-0.60522389) q[3];
sx q[3];
rz(0.38133347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8752276) q[0];
sx q[0];
rz(-1.9748283) q[0];
sx q[0];
rz(-0.27045989) q[0];
rz(0.71331435) q[1];
sx q[1];
rz(-1.0353054) q[1];
sx q[1];
rz(1.6289904) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0892031) q[0];
sx q[0];
rz(-2.4501778) q[0];
sx q[0];
rz(-2.0877439) q[0];
rz(-pi) q[1];
rz(-2.8578051) q[2];
sx q[2];
rz(-1.5675401) q[2];
sx q[2];
rz(-0.46797215) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0257033) q[1];
sx q[1];
rz(-1.70417) q[1];
sx q[1];
rz(-0.89110903) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2849502) q[3];
sx q[3];
rz(-1.8010406) q[3];
sx q[3];
rz(-2.2765991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9397395) q[2];
sx q[2];
rz(-2.8994603) q[2];
sx q[2];
rz(0.80336037) q[2];
rz(-2.0837636) q[3];
sx q[3];
rz(-1.4927031) q[3];
sx q[3];
rz(-0.025618205) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8818883) q[0];
sx q[0];
rz(-1.0972247) q[0];
sx q[0];
rz(-2.846068) q[0];
rz(2.9064536) q[1];
sx q[1];
rz(-1.7087015) q[1];
sx q[1];
rz(2.3957516) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64622067) q[0];
sx q[0];
rz(-1.6618068) q[0];
sx q[0];
rz(-1.3489086) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27772851) q[2];
sx q[2];
rz(-0.54364294) q[2];
sx q[2];
rz(-3.1090528) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6927166) q[1];
sx q[1];
rz(-1.2926896) q[1];
sx q[1];
rz(2.6443308) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3219222) q[3];
sx q[3];
rz(-1.1296774) q[3];
sx q[3];
rz(-0.53019023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0535584) q[2];
sx q[2];
rz(-3.1159248) q[2];
sx q[2];
rz(2.4528465) q[2];
rz(-3.0912494) q[3];
sx q[3];
rz(-2.2294932) q[3];
sx q[3];
rz(-1.5215727) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21928366) q[0];
sx q[0];
rz(-1.0204717) q[0];
sx q[0];
rz(3.0396089) q[0];
rz(0.12022262) q[1];
sx q[1];
rz(-0.52821237) q[1];
sx q[1];
rz(2.8682958) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2647588) q[0];
sx q[0];
rz(-0.33565531) q[0];
sx q[0];
rz(1.8434974) q[0];
rz(-pi) q[1];
rz(0.64812135) q[2];
sx q[2];
rz(-2.0378049) q[2];
sx q[2];
rz(-1.8304706) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8371493) q[1];
sx q[1];
rz(-1.4238365) q[1];
sx q[1];
rz(1.8069581) q[1];
rz(1.7498383) q[3];
sx q[3];
rz(-0.92902196) q[3];
sx q[3];
rz(0.87953506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3499202) q[2];
sx q[2];
rz(-0.73799729) q[2];
sx q[2];
rz(-2.6376574) q[2];
rz(3.062011) q[3];
sx q[3];
rz(-1.1527529) q[3];
sx q[3];
rz(-2.7664405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85916096) q[0];
sx q[0];
rz(-2.0895884) q[0];
sx q[0];
rz(0.082745634) q[0];
rz(-2.4619608) q[1];
sx q[1];
rz(-1.6487164) q[1];
sx q[1];
rz(0.98714978) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3521096) q[0];
sx q[0];
rz(-1.1687359) q[0];
sx q[0];
rz(2.2228918) q[0];
x q[1];
rz(1.5613902) q[2];
sx q[2];
rz(-0.9270037) q[2];
sx q[2];
rz(-0.26088342) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8122711) q[1];
sx q[1];
rz(-1.2580039) q[1];
sx q[1];
rz(-0.66004628) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4694674) q[3];
sx q[3];
rz(-0.39468995) q[3];
sx q[3];
rz(-2.8916388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2510898) q[2];
sx q[2];
rz(-0.7901929) q[2];
sx q[2];
rz(0.21128543) q[2];
rz(0.42090297) q[3];
sx q[3];
rz(-0.55287164) q[3];
sx q[3];
rz(3.0781854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(-1.485065) q[0];
sx q[0];
rz(-1.2542897) q[0];
sx q[0];
rz(2.3497537) q[0];
rz(0.99545288) q[1];
sx q[1];
rz(-0.96770006) q[1];
sx q[1];
rz(1.2794367) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0148221) q[0];
sx q[0];
rz(-1.556067) q[0];
sx q[0];
rz(-1.8526088) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5336669) q[2];
sx q[2];
rz(-1.330901) q[2];
sx q[2];
rz(-1.243967) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.8705604) q[1];
sx q[1];
rz(-0.522627) q[1];
sx q[1];
rz(-1.6305627) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7518696) q[3];
sx q[3];
rz(-1.7098134) q[3];
sx q[3];
rz(-2.3982323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.25097686) q[2];
sx q[2];
rz(-1.0303409) q[2];
sx q[2];
rz(0.77077579) q[2];
rz(1.6714913) q[3];
sx q[3];
rz(-0.41430587) q[3];
sx q[3];
rz(-2.9582086) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32271785) q[0];
sx q[0];
rz(-0.6495496) q[0];
sx q[0];
rz(0.055667002) q[0];
rz(0.21559134) q[1];
sx q[1];
rz(-0.76342738) q[1];
sx q[1];
rz(-0.0035704426) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17655003) q[0];
sx q[0];
rz(-1.9854443) q[0];
sx q[0];
rz(1.167017) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7262906) q[2];
sx q[2];
rz(-0.79681444) q[2];
sx q[2];
rz(2.08564) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0842972) q[1];
sx q[1];
rz(-2.1665386) q[1];
sx q[1];
rz(2.7792395) q[1];
x q[2];
rz(3.0655541) q[3];
sx q[3];
rz(-1.7685316) q[3];
sx q[3];
rz(-2.7969489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5443762) q[2];
sx q[2];
rz(-0.95149136) q[2];
sx q[2];
rz(-0.92010951) q[2];
rz(0.19872935) q[3];
sx q[3];
rz(-1.9072429) q[3];
sx q[3];
rz(-0.41771093) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51628095) q[0];
sx q[0];
rz(-1.6315062) q[0];
sx q[0];
rz(-1.7238808) q[0];
rz(2.7334546) q[1];
sx q[1];
rz(-1.0393655) q[1];
sx q[1];
rz(2.4628941) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1269826) q[0];
sx q[0];
rz(-1.8982366) q[0];
sx q[0];
rz(-1.8493269) q[0];
rz(-pi) q[1];
rz(2.3085824) q[2];
sx q[2];
rz(-1.2150803) q[2];
sx q[2];
rz(-0.21258159) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7486836) q[1];
sx q[1];
rz(-0.52782413) q[1];
sx q[1];
rz(1.4455568) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4123165) q[3];
sx q[3];
rz(-1.1500119) q[3];
sx q[3];
rz(-2.5077523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.17710182) q[2];
sx q[2];
rz(-1.0937546) q[2];
sx q[2];
rz(-0.91782451) q[2];
rz(-1.9994036) q[3];
sx q[3];
rz(-0.85421383) q[3];
sx q[3];
rz(-0.93723047) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98638242) q[0];
sx q[0];
rz(-0.6739524) q[0];
sx q[0];
rz(0.25979364) q[0];
rz(2.4329176) q[1];
sx q[1];
rz(-2.8627113) q[1];
sx q[1];
rz(-0.52694595) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.066756847) q[0];
sx q[0];
rz(-1.2196676) q[0];
sx q[0];
rz(-1.3394651) q[0];
x q[1];
rz(0.2692659) q[2];
sx q[2];
rz(-2.2109647) q[2];
sx q[2];
rz(-2.7562041) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7582015) q[1];
sx q[1];
rz(-1.2609298) q[1];
sx q[1];
rz(1.8499225) q[1];
rz(-pi) q[2];
x q[2];
rz(1.280904) q[3];
sx q[3];
rz(-1.1416417) q[3];
sx q[3];
rz(3.1115301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8156585) q[2];
sx q[2];
rz(-0.84647536) q[2];
sx q[2];
rz(2.3507067) q[2];
rz(-2.7549426) q[3];
sx q[3];
rz(-1.0145885) q[3];
sx q[3];
rz(2.8044243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40795046) q[0];
sx q[0];
rz(-2.9691417) q[0];
sx q[0];
rz(-2.1561484) q[0];
rz(-2.573029) q[1];
sx q[1];
rz(-2.0303346) q[1];
sx q[1];
rz(-2.7808166) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11689582) q[0];
sx q[0];
rz(-1.403406) q[0];
sx q[0];
rz(0.025339729) q[0];
x q[1];
rz(-0.12014328) q[2];
sx q[2];
rz(-1.7842245) q[2];
sx q[2];
rz(-0.10211589) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1211179) q[1];
sx q[1];
rz(-1.7133683) q[1];
sx q[1];
rz(-2.3974182) q[1];
rz(1.4233227) q[3];
sx q[3];
rz(-1.7383988) q[3];
sx q[3];
rz(-0.22351219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0976022) q[2];
sx q[2];
rz(-2.4520935) q[2];
sx q[2];
rz(2.1981751) q[2];
rz(-0.57389456) q[3];
sx q[3];
rz(-0.4807764) q[3];
sx q[3];
rz(1.3963612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5744793) q[0];
sx q[0];
rz(-1.4470826) q[0];
sx q[0];
rz(-0.8599109) q[0];
rz(-1.3600596) q[1];
sx q[1];
rz(-2.139745) q[1];
sx q[1];
rz(2.3812961) q[1];
rz(1.3654937) q[2];
sx q[2];
rz(-1.6663972) q[2];
sx q[2];
rz(-1.321928) q[2];
rz(-2.4102224) q[3];
sx q[3];
rz(-0.89805713) q[3];
sx q[3];
rz(-3.0039207) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];