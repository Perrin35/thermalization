OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.31792274) q[0];
sx q[0];
rz(-1.5225141) q[0];
sx q[0];
rz(3.0784399) q[0];
rz(1.8237279) q[1];
sx q[1];
rz(-0.39931077) q[1];
sx q[1];
rz(2.3656486) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5781893) q[0];
sx q[0];
rz(-3.0886123) q[0];
sx q[0];
rz(0.53081675) q[0];
rz(-pi) q[1];
rz(2.2047441) q[2];
sx q[2];
rz(-2.4000451) q[2];
sx q[2];
rz(-0.28683603) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9777955) q[1];
sx q[1];
rz(-1.3936636) q[1];
sx q[1];
rz(2.2832485) q[1];
rz(-pi) q[2];
rz(-1.8100061) q[3];
sx q[3];
rz(-2.0612217) q[3];
sx q[3];
rz(-1.8338721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.45106384) q[2];
sx q[2];
rz(-1.3961926) q[2];
sx q[2];
rz(-2.6032791) q[2];
rz(-2.8484143) q[3];
sx q[3];
rz(-1.5436951) q[3];
sx q[3];
rz(-1.258491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9341768) q[0];
sx q[0];
rz(-0.84686142) q[0];
sx q[0];
rz(2.9127981) q[0];
rz(-3.0711807) q[1];
sx q[1];
rz(-1.8682624) q[1];
sx q[1];
rz(2.0796897) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7477729) q[0];
sx q[0];
rz(-2.5508587) q[0];
sx q[0];
rz(0.48818199) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47669551) q[2];
sx q[2];
rz(-1.9856493) q[2];
sx q[2];
rz(1.7384993) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6390266) q[1];
sx q[1];
rz(-0.4430534) q[1];
sx q[1];
rz(0.73213099) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.44159378) q[3];
sx q[3];
rz(-0.46735763) q[3];
sx q[3];
rz(-0.87477126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.215302) q[2];
sx q[2];
rz(-2.6965202) q[2];
sx q[2];
rz(-0.11623795) q[2];
rz(2.2262573) q[3];
sx q[3];
rz(-1.4696308) q[3];
sx q[3];
rz(0.62469971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4794469) q[0];
sx q[0];
rz(-1.5597458) q[0];
sx q[0];
rz(-0.33552718) q[0];
rz(-1.7299995) q[1];
sx q[1];
rz(-0.84452191) q[1];
sx q[1];
rz(1.1988877) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6035164) q[0];
sx q[0];
rz(-1.8545574) q[0];
sx q[0];
rz(2.1233344) q[0];
x q[1];
rz(-2.6012318) q[2];
sx q[2];
rz(-0.66695222) q[2];
sx q[2];
rz(-1.6467384) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4863805) q[1];
sx q[1];
rz(-0.81334121) q[1];
sx q[1];
rz(1.4105303) q[1];
rz(-pi) q[2];
rz(-0.24925225) q[3];
sx q[3];
rz(-2.0468759) q[3];
sx q[3];
rz(-2.2152546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.18232839) q[2];
sx q[2];
rz(-0.79572833) q[2];
sx q[2];
rz(-2.0286593) q[2];
rz(-2.6304604) q[3];
sx q[3];
rz(-1.7685578) q[3];
sx q[3];
rz(-2.6326877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6540601) q[0];
sx q[0];
rz(-0.075981058) q[0];
sx q[0];
rz(-0.37387601) q[0];
rz(-1.182425) q[1];
sx q[1];
rz(-1.9805209) q[1];
sx q[1];
rz(2.9808796) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5406093) q[0];
sx q[0];
rz(-1.4330787) q[0];
sx q[0];
rz(1.4871554) q[0];
x q[1];
rz(0.66443759) q[2];
sx q[2];
rz(-0.91240806) q[2];
sx q[2];
rz(0.13641549) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0707145) q[1];
sx q[1];
rz(-1.9556112) q[1];
sx q[1];
rz(-1.0139731) q[1];
rz(2.3869208) q[3];
sx q[3];
rz(-0.42413482) q[3];
sx q[3];
rz(1.4671668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.60761991) q[2];
sx q[2];
rz(-2.0325568) q[2];
sx q[2];
rz(-0.71471659) q[2];
rz(2.886582) q[3];
sx q[3];
rz(-1.8111572) q[3];
sx q[3];
rz(-2.3431006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.052499972) q[0];
sx q[0];
rz(-0.30464259) q[0];
sx q[0];
rz(0.76328817) q[0];
rz(2.8870562) q[1];
sx q[1];
rz(-1.3724047) q[1];
sx q[1];
rz(1.6932142) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7970552) q[0];
sx q[0];
rz(-1.9112003) q[0];
sx q[0];
rz(-1.1974687) q[0];
rz(0.82008597) q[2];
sx q[2];
rz(-1.4561355) q[2];
sx q[2];
rz(-1.3263758) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9824286) q[1];
sx q[1];
rz(-0.64627534) q[1];
sx q[1];
rz(0.29781945) q[1];
rz(-0.71131018) q[3];
sx q[3];
rz(-1.8790725) q[3];
sx q[3];
rz(3.0570249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.86108565) q[2];
sx q[2];
rz(-2.0528767) q[2];
sx q[2];
rz(-2.1121934) q[2];
rz(1.5329125) q[3];
sx q[3];
rz(-2.2423988) q[3];
sx q[3];
rz(-0.6338343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3935881) q[0];
sx q[0];
rz(-0.91359502) q[0];
sx q[0];
rz(-2.9314801) q[0];
rz(-0.70136079) q[1];
sx q[1];
rz(-1.7751834) q[1];
sx q[1];
rz(-1.1022386) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022033545) q[0];
sx q[0];
rz(-2.4071472) q[0];
sx q[0];
rz(2.2950315) q[0];
x q[1];
rz(0.66221018) q[2];
sx q[2];
rz(-1.5689611) q[2];
sx q[2];
rz(-2.3037132) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9241198) q[1];
sx q[1];
rz(-0.77581166) q[1];
sx q[1];
rz(2.6865143) q[1];
x q[2];
rz(-2.5999448) q[3];
sx q[3];
rz(-1.7078425) q[3];
sx q[3];
rz(-2.8769719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.29409757) q[2];
sx q[2];
rz(-1.5832486) q[2];
sx q[2];
rz(1.6451277) q[2];
rz(1.3930813) q[3];
sx q[3];
rz(-2.2025509) q[3];
sx q[3];
rz(2.4699672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74712718) q[0];
sx q[0];
rz(-1.8165996) q[0];
sx q[0];
rz(2.6649244) q[0];
rz(1.3564302) q[1];
sx q[1];
rz(-1.2973659) q[1];
sx q[1];
rz(-0.93723255) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9324324) q[0];
sx q[0];
rz(-1.5346955) q[0];
sx q[0];
rz(-0.88113992) q[0];
rz(-pi) q[1];
x q[1];
rz(1.632948) q[2];
sx q[2];
rz(-0.85374933) q[2];
sx q[2];
rz(0.99822068) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7644219) q[1];
sx q[1];
rz(-2.563919) q[1];
sx q[1];
rz(2.38297) q[1];
rz(2.9415575) q[3];
sx q[3];
rz(-0.45645255) q[3];
sx q[3];
rz(2.4394622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.53200191) q[2];
sx q[2];
rz(-1.8014182) q[2];
sx q[2];
rz(-2.0659921) q[2];
rz(-2.7514451) q[3];
sx q[3];
rz(-1.3107927) q[3];
sx q[3];
rz(2.3384317) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6049062) q[0];
sx q[0];
rz(-2.0669879) q[0];
sx q[0];
rz(-0.54501504) q[0];
rz(-0.99264985) q[1];
sx q[1];
rz(-0.98572171) q[1];
sx q[1];
rz(-1.6880021) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6777991) q[0];
sx q[0];
rz(-2.2132769) q[0];
sx q[0];
rz(0.55046659) q[0];
x q[1];
rz(-3.0849669) q[2];
sx q[2];
rz(-1.4384965) q[2];
sx q[2];
rz(-0.43435979) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7632227) q[1];
sx q[1];
rz(-1.4577796) q[1];
sx q[1];
rz(2.3341353) q[1];
x q[2];
rz(-0.060096459) q[3];
sx q[3];
rz(-1.5240476) q[3];
sx q[3];
rz(-3.0936733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.40889007) q[2];
sx q[2];
rz(-1.3971034) q[2];
sx q[2];
rz(0.2962386) q[2];
rz(-0.57684165) q[3];
sx q[3];
rz(-2.7448765) q[3];
sx q[3];
rz(-1.860994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7883564) q[0];
sx q[0];
rz(-0.78859538) q[0];
sx q[0];
rz(-0.22407918) q[0];
rz(-0.60399404) q[1];
sx q[1];
rz(-0.9915587) q[1];
sx q[1];
rz(1.6557065) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0191658) q[0];
sx q[0];
rz(-0.44394025) q[0];
sx q[0];
rz(-2.5371518) q[0];
x q[1];
rz(0.9233711) q[2];
sx q[2];
rz(-2.8956872) q[2];
sx q[2];
rz(0.78215037) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8454885) q[1];
sx q[1];
rz(-1.316324) q[1];
sx q[1];
rz(2.5135882) q[1];
rz(-pi) q[2];
rz(1.2717683) q[3];
sx q[3];
rz(-2.0439897) q[3];
sx q[3];
rz(2.2132471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.47933444) q[2];
sx q[2];
rz(-1.0666288) q[2];
sx q[2];
rz(2.7545641) q[2];
rz(-2.8582063) q[3];
sx q[3];
rz(-0.75205386) q[3];
sx q[3];
rz(-0.40248218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1166444) q[0];
sx q[0];
rz(-0.270917) q[0];
sx q[0];
rz(1.1234294) q[0];
rz(1.6472389) q[1];
sx q[1];
rz(-1.6554183) q[1];
sx q[1];
rz(2.2656238) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94695222) q[0];
sx q[0];
rz(-2.2966697) q[0];
sx q[0];
rz(-1.2256757) q[0];
rz(-pi) q[1];
rz(-1.0303622) q[2];
sx q[2];
rz(-1.9649995) q[2];
sx q[2];
rz(2.9351007) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3800602) q[1];
sx q[1];
rz(-2.426154) q[1];
sx q[1];
rz(-1.2334633) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.749529) q[3];
sx q[3];
rz(-1.9944085) q[3];
sx q[3];
rz(-2.7152747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7312077) q[2];
sx q[2];
rz(-1.3599675) q[2];
sx q[2];
rz(-3.0044921) q[2];
rz(1.7588663) q[3];
sx q[3];
rz(-0.9226678) q[3];
sx q[3];
rz(1.2285129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41505861) q[0];
sx q[0];
rz(-2.6714323) q[0];
sx q[0];
rz(2.5191125) q[0];
rz(-0.38901916) q[1];
sx q[1];
rz(-0.35094378) q[1];
sx q[1];
rz(2.6019179) q[1];
rz(2.0034267) q[2];
sx q[2];
rz(-2.0180704) q[2];
sx q[2];
rz(1.13712) q[2];
rz(2.9403654) q[3];
sx q[3];
rz(-1.22898) q[3];
sx q[3];
rz(-2.1635319) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
