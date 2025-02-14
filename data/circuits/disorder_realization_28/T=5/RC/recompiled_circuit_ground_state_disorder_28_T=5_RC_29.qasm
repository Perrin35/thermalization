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
rz(-0.063152753) q[0];
rz(1.8237279) q[1];
sx q[1];
rz(-0.39931077) q[1];
sx q[1];
rz(2.3656486) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5781893) q[0];
sx q[0];
rz(-0.052980352) q[0];
sx q[0];
rz(2.6107759) q[0];
rz(2.2065483) q[2];
sx q[2];
rz(-1.9823977) q[2];
sx q[2];
rz(1.7807478) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9358238) q[1];
sx q[1];
rz(-2.4112114) q[1];
sx q[1];
rz(-1.8380866) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3315866) q[3];
sx q[3];
rz(-1.080371) q[3];
sx q[3];
rz(-1.8338721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6905288) q[2];
sx q[2];
rz(-1.3961926) q[2];
sx q[2];
rz(2.6032791) q[2];
rz(-2.8484143) q[3];
sx q[3];
rz(-1.5978975) q[3];
sx q[3];
rz(-1.8831016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20741589) q[0];
sx q[0];
rz(-2.2947312) q[0];
sx q[0];
rz(2.9127981) q[0];
rz(3.0711807) q[1];
sx q[1];
rz(-1.8682624) q[1];
sx q[1];
rz(1.061903) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.966514) q[0];
sx q[0];
rz(-1.0565149) q[0];
sx q[0];
rz(-1.8755336) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0309762) q[2];
sx q[2];
rz(-2.0041582) q[2];
sx q[2];
rz(-3.1040807) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3913566) q[1];
sx q[1];
rz(-1.8614381) q[1];
sx q[1];
rz(2.802317) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44159378) q[3];
sx q[3];
rz(-2.674235) q[3];
sx q[3];
rz(2.2668214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9262907) q[2];
sx q[2];
rz(-0.44507241) q[2];
sx q[2];
rz(-0.11623795) q[2];
rz(-0.91533533) q[3];
sx q[3];
rz(-1.4696308) q[3];
sx q[3];
rz(0.62469971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-1.4794469) q[0];
sx q[0];
rz(-1.5597458) q[0];
sx q[0];
rz(-0.33552718) q[0];
rz(-1.7299995) q[1];
sx q[1];
rz(-2.2970707) q[1];
sx q[1];
rz(1.9427049) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20366448) q[0];
sx q[0];
rz(-2.0988582) q[0];
sx q[0];
rz(-2.8115134) q[0];
rz(-2.6012318) q[2];
sx q[2];
rz(-0.66695222) q[2];
sx q[2];
rz(1.4948542) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80497031) q[1];
sx q[1];
rz(-1.687007) q[1];
sx q[1];
rz(2.3776965) q[1];
rz(-pi) q[2];
rz(2.8923404) q[3];
sx q[3];
rz(-1.0947168) q[3];
sx q[3];
rz(2.2152546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.18232839) q[2];
sx q[2];
rz(-2.3458643) q[2];
sx q[2];
rz(1.1129334) q[2];
rz(-0.5111323) q[3];
sx q[3];
rz(-1.7685578) q[3];
sx q[3];
rz(2.6326877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48753259) q[0];
sx q[0];
rz(-0.075981058) q[0];
sx q[0];
rz(0.37387601) q[0];
rz(1.182425) q[1];
sx q[1];
rz(-1.9805209) q[1];
sx q[1];
rz(-2.9808796) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.041695853) q[0];
sx q[0];
rz(-1.6536435) q[0];
sx q[0];
rz(3.003398) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4771551) q[2];
sx q[2];
rz(-2.2291846) q[2];
sx q[2];
rz(0.13641549) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4121033) q[1];
sx q[1];
rz(-2.0827052) q[1];
sx q[1];
rz(0.44513925) q[1];
rz(-0.31779685) q[3];
sx q[3];
rz(-1.2850015) q[3];
sx q[3];
rz(0.60493776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5339727) q[2];
sx q[2];
rz(-1.1090358) q[2];
sx q[2];
rz(-2.4268761) q[2];
rz(2.886582) q[3];
sx q[3];
rz(-1.3304354) q[3];
sx q[3];
rz(2.3431006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(3.0890927) q[0];
sx q[0];
rz(-0.30464259) q[0];
sx q[0];
rz(-2.3783045) q[0];
rz(-0.25453645) q[1];
sx q[1];
rz(-1.7691879) q[1];
sx q[1];
rz(-1.6932142) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096222046) q[0];
sx q[0];
rz(-1.9217446) q[0];
sx q[0];
rz(-2.778102) q[0];
rz(1.4035475) q[2];
sx q[2];
rz(-2.3838701) q[2];
sx q[2];
rz(3.0191985) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9336945) q[1];
sx q[1];
rz(-0.95732821) q[1];
sx q[1];
rz(1.3529569) q[1];
x q[2];
rz(-0.71131018) q[3];
sx q[3];
rz(-1.2625202) q[3];
sx q[3];
rz(0.084567794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.86108565) q[2];
sx q[2];
rz(-2.0528767) q[2];
sx q[2];
rz(-1.0293993) q[2];
rz(-1.5329125) q[3];
sx q[3];
rz(-0.89919388) q[3];
sx q[3];
rz(2.5077584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7480046) q[0];
sx q[0];
rz(-2.2279976) q[0];
sx q[0];
rz(0.21011259) q[0];
rz(2.4402319) q[1];
sx q[1];
rz(-1.3664093) q[1];
sx q[1];
rz(-2.0393541) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85066807) q[0];
sx q[0];
rz(-2.0966665) q[0];
sx q[0];
rz(-2.6024502) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66221018) q[2];
sx q[2];
rz(-1.5689611) q[2];
sx q[2];
rz(-2.3037132) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6893856) q[1];
sx q[1];
rz(-1.2579134) q[1];
sx q[1];
rz(-2.41928) q[1];
rz(0.54164782) q[3];
sx q[3];
rz(-1.7078425) q[3];
sx q[3];
rz(0.26462072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.29409757) q[2];
sx q[2];
rz(-1.558344) q[2];
sx q[2];
rz(1.6451277) q[2];
rz(1.3930813) q[3];
sx q[3];
rz(-0.93904177) q[3];
sx q[3];
rz(0.67162544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74712718) q[0];
sx q[0];
rz(-1.8165996) q[0];
sx q[0];
rz(-2.6649244) q[0];
rz(-1.3564302) q[1];
sx q[1];
rz(-1.2973659) q[1];
sx q[1];
rz(-2.2043601) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3318769) q[0];
sx q[0];
rz(-2.2599155) q[0];
sx q[0];
rz(-3.0948113) q[0];
rz(-pi) q[1];
x q[1];
rz(1.632948) q[2];
sx q[2];
rz(-2.2878433) q[2];
sx q[2];
rz(2.143372) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.46980935) q[1];
sx q[1];
rz(-1.1632763) q[1];
sx q[1];
rz(1.1492561) q[1];
rz(-2.6930599) q[3];
sx q[3];
rz(-1.6584907) q[3];
sx q[3];
rz(-2.0929071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6095907) q[2];
sx q[2];
rz(-1.8014182) q[2];
sx q[2];
rz(2.0659921) q[2];
rz(0.39014751) q[3];
sx q[3];
rz(-1.3107927) q[3];
sx q[3];
rz(-0.80316097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53668642) q[0];
sx q[0];
rz(-2.0669879) q[0];
sx q[0];
rz(0.54501504) q[0];
rz(-2.1489428) q[1];
sx q[1];
rz(-2.1558709) q[1];
sx q[1];
rz(-1.6880021) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75460583) q[0];
sx q[0];
rz(-2.0029066) q[0];
sx q[0];
rz(2.291392) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7033061) q[2];
sx q[2];
rz(-1.514666) q[2];
sx q[2];
rz(-2.0126337) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8416748) q[1];
sx q[1];
rz(-2.3280488) q[1];
sx q[1];
rz(0.15581375) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66176118) q[3];
sx q[3];
rz(-0.076120928) q[3];
sx q[3];
rz(0.95835987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7327026) q[2];
sx q[2];
rz(-1.7444892) q[2];
sx q[2];
rz(2.8453541) q[2];
rz(2.564751) q[3];
sx q[3];
rz(-0.39671612) q[3];
sx q[3];
rz(1.860994) q[3];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7883564) q[0];
sx q[0];
rz(-0.78859538) q[0];
sx q[0];
rz(2.9175135) q[0];
rz(2.5375986) q[1];
sx q[1];
rz(-0.9915587) q[1];
sx q[1];
rz(1.6557065) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1224269) q[0];
sx q[0];
rz(-0.44394025) q[0];
sx q[0];
rz(-2.5371518) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2182216) q[2];
sx q[2];
rz(-0.24590547) q[2];
sx q[2];
rz(-0.78215037) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.60843492) q[1];
sx q[1];
rz(-0.67109334) q[1];
sx q[1];
rz(0.41678269) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6497439) q[3];
sx q[3];
rz(-1.8361194) q[3];
sx q[3];
rz(-0.50287528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6622582) q[2];
sx q[2];
rz(-1.0666288) q[2];
sx q[2];
rz(0.38702854) q[2];
rz(-0.28338638) q[3];
sx q[3];
rz(-2.3895388) q[3];
sx q[3];
rz(2.7391105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0249483) q[0];
sx q[0];
rz(-2.8706757) q[0];
sx q[0];
rz(-1.1234294) q[0];
rz(1.6472389) q[1];
sx q[1];
rz(-1.4861743) q[1];
sx q[1];
rz(-2.2656238) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94695222) q[0];
sx q[0];
rz(-2.2966697) q[0];
sx q[0];
rz(1.915917) q[0];
rz(-pi) q[1];
x q[1];
rz(0.89089762) q[2];
sx q[2];
rz(-2.4844137) q[2];
sx q[2];
rz(-0.79492043) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3265423) q[1];
sx q[1];
rz(-0.90335323) q[1];
sx q[1];
rz(-0.28006458) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7660329) q[3];
sx q[3];
rz(-0.45766214) q[3];
sx q[3];
rz(0.012227146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4103849) q[2];
sx q[2];
rz(-1.7816252) q[2];
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
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41505861) q[0];
sx q[0];
rz(-0.4701604) q[0];
sx q[0];
rz(-0.62248019) q[0];
rz(-0.38901916) q[1];
sx q[1];
rz(-0.35094378) q[1];
sx q[1];
rz(2.6019179) q[1];
rz(1.1381659) q[2];
sx q[2];
rz(-1.1235222) q[2];
sx q[2];
rz(-2.0044727) q[2];
rz(-1.9191011) q[3];
sx q[3];
rz(-1.7602362) q[3];
sx q[3];
rz(2.6171274) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
