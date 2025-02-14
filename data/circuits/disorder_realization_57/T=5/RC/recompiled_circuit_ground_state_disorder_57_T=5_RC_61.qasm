OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.39855555) q[0];
sx q[0];
rz(-0.86657137) q[0];
sx q[0];
rz(0.33696365) q[0];
rz(0.26951867) q[1];
sx q[1];
rz(-2.1989006) q[1];
sx q[1];
rz(1.2595133) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90010479) q[0];
sx q[0];
rz(-0.83118992) q[0];
sx q[0];
rz(0.8054975) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7808229) q[2];
sx q[2];
rz(-0.84473306) q[2];
sx q[2];
rz(3.0371916) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9971157) q[1];
sx q[1];
rz(-1.0273178) q[1];
sx q[1];
rz(-1.687084) q[1];
rz(-0.99997493) q[3];
sx q[3];
rz(-0.68203841) q[3];
sx q[3];
rz(1.2004971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.17244615) q[2];
sx q[2];
rz(-1.5772051) q[2];
sx q[2];
rz(-1.9077612) q[2];
rz(-1.6110169) q[3];
sx q[3];
rz(-1.4065892) q[3];
sx q[3];
rz(-0.56510258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81561404) q[0];
sx q[0];
rz(-2.1765206) q[0];
sx q[0];
rz(-0.24398971) q[0];
rz(-1.1583534) q[1];
sx q[1];
rz(-2.127425) q[1];
sx q[1];
rz(2.6932531) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2650484) q[0];
sx q[0];
rz(-2.6175313) q[0];
sx q[0];
rz(2.7830809) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82996394) q[2];
sx q[2];
rz(-1.1051851) q[2];
sx q[2];
rz(-0.59620406) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.92384752) q[1];
sx q[1];
rz(-0.73647803) q[1];
sx q[1];
rz(-0.86935161) q[1];
rz(-0.48359343) q[3];
sx q[3];
rz(-1.075287) q[3];
sx q[3];
rz(-2.7413975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9490732) q[2];
sx q[2];
rz(-0.039704617) q[2];
sx q[2];
rz(2.330759) q[2];
rz(-3.132498) q[3];
sx q[3];
rz(-1.5261212) q[3];
sx q[3];
rz(-0.31753376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4119165) q[0];
sx q[0];
rz(-1.3443953) q[0];
sx q[0];
rz(-0.94773951) q[0];
rz(0.89836994) q[1];
sx q[1];
rz(-0.71304524) q[1];
sx q[1];
rz(0.20634849) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8074051) q[0];
sx q[0];
rz(-1.5631521) q[0];
sx q[0];
rz(-0.81595465) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1694109) q[2];
sx q[2];
rz(-2.3222409) q[2];
sx q[2];
rz(-1.8470905) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0230081) q[1];
sx q[1];
rz(-2.6381341) q[1];
sx q[1];
rz(-2.0082824) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6177267) q[3];
sx q[3];
rz(-1.2447671) q[3];
sx q[3];
rz(-2.8992931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4121033) q[2];
sx q[2];
rz(-0.99200839) q[2];
sx q[2];
rz(-2.1573529) q[2];
rz(-0.49736831) q[3];
sx q[3];
rz(-1.2593185) q[3];
sx q[3];
rz(2.3965912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7797778) q[0];
sx q[0];
rz(-3.0497157) q[0];
sx q[0];
rz(0.20633695) q[0];
rz(-1.4641209) q[1];
sx q[1];
rz(-2.3787777) q[1];
sx q[1];
rz(2.5154571) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0447606) q[0];
sx q[0];
rz(-0.87327582) q[0];
sx q[0];
rz(0.4647073) q[0];
rz(-pi) q[1];
rz(3.1244833) q[2];
sx q[2];
rz(-0.10280156) q[2];
sx q[2];
rz(-2.4633165) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5474769) q[1];
sx q[1];
rz(-2.0491338) q[1];
sx q[1];
rz(-2.5376471) q[1];
rz(3.042661) q[3];
sx q[3];
rz(-1.9557992) q[3];
sx q[3];
rz(-1.0855293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3760486) q[2];
sx q[2];
rz(-0.77771336) q[2];
sx q[2];
rz(-2.1518339) q[2];
rz(-1.0342213) q[3];
sx q[3];
rz(-1.9725622) q[3];
sx q[3];
rz(2.6160252) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.702221) q[0];
sx q[0];
rz(-1.419743) q[0];
sx q[0];
rz(1.1628344) q[0];
rz(-2.974466) q[1];
sx q[1];
rz(-1.6163328) q[1];
sx q[1];
rz(-0.67536813) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5655095) q[0];
sx q[0];
rz(-1.352549) q[0];
sx q[0];
rz(-0.670114) q[0];
x q[1];
rz(-2.295601) q[2];
sx q[2];
rz(-0.48946871) q[2];
sx q[2];
rz(2.6014525) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8307504) q[1];
sx q[1];
rz(-1.1427587) q[1];
sx q[1];
rz(0.38106783) q[1];
rz(-2.1618103) q[3];
sx q[3];
rz(-1.3434935) q[3];
sx q[3];
rz(1.0349555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2165788) q[2];
sx q[2];
rz(-1.4224956) q[2];
sx q[2];
rz(-0.50755802) q[2];
rz(0.01928586) q[3];
sx q[3];
rz(-0.79500335) q[3];
sx q[3];
rz(1.9222586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60292202) q[0];
sx q[0];
rz(-0.62748533) q[0];
sx q[0];
rz(2.4399309) q[0];
rz(-2.498846) q[1];
sx q[1];
rz(-1.9822491) q[1];
sx q[1];
rz(-1.8205551) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76687251) q[0];
sx q[0];
rz(-1.5581616) q[0];
sx q[0];
rz(-2.1983474) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3568383) q[2];
sx q[2];
rz(-2.2388487) q[2];
sx q[2];
rz(0.84260637) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.031547) q[1];
sx q[1];
rz(-2.4737691) q[1];
sx q[1];
rz(-0.948443) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.28863971) q[3];
sx q[3];
rz(-0.83191365) q[3];
sx q[3];
rz(-2.6872241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8464437) q[2];
sx q[2];
rz(-2.6602793) q[2];
sx q[2];
rz(-0.63977891) q[2];
rz(-2.4844737) q[3];
sx q[3];
rz(-1.492447) q[3];
sx q[3];
rz(-0.28321987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.101864) q[0];
sx q[0];
rz(-1.2861847) q[0];
sx q[0];
rz(-0.78045994) q[0];
rz(1.9038433) q[1];
sx q[1];
rz(-1.2982118) q[1];
sx q[1];
rz(-0.16403988) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0152978) q[0];
sx q[0];
rz(-1.3336714) q[0];
sx q[0];
rz(2.1046647) q[0];
rz(-pi) q[1];
rz(-2.3140644) q[2];
sx q[2];
rz(-2.568576) q[2];
sx q[2];
rz(0.94486754) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.32719993) q[1];
sx q[1];
rz(-0.76451028) q[1];
sx q[1];
rz(-2.3491377) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0290506) q[3];
sx q[3];
rz(-2.3083105) q[3];
sx q[3];
rz(2.737711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.26874545) q[2];
sx q[2];
rz(-1.9445323) q[2];
sx q[2];
rz(0.8806814) q[2];
rz(-1.667048) q[3];
sx q[3];
rz(-1.0677974) q[3];
sx q[3];
rz(1.4145981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47377652) q[0];
sx q[0];
rz(-0.18874636) q[0];
sx q[0];
rz(-2.012398) q[0];
rz(0.95343626) q[1];
sx q[1];
rz(-2.2013181) q[1];
sx q[1];
rz(0.58852351) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55710965) q[0];
sx q[0];
rz(-2.5059688) q[0];
sx q[0];
rz(2.0899713) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2382201) q[2];
sx q[2];
rz(-1.3540478) q[2];
sx q[2];
rz(-0.34453604) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5626763) q[1];
sx q[1];
rz(-0.92797083) q[1];
sx q[1];
rz(-1.972727) q[1];
x q[2];
rz(1.7950012) q[3];
sx q[3];
rz(-1.5191169) q[3];
sx q[3];
rz(1.3286852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9863161) q[2];
sx q[2];
rz(-1.0078112) q[2];
sx q[2];
rz(-2.0493832) q[2];
rz(-2.3327667) q[3];
sx q[3];
rz(-1.0816962) q[3];
sx q[3];
rz(-1.4441747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5539733) q[0];
sx q[0];
rz(-2.0646136) q[0];
sx q[0];
rz(-0.67163604) q[0];
rz(0.49297586) q[1];
sx q[1];
rz(-2.2608345) q[1];
sx q[1];
rz(-1.2395476) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9705479) q[0];
sx q[0];
rz(-1.8080304) q[0];
sx q[0];
rz(3.1389272) q[0];
x q[1];
rz(3.0918248) q[2];
sx q[2];
rz(-2.9441212) q[2];
sx q[2];
rz(1.2021499) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0860492) q[1];
sx q[1];
rz(-2.1430839) q[1];
sx q[1];
rz(-2.0949165) q[1];
rz(2.3726241) q[3];
sx q[3];
rz(-1.4533224) q[3];
sx q[3];
rz(-0.75107915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0121743) q[2];
sx q[2];
rz(-1.47379) q[2];
sx q[2];
rz(-2.5950281) q[2];
rz(-0.92285815) q[3];
sx q[3];
rz(-2.512629) q[3];
sx q[3];
rz(-1.9950689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4839812) q[0];
sx q[0];
rz(-0.46009362) q[0];
sx q[0];
rz(2.7125603) q[0];
rz(-1.2826762) q[1];
sx q[1];
rz(-1.3653283) q[1];
sx q[1];
rz(0.081258953) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0711533) q[0];
sx q[0];
rz(-0.18658802) q[0];
sx q[0];
rz(-0.43796118) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0308915) q[2];
sx q[2];
rz(-1.541637) q[2];
sx q[2];
rz(-2.3612446) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.024549896) q[1];
sx q[1];
rz(-0.83210841) q[1];
sx q[1];
rz(2.1348272) q[1];
x q[2];
rz(-0.83286442) q[3];
sx q[3];
rz(-0.64117764) q[3];
sx q[3];
rz(2.8744551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.405534) q[2];
sx q[2];
rz(-2.2106876) q[2];
sx q[2];
rz(-1.1116568) q[2];
rz(1.9723655) q[3];
sx q[3];
rz(-0.71318316) q[3];
sx q[3];
rz(0.14298239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1375167) q[0];
sx q[0];
rz(-2.5181894) q[0];
sx q[0];
rz(-0.85859437) q[0];
rz(-0.89523347) q[1];
sx q[1];
rz(-1.8276855) q[1];
sx q[1];
rz(-3.102416) q[1];
rz(1.5262114) q[2];
sx q[2];
rz(-2.8845061) q[2];
sx q[2];
rz(1.3731879) q[2];
rz(-1.7029892) q[3];
sx q[3];
rz(-0.93441603) q[3];
sx q[3];
rz(-0.84032755) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
