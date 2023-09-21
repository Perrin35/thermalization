OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.70513201) q[0];
sx q[0];
rz(-2.5897265) q[0];
sx q[0];
rz(-0.021835672) q[0];
rz(-0.39437374) q[1];
sx q[1];
rz(-1.6819277) q[1];
sx q[1];
rz(0.2149166) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41479933) q[0];
sx q[0];
rz(-2.7601295) q[0];
sx q[0];
rz(0.8154072) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29834892) q[2];
sx q[2];
rz(-2.6539408) q[2];
sx q[2];
rz(1.8698486) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.51337459) q[1];
sx q[1];
rz(-2.3623423) q[1];
sx q[1];
rz(0.90374225) q[1];
rz(-pi) q[2];
rz(-2.5932556) q[3];
sx q[3];
rz(-0.79474802) q[3];
sx q[3];
rz(2.7693975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4102143) q[2];
sx q[2];
rz(-1.6822858) q[2];
sx q[2];
rz(2.5773876) q[2];
rz(-1.7764067) q[3];
sx q[3];
rz(-0.44962883) q[3];
sx q[3];
rz(-1.8723429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0974225) q[0];
sx q[0];
rz(-1.2133657) q[0];
sx q[0];
rz(-0.92798293) q[0];
rz(1.1652975) q[1];
sx q[1];
rz(-1.6033283) q[1];
sx q[1];
rz(-0.89675084) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3551536) q[0];
sx q[0];
rz(-2.4053898) q[0];
sx q[0];
rz(1.2765221) q[0];
rz(-pi) q[1];
rz(1.4977658) q[2];
sx q[2];
rz(-1.4023997) q[2];
sx q[2];
rz(-1.3161236) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5337199) q[1];
sx q[1];
rz(-0.80104154) q[1];
sx q[1];
rz(2.968722) q[1];
x q[2];
rz(-2.7249343) q[3];
sx q[3];
rz(-2.2778802) q[3];
sx q[3];
rz(-1.3980896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.26560489) q[2];
sx q[2];
rz(-2.6066055) q[2];
sx q[2];
rz(-2.1014452) q[2];
rz(-1.6863719) q[3];
sx q[3];
rz(-1.8486332) q[3];
sx q[3];
rz(-2.7868328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74137694) q[0];
sx q[0];
rz(-0.57755661) q[0];
sx q[0];
rz(-1.0282015) q[0];
rz(-1.0785412) q[1];
sx q[1];
rz(-0.56285793) q[1];
sx q[1];
rz(-2.7064586) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2745278) q[0];
sx q[0];
rz(-1.2194249) q[0];
sx q[0];
rz(-1.4562796) q[0];
rz(-0.52559678) q[2];
sx q[2];
rz(-0.28738775) q[2];
sx q[2];
rz(-1.692786) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5907026) q[1];
sx q[1];
rz(-2.0520376) q[1];
sx q[1];
rz(-0.18443702) q[1];
x q[2];
rz(-2.3782303) q[3];
sx q[3];
rz(-1.3171139) q[3];
sx q[3];
rz(0.34624472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6083287) q[2];
sx q[2];
rz(-1.8412795) q[2];
sx q[2];
rz(-2.8386774) q[2];
rz(1.8164002) q[3];
sx q[3];
rz(-1.9830827) q[3];
sx q[3];
rz(3.0505676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19642297) q[0];
sx q[0];
rz(-1.7315995) q[0];
sx q[0];
rz(2.2241425) q[0];
rz(-2.4687185) q[1];
sx q[1];
rz(-1.0854951) q[1];
sx q[1];
rz(-0.26487574) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0850071) q[0];
sx q[0];
rz(-1.5406973) q[0];
sx q[0];
rz(-2.4823275) q[0];
x q[1];
rz(1.869309) q[2];
sx q[2];
rz(-1.8372756) q[2];
sx q[2];
rz(2.4556015) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5323822) q[1];
sx q[1];
rz(-0.86984837) q[1];
sx q[1];
rz(-2.1125395) q[1];
rz(0.81569205) q[3];
sx q[3];
rz(-0.94592735) q[3];
sx q[3];
rz(-0.16251414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.36310568) q[2];
sx q[2];
rz(-2.6532756) q[2];
sx q[2];
rz(1.5650361) q[2];
rz(-1.0270843) q[3];
sx q[3];
rz(-0.74147195) q[3];
sx q[3];
rz(-2.0402133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.251579) q[0];
sx q[0];
rz(-3.0047834) q[0];
sx q[0];
rz(2.662861) q[0];
rz(2.1084673) q[1];
sx q[1];
rz(-2.1703576) q[1];
sx q[1];
rz(-0.95265257) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.032741) q[0];
sx q[0];
rz(-1.6881588) q[0];
sx q[0];
rz(-0.54876859) q[0];
x q[1];
rz(1.9031992) q[2];
sx q[2];
rz(-1.1754416) q[2];
sx q[2];
rz(-2.7796641) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2561803) q[1];
sx q[1];
rz(-0.46426877) q[1];
sx q[1];
rz(-2.1229565) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5289375) q[3];
sx q[3];
rz(-2.0484945) q[3];
sx q[3];
rz(2.4822513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4218563) q[2];
sx q[2];
rz(-2.77878) q[2];
sx q[2];
rz(2.6110113) q[2];
rz(1.7355708) q[3];
sx q[3];
rz(-1.1281745) q[3];
sx q[3];
rz(-0.83166844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.56753165) q[0];
sx q[0];
rz(-1.644779) q[0];
sx q[0];
rz(1.6249599) q[0];
rz(-1.8364871) q[1];
sx q[1];
rz(-1.790698) q[1];
sx q[1];
rz(-2.9690202) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0046376) q[0];
sx q[0];
rz(-1.9872268) q[0];
sx q[0];
rz(0.01491551) q[0];
x q[1];
rz(2.0270945) q[2];
sx q[2];
rz(-1.8048394) q[2];
sx q[2];
rz(2.8161088) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1291618) q[1];
sx q[1];
rz(-0.42814246) q[1];
sx q[1];
rz(-1.3151602) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8014088) q[3];
sx q[3];
rz(-0.95313822) q[3];
sx q[3];
rz(-0.85268439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.83795786) q[2];
sx q[2];
rz(-1.4540318) q[2];
sx q[2];
rz(-1.1266358) q[2];
rz(-2.3593694) q[3];
sx q[3];
rz(-1.9061079) q[3];
sx q[3];
rz(-1.3379898) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28850266) q[0];
sx q[0];
rz(-2.8350916) q[0];
sx q[0];
rz(-2.4801168) q[0];
rz(2.181197) q[1];
sx q[1];
rz(-1.7405225) q[1];
sx q[1];
rz(-0.75659928) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21515439) q[0];
sx q[0];
rz(-1.7225791) q[0];
sx q[0];
rz(-3.0477235) q[0];
rz(-2.7034764) q[2];
sx q[2];
rz(-0.44898673) q[2];
sx q[2];
rz(-2.9012836) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7927671) q[1];
sx q[1];
rz(-1.4139237) q[1];
sx q[1];
rz(-1.8030333) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5114622) q[3];
sx q[3];
rz(-1.7620798) q[3];
sx q[3];
rz(-2.4997366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1371655) q[2];
sx q[2];
rz(-2.9512773) q[2];
sx q[2];
rz(2.9471617) q[2];
rz(2.2284609) q[3];
sx q[3];
rz(-1.7539932) q[3];
sx q[3];
rz(2.156179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59654355) q[0];
sx q[0];
rz(-2.5248435) q[0];
sx q[0];
rz(-0.066666691) q[0];
rz(-2.8170259) q[1];
sx q[1];
rz(-1.5044183) q[1];
sx q[1];
rz(2.1527122) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37799997) q[0];
sx q[0];
rz(-2.6110296) q[0];
sx q[0];
rz(-2.2647122) q[0];
rz(-pi) q[1];
rz(-1.5547995) q[2];
sx q[2];
rz(-1.664186) q[2];
sx q[2];
rz(-2.8430251) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0162184) q[1];
sx q[1];
rz(-1.6375293) q[1];
sx q[1];
rz(1.9858422) q[1];
rz(-pi) q[2];
rz(2.3905121) q[3];
sx q[3];
rz(-0.6595279) q[3];
sx q[3];
rz(-2.2758323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6797592) q[2];
sx q[2];
rz(-2.2448848) q[2];
sx q[2];
rz(-2.7339593) q[2];
rz(-0.76861012) q[3];
sx q[3];
rz(-1.3137484) q[3];
sx q[3];
rz(2.2176946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3826564) q[0];
sx q[0];
rz(-1.8396682) q[0];
sx q[0];
rz(-0.60920238) q[0];
rz(0.095104782) q[1];
sx q[1];
rz(-1.8895878) q[1];
sx q[1];
rz(0.87337714) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4163923) q[0];
sx q[0];
rz(-1.3915477) q[0];
sx q[0];
rz(-0.068530131) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74553211) q[2];
sx q[2];
rz(-0.19956707) q[2];
sx q[2];
rz(2.0419288) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5492591) q[1];
sx q[1];
rz(-0.98456406) q[1];
sx q[1];
rz(1.9418282) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94536762) q[3];
sx q[3];
rz(-2.4646467) q[3];
sx q[3];
rz(-0.35811801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0218899) q[2];
sx q[2];
rz(-0.78803524) q[2];
sx q[2];
rz(0.68230391) q[2];
rz(-0.37426379) q[3];
sx q[3];
rz(-1.5687317) q[3];
sx q[3];
rz(0.16690978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4436214) q[0];
sx q[0];
rz(-1.781783) q[0];
sx q[0];
rz(2.1886254) q[0];
rz(2.3151746) q[1];
sx q[1];
rz(-0.73917878) q[1];
sx q[1];
rz(1.7451161) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1288426) q[0];
sx q[0];
rz(-1.8930952) q[0];
sx q[0];
rz(1.661257) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65812494) q[2];
sx q[2];
rz(-1.8089559) q[2];
sx q[2];
rz(2.4326774) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.64618387) q[1];
sx q[1];
rz(-2.7135661) q[1];
sx q[1];
rz(0.82823786) q[1];
rz(-pi) q[2];
rz(2.8793328) q[3];
sx q[3];
rz(-1.5159303) q[3];
sx q[3];
rz(-0.39615397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.46618) q[2];
sx q[2];
rz(-2.7853577) q[2];
sx q[2];
rz(-2.9818025) q[2];
rz(-2.8397078) q[3];
sx q[3];
rz(-0.92697898) q[3];
sx q[3];
rz(-0.2872428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0970584) q[0];
sx q[0];
rz(-2.4659768) q[0];
sx q[0];
rz(1.5855047) q[0];
rz(3.0083169) q[1];
sx q[1];
rz(-1.517308) q[1];
sx q[1];
rz(3.0130253) q[1];
rz(-1.3143905) q[2];
sx q[2];
rz(-2.495043) q[2];
sx q[2];
rz(-1.3174353) q[2];
rz(0.084074323) q[3];
sx q[3];
rz(-2.6064059) q[3];
sx q[3];
rz(2.4168766) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];