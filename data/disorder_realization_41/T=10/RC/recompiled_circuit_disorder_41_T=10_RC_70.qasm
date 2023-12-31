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
rz(3.119757) q[0];
rz(-0.39437374) q[1];
sx q[1];
rz(-1.6819277) q[1];
sx q[1];
rz(0.2149166) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7267933) q[0];
sx q[0];
rz(-2.7601295) q[0];
sx q[0];
rz(0.8154072) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6724042) q[2];
sx q[2];
rz(-1.4326296) q[2];
sx q[2];
rz(-0.56433041) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.54675198) q[1];
sx q[1];
rz(-1.1210124) q[1];
sx q[1];
rz(0.91083281) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5932556) q[3];
sx q[3];
rz(-2.3468446) q[3];
sx q[3];
rz(-2.7693975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4102143) q[2];
sx q[2];
rz(-1.4593068) q[2];
sx q[2];
rz(2.5773876) q[2];
rz(-1.365186) q[3];
sx q[3];
rz(-0.44962883) q[3];
sx q[3];
rz(-1.2692497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0974225) q[0];
sx q[0];
rz(-1.2133657) q[0];
sx q[0];
rz(-2.2136097) q[0];
rz(1.9762951) q[1];
sx q[1];
rz(-1.5382643) q[1];
sx q[1];
rz(2.2448418) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1363163) q[0];
sx q[0];
rz(-1.3747842) q[0];
sx q[0];
rz(0.85640237) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4977658) q[2];
sx q[2];
rz(-1.739193) q[2];
sx q[2];
rz(-1.3161236) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0838544) q[1];
sx q[1];
rz(-1.4469622) q[1];
sx q[1];
rz(0.79353516) q[1];
x q[2];
rz(-2.7249343) q[3];
sx q[3];
rz(-2.2778802) q[3];
sx q[3];
rz(1.7435031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8759878) q[2];
sx q[2];
rz(-2.6066055) q[2];
sx q[2];
rz(-2.1014452) q[2];
rz(1.4552207) q[3];
sx q[3];
rz(-1.2929595) q[3];
sx q[3];
rz(-0.35475981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4002157) q[0];
sx q[0];
rz(-2.564036) q[0];
sx q[0];
rz(-1.0282015) q[0];
rz(-2.0630515) q[1];
sx q[1];
rz(-2.5787347) q[1];
sx q[1];
rz(-2.7064586) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9520156) q[0];
sx q[0];
rz(-0.36882419) q[0];
sx q[0];
rz(-0.30216218) q[0];
rz(-1.4235731) q[2];
sx q[2];
rz(-1.8185116) q[2];
sx q[2];
rz(1.148828) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1060461) q[1];
sx q[1];
rz(-1.7340845) q[1];
sx q[1];
rz(1.0825023) q[1];
rz(2.3782303) q[3];
sx q[3];
rz(-1.8244787) q[3];
sx q[3];
rz(-2.7953479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6083287) q[2];
sx q[2];
rz(-1.3003131) q[2];
sx q[2];
rz(2.8386774) q[2];
rz(-1.3251925) q[3];
sx q[3];
rz(-1.9830827) q[3];
sx q[3];
rz(3.0505676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9451697) q[0];
sx q[0];
rz(-1.7315995) q[0];
sx q[0];
rz(2.2241425) q[0];
rz(2.4687185) q[1];
sx q[1];
rz(-2.0560975) q[1];
sx q[1];
rz(-0.26487574) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0565856) q[0];
sx q[0];
rz(-1.6008953) q[0];
sx q[0];
rz(2.4823275) q[0];
rz(-pi) q[1];
x q[1];
rz(1.869309) q[2];
sx q[2];
rz(-1.8372756) q[2];
sx q[2];
rz(2.4556015) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5502364) q[1];
sx q[1];
rz(-1.9758421) q[1];
sx q[1];
rz(2.3637799) q[1];
rz(-2.360965) q[3];
sx q[3];
rz(-0.98140162) q[3];
sx q[3];
rz(1.9115703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.36310568) q[2];
sx q[2];
rz(-2.6532756) q[2];
sx q[2];
rz(-1.5765566) q[2];
rz(-1.0270843) q[3];
sx q[3];
rz(-0.74147195) q[3];
sx q[3];
rz(-2.0402133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89001369) q[0];
sx q[0];
rz(-3.0047834) q[0];
sx q[0];
rz(-0.47873163) q[0];
rz(2.1084673) q[1];
sx q[1];
rz(-0.9712351) q[1];
sx q[1];
rz(-2.1889401) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6081776) q[0];
sx q[0];
rz(-2.1153643) q[0];
sx q[0];
rz(1.4334701) q[0];
rz(2.7258337) q[2];
sx q[2];
rz(-1.8766878) q[2];
sx q[2];
rz(1.3410459) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.18187411) q[1];
sx q[1];
rz(-1.8078783) q[1];
sx q[1];
rz(1.167776) q[1];
rz(1.5289375) q[3];
sx q[3];
rz(-2.0484945) q[3];
sx q[3];
rz(-2.4822513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4218563) q[2];
sx q[2];
rz(-0.36281261) q[2];
sx q[2];
rz(2.6110113) q[2];
rz(1.4060219) q[3];
sx q[3];
rz(-2.0134182) q[3];
sx q[3];
rz(2.3099242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56753165) q[0];
sx q[0];
rz(-1.644779) q[0];
sx q[0];
rz(-1.6249599) q[0];
rz(1.3051055) q[1];
sx q[1];
rz(-1.790698) q[1];
sx q[1];
rz(0.17257246) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1369551) q[0];
sx q[0];
rz(-1.1543659) q[0];
sx q[0];
rz(-3.1266771) q[0];
x q[1];
rz(0.2595915) q[2];
sx q[2];
rz(-1.1278369) q[2];
sx q[2];
rz(1.358658) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1291618) q[1];
sx q[1];
rz(-0.42814246) q[1];
sx q[1];
rz(-1.3151602) q[1];
rz(0.63038007) q[3];
sx q[3];
rz(-1.7582338) q[3];
sx q[3];
rz(0.58296766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.83795786) q[2];
sx q[2];
rz(-1.4540318) q[2];
sx q[2];
rz(1.1266358) q[2];
rz(0.78222328) q[3];
sx q[3];
rz(-1.9061079) q[3];
sx q[3];
rz(-1.3379898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
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
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7717168) q[0];
sx q[0];
rz(-1.4780095) q[0];
sx q[0];
rz(-1.4183527) q[0];
x q[1];
rz(-0.41140326) q[2];
sx q[2];
rz(-1.7559933) q[2];
sx q[2];
rz(1.411737) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7927671) q[1];
sx q[1];
rz(-1.7276689) q[1];
sx q[1];
rz(1.3385593) q[1];
rz(-pi) q[2];
rz(1.5114622) q[3];
sx q[3];
rz(-1.7620798) q[3];
sx q[3];
rz(-0.64185601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0044272) q[2];
sx q[2];
rz(-2.9512773) q[2];
sx q[2];
rz(2.9471617) q[2];
rz(0.91313177) q[3];
sx q[3];
rz(-1.3875995) q[3];
sx q[3];
rz(2.156179) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5450491) q[0];
sx q[0];
rz(-0.61674917) q[0];
sx q[0];
rz(0.066666691) q[0];
rz(-0.32456675) q[1];
sx q[1];
rz(-1.5044183) q[1];
sx q[1];
rz(-2.1527122) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3263771) q[0];
sx q[0];
rz(-1.2412374) q[0];
sx q[0];
rz(1.9944847) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5547995) q[2];
sx q[2];
rz(-1.4774067) q[2];
sx q[2];
rz(2.8430251) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7372072) q[1];
sx q[1];
rz(-2.7215241) q[1];
sx q[1];
rz(-1.7350446) q[1];
rz(-pi) q[2];
rz(0.75108053) q[3];
sx q[3];
rz(-0.6595279) q[3];
sx q[3];
rz(-0.86576033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4618335) q[2];
sx q[2];
rz(-2.2448848) q[2];
sx q[2];
rz(-2.7339593) q[2];
rz(0.76861012) q[3];
sx q[3];
rz(-1.3137484) q[3];
sx q[3];
rz(0.9238981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75893629) q[0];
sx q[0];
rz(-1.3019245) q[0];
sx q[0];
rz(2.5323903) q[0];
rz(0.095104782) q[1];
sx q[1];
rz(-1.2520049) q[1];
sx q[1];
rz(2.2682155) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0926889) q[0];
sx q[0];
rz(-2.9498219) q[0];
sx q[0];
rz(1.9321241) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4344425) q[2];
sx q[2];
rz(-1.4246203) q[2];
sx q[2];
rz(-1.855195) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1904618) q[1];
sx q[1];
rz(-1.8776263) q[1];
sx q[1];
rz(2.5224586) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99336266) q[3];
sx q[3];
rz(-1.1953029) q[3];
sx q[3];
rz(1.4162228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1197027) q[2];
sx q[2];
rz(-2.3535574) q[2];
sx q[2];
rz(-2.4592887) q[2];
rz(0.37426379) q[3];
sx q[3];
rz(-1.572861) q[3];
sx q[3];
rz(-2.9746829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4436214) q[0];
sx q[0];
rz(-1.781783) q[0];
sx q[0];
rz(0.95296729) q[0];
rz(-2.3151746) q[1];
sx q[1];
rz(-0.73917878) q[1];
sx q[1];
rz(1.3964765) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1288426) q[0];
sx q[0];
rz(-1.8930952) q[0];
sx q[0];
rz(1.661257) q[0];
x q[1];
rz(-2.7637485) q[2];
sx q[2];
rz(-2.4477738) q[2];
sx q[2];
rz(-2.5758884) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4954088) q[1];
sx q[1];
rz(-0.42802654) q[1];
sx q[1];
rz(0.82823786) q[1];
rz(0.20874899) q[3];
sx q[3];
rz(-2.8737846) q[3];
sx q[3];
rz(0.97313125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6754127) q[2];
sx q[2];
rz(-2.7853577) q[2];
sx q[2];
rz(-2.9818025) q[2];
rz(-2.8397078) q[3];
sx q[3];
rz(-2.2146137) q[3];
sx q[3];
rz(0.2872428) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044534279) q[0];
sx q[0];
rz(-2.4659768) q[0];
sx q[0];
rz(1.5855047) q[0];
rz(3.0083169) q[1];
sx q[1];
rz(-1.517308) q[1];
sx q[1];
rz(3.0130253) q[1];
rz(-0.940154) q[2];
sx q[2];
rz(-1.7241782) q[2];
sx q[2];
rz(0.45964514) q[2];
rz(-0.53363579) q[3];
sx q[3];
rz(-1.5279557) q[3];
sx q[3];
rz(0.91844311) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
