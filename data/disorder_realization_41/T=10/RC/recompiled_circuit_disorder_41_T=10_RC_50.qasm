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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2673187) q[0];
sx q[0];
rz(-1.8288757) q[0];
sx q[0];
rz(-1.8549071) q[0];
rz(-pi) q[1];
rz(-1.7254513) q[2];
sx q[2];
rz(-2.0351595) q[2];
sx q[2];
rz(2.2048339) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.79214) q[1];
sx q[1];
rz(-2.1556902) q[1];
sx q[1];
rz(2.5930415) q[1];
rz(1.0825726) q[3];
sx q[3];
rz(-2.2256652) q[3];
sx q[3];
rz(-2.0522576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4102143) q[2];
sx q[2];
rz(-1.6822858) q[2];
sx q[2];
rz(0.56420502) q[2];
rz(1.7764067) q[3];
sx q[3];
rz(-2.6919638) q[3];
sx q[3];
rz(1.2692497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-1.0974225) q[0];
sx q[0];
rz(-1.928227) q[0];
sx q[0];
rz(0.92798293) q[0];
rz(1.9762951) q[1];
sx q[1];
rz(-1.6033283) q[1];
sx q[1];
rz(-2.2448418) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78643909) q[0];
sx q[0];
rz(-0.73620287) q[0];
sx q[0];
rz(-1.2765221) q[0];
rz(-pi) q[1];
rz(0.40541655) q[2];
sx q[2];
rz(-2.9581796) q[2];
sx q[2];
rz(0.90454067) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7795418) q[1];
sx q[1];
rz(-2.3565787) q[1];
sx q[1];
rz(1.3951468) q[1];
rz(2.0131301) q[3];
sx q[3];
rz(-2.3395174) q[3];
sx q[3];
rz(-2.341552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8759878) q[2];
sx q[2];
rz(-0.53498712) q[2];
sx q[2];
rz(1.0401475) q[2];
rz(1.6863719) q[3];
sx q[3];
rz(-1.2929595) q[3];
sx q[3];
rz(-2.7868328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74137694) q[0];
sx q[0];
rz(-0.57755661) q[0];
sx q[0];
rz(1.0282015) q[0];
rz(1.0785412) q[1];
sx q[1];
rz(-0.56285793) q[1];
sx q[1];
rz(-0.43513402) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66416392) q[0];
sx q[0];
rz(-1.6782883) q[0];
sx q[0];
rz(0.35350032) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6159959) q[2];
sx q[2];
rz(-2.8542049) q[2];
sx q[2];
rz(1.4488066) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1060461) q[1];
sx q[1];
rz(-1.4075081) q[1];
sx q[1];
rz(-2.0590904) q[1];
rz(2.3782303) q[3];
sx q[3];
rz(-1.3171139) q[3];
sx q[3];
rz(2.7953479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.53326398) q[2];
sx q[2];
rz(-1.3003131) q[2];
sx q[2];
rz(-0.30291525) q[2];
rz(1.3251925) q[3];
sx q[3];
rz(-1.9830827) q[3];
sx q[3];
rz(0.091025092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
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
rz(0.67287412) q[1];
sx q[1];
rz(-2.0560975) q[1];
sx q[1];
rz(-2.8767169) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0565856) q[0];
sx q[0];
rz(-1.5406973) q[0];
sx q[0];
rz(-0.65926512) q[0];
rz(-2.3189544) q[2];
sx q[2];
rz(-2.7441141) q[2];
sx q[2];
rz(2.9646404) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5502364) q[1];
sx q[1];
rz(-1.9758421) q[1];
sx q[1];
rz(0.77781271) q[1];
x q[2];
rz(-0.75986741) q[3];
sx q[3];
rz(-2.2025975) q[3];
sx q[3];
rz(0.85216537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.36310568) q[2];
sx q[2];
rz(-0.48831707) q[2];
sx q[2];
rz(-1.5650361) q[2];
rz(2.1145084) q[3];
sx q[3];
rz(-0.74147195) q[3];
sx q[3];
rz(1.1013793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(-2.1703576) q[1];
sx q[1];
rz(-0.95265257) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53341502) q[0];
sx q[0];
rz(-2.1153643) q[0];
sx q[0];
rz(-1.4334701) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9031992) q[2];
sx q[2];
rz(-1.9661511) q[2];
sx q[2];
rz(2.7796641) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.488727) q[1];
sx q[1];
rz(-1.9619202) q[1];
sx q[1];
rz(-2.8847242) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6126552) q[3];
sx q[3];
rz(-1.0930982) q[3];
sx q[3];
rz(0.65934138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4218563) q[2];
sx q[2];
rz(-2.77878) q[2];
sx q[2];
rz(-0.53058132) q[2];
rz(1.4060219) q[3];
sx q[3];
rz(-2.0134182) q[3];
sx q[3];
rz(2.3099242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56753165) q[0];
sx q[0];
rz(-1.4968137) q[0];
sx q[0];
rz(1.5166327) q[0];
rz(1.3051055) q[1];
sx q[1];
rz(-1.790698) q[1];
sx q[1];
rz(0.17257246) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5694002) q[0];
sx q[0];
rz(-1.5571556) q[0];
sx q[0];
rz(-1.9872679) q[0];
rz(-pi) q[1];
rz(-2.8820011) q[2];
sx q[2];
rz(-1.1278369) q[2];
sx q[2];
rz(-1.7829347) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.012430819) q[1];
sx q[1];
rz(-2.7134502) q[1];
sx q[1];
rz(1.8264324) q[1];
rz(-pi) q[2];
rz(-0.31130143) q[3];
sx q[3];
rz(-2.4875896) q[3];
sx q[3];
rz(1.2378539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.83795786) q[2];
sx q[2];
rz(-1.4540318) q[2];
sx q[2];
rz(2.0149569) q[2];
rz(2.3593694) q[3];
sx q[3];
rz(-1.9061079) q[3];
sx q[3];
rz(1.3379898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.85309) q[0];
sx q[0];
rz(-0.30650109) q[0];
sx q[0];
rz(-2.4801168) q[0];
rz(-2.181197) q[1];
sx q[1];
rz(-1.4010701) q[1];
sx q[1];
rz(2.3849934) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3698759) q[0];
sx q[0];
rz(-1.4780095) q[0];
sx q[0];
rz(-1.4183527) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7034764) q[2];
sx q[2];
rz(-0.44898673) q[2];
sx q[2];
rz(2.9012836) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3488256) q[1];
sx q[1];
rz(-1.7276689) q[1];
sx q[1];
rz(1.8030333) q[1];
x q[2];
rz(0.29715085) q[3];
sx q[3];
rz(-0.20016709) q[3];
sx q[3];
rz(-2.196892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1371655) q[2];
sx q[2];
rz(-0.1903154) q[2];
sx q[2];
rz(2.9471617) q[2];
rz(-0.91313177) q[3];
sx q[3];
rz(-1.3875995) q[3];
sx q[3];
rz(-2.156179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5450491) q[0];
sx q[0];
rz(-2.5248435) q[0];
sx q[0];
rz(-0.066666691) q[0];
rz(2.8170259) q[1];
sx q[1];
rz(-1.5044183) q[1];
sx q[1];
rz(-2.1527122) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3263771) q[0];
sx q[0];
rz(-1.9003552) q[0];
sx q[0];
rz(-1.1471079) q[0];
rz(2.9724389) q[2];
sx q[2];
rz(-3.0468468) q[2];
sx q[2];
rz(0.46846889) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.40438548) q[1];
sx q[1];
rz(-2.7215241) q[1];
sx q[1];
rz(-1.7350446) q[1];
rz(-pi) q[2];
rz(2.057468) q[3];
sx q[3];
rz(-2.0351962) q[3];
sx q[3];
rz(1.407479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4618335) q[2];
sx q[2];
rz(-2.2448848) q[2];
sx q[2];
rz(2.7339593) q[2];
rz(-2.3729825) q[3];
sx q[3];
rz(-1.8278443) q[3];
sx q[3];
rz(-0.9238981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3826564) q[0];
sx q[0];
rz(-1.3019245) q[0];
sx q[0];
rz(0.60920238) q[0];
rz(0.095104782) q[1];
sx q[1];
rz(-1.8895878) q[1];
sx q[1];
rz(0.87337714) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0489037) q[0];
sx q[0];
rz(-0.19177076) q[0];
sx q[0];
rz(-1.9321241) q[0];
x q[1];
rz(1.7071502) q[2];
sx q[2];
rz(-1.7169723) q[2];
sx q[2];
rz(1.2863976) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9511309) q[1];
sx q[1];
rz(-1.8776263) q[1];
sx q[1];
rz(-2.5224586) q[1];
rz(-pi) q[2];
x q[2];
rz(2.196225) q[3];
sx q[3];
rz(-2.4646467) q[3];
sx q[3];
rz(2.7834746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0218899) q[2];
sx q[2];
rz(-2.3535574) q[2];
sx q[2];
rz(2.4592887) q[2];
rz(-0.37426379) q[3];
sx q[3];
rz(-1.5687317) q[3];
sx q[3];
rz(0.16690978) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69797126) q[0];
sx q[0];
rz(-1.3598096) q[0];
sx q[0];
rz(-0.95296729) q[0];
rz(2.3151746) q[1];
sx q[1];
rz(-2.4024139) q[1];
sx q[1];
rz(1.3964765) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1288426) q[0];
sx q[0];
rz(-1.2484974) q[0];
sx q[0];
rz(-1.4803356) q[0];
rz(2.4834677) q[2];
sx q[2];
rz(-1.8089559) q[2];
sx q[2];
rz(-2.4326774) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9979981) q[1];
sx q[1];
rz(-1.2600139) q[1];
sx q[1];
rz(0.29923156) q[1];
rz(-pi) q[2];
rz(0.26225984) q[3];
sx q[3];
rz(-1.5159303) q[3];
sx q[3];
rz(-2.7454387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.46618) q[2];
sx q[2];
rz(-0.35623494) q[2];
sx q[2];
rz(-0.15979016) q[2];
rz(-0.30188489) q[3];
sx q[3];
rz(-2.2146137) q[3];
sx q[3];
rz(-0.2872428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(-0.940154) q[2];
sx q[2];
rz(-1.7241782) q[2];
sx q[2];
rz(0.45964514) q[2];
rz(3.0575183) q[3];
sx q[3];
rz(-0.53518674) q[3];
sx q[3];
rz(-0.72471602) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
