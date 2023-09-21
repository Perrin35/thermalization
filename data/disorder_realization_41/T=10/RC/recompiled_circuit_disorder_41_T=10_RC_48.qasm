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
rz(2.7472189) q[1];
sx q[1];
rz(-1.4596649) q[1];
sx q[1];
rz(2.9266761) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7267933) q[0];
sx q[0];
rz(-0.3814632) q[0];
sx q[0];
rz(-2.3261855) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6724042) q[2];
sx q[2];
rz(-1.708963) q[2];
sx q[2];
rz(-0.56433041) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.51337459) q[1];
sx q[1];
rz(-2.3623423) q[1];
sx q[1];
rz(-2.2378504) q[1];
x q[2];
rz(1.0825726) q[3];
sx q[3];
rz(-2.2256652) q[3];
sx q[3];
rz(1.0893351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.73137838) q[2];
sx q[2];
rz(-1.4593068) q[2];
sx q[2];
rz(-0.56420502) q[2];
rz(-1.365186) q[3];
sx q[3];
rz(-2.6919638) q[3];
sx q[3];
rz(-1.8723429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0441701) q[0];
sx q[0];
rz(-1.2133657) q[0];
sx q[0];
rz(0.92798293) q[0];
rz(1.1652975) q[1];
sx q[1];
rz(-1.6033283) q[1];
sx q[1];
rz(2.2448418) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7433886) q[0];
sx q[0];
rz(-0.87289116) q[0];
sx q[0];
rz(0.25701216) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7361761) q[2];
sx q[2];
rz(-0.18341309) q[2];
sx q[2];
rz(0.90454067) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0577382) q[1];
sx q[1];
rz(-1.4469622) q[1];
sx q[1];
rz(2.3480575) q[1];
rz(2.3223022) q[3];
sx q[3];
rz(-1.8835526) q[3];
sx q[3];
rz(-0.4526588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.26560489) q[2];
sx q[2];
rz(-0.53498712) q[2];
sx q[2];
rz(-2.1014452) q[2];
rz(1.6863719) q[3];
sx q[3];
rz(-1.8486332) q[3];
sx q[3];
rz(2.7868328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4002157) q[0];
sx q[0];
rz(-2.564036) q[0];
sx q[0];
rz(2.1133912) q[0];
rz(-1.0785412) q[1];
sx q[1];
rz(-0.56285793) q[1];
sx q[1];
rz(0.43513402) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4774287) q[0];
sx q[0];
rz(-1.6782883) q[0];
sx q[0];
rz(0.35350032) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8912796) q[2];
sx q[2];
rz(-1.4280983) q[2];
sx q[2];
rz(-2.755969) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.55089009) q[1];
sx q[1];
rz(-2.0520376) q[1];
sx q[1];
rz(0.18443702) q[1];
x q[2];
rz(2.3782303) q[3];
sx q[3];
rz(-1.3171139) q[3];
sx q[3];
rz(-0.34624472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6083287) q[2];
sx q[2];
rz(-1.8412795) q[2];
sx q[2];
rz(2.8386774) q[2];
rz(-1.8164002) q[3];
sx q[3];
rz(-1.15851) q[3];
sx q[3];
rz(3.0505676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9451697) q[0];
sx q[0];
rz(-1.7315995) q[0];
sx q[0];
rz(0.91745013) q[0];
rz(2.4687185) q[1];
sx q[1];
rz(-2.0560975) q[1];
sx q[1];
rz(2.8767169) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6040651) q[0];
sx q[0];
rz(-2.2297105) q[0];
sx q[0];
rz(-1.6088681) q[0];
rz(-pi) q[1];
rz(2.3189544) q[2];
sx q[2];
rz(-0.3974786) q[2];
sx q[2];
rz(2.9646404) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.6092104) q[1];
sx q[1];
rz(-0.86984837) q[1];
sx q[1];
rz(1.0290531) q[1];
rz(-2.3817252) q[3];
sx q[3];
rz(-2.2025975) q[3];
sx q[3];
rz(2.2894273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.36310568) q[2];
sx q[2];
rz(-2.6532756) q[2];
sx q[2];
rz(-1.5765566) q[2];
rz(1.0270843) q[3];
sx q[3];
rz(-2.4001207) q[3];
sx q[3];
rz(-2.0402133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.251579) q[0];
sx q[0];
rz(-0.13680923) q[0];
sx q[0];
rz(0.47873163) q[0];
rz(1.0331253) q[1];
sx q[1];
rz(-2.1703576) q[1];
sx q[1];
rz(-2.1889401) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.032741) q[0];
sx q[0];
rz(-1.4534338) q[0];
sx q[0];
rz(-2.5928241) q[0];
rz(-pi) q[1];
rz(-0.41575899) q[2];
sx q[2];
rz(-1.8766878) q[2];
sx q[2];
rz(1.3410459) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2561803) q[1];
sx q[1];
rz(-2.6773239) q[1];
sx q[1];
rz(-2.1229565) q[1];
rz(1.6126552) q[3];
sx q[3];
rz(-1.0930982) q[3];
sx q[3];
rz(0.65934138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4218563) q[2];
sx q[2];
rz(-0.36281261) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.574061) q[0];
sx q[0];
rz(-1.644779) q[0];
sx q[0];
rz(-1.5166327) q[0];
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
rz(2.1369551) q[0];
sx q[0];
rz(-1.9872268) q[0];
sx q[0];
rz(0.01491551) q[0];
rz(-pi) q[1];
rz(-0.2595915) q[2];
sx q[2];
rz(-2.0137557) q[2];
sx q[2];
rz(-1.7829347) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8166568) q[1];
sx q[1];
rz(-1.6759733) q[1];
sx q[1];
rz(1.9865958) q[1];
rz(1.8014088) q[3];
sx q[3];
rz(-2.1884544) q[3];
sx q[3];
rz(-2.2889083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3036348) q[2];
sx q[2];
rz(-1.6875608) q[2];
sx q[2];
rz(1.1266358) q[2];
rz(2.3593694) q[3];
sx q[3];
rz(-1.9061079) q[3];
sx q[3];
rz(1.3379898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.85309) q[0];
sx q[0];
rz(-2.8350916) q[0];
sx q[0];
rz(-0.66147584) q[0];
rz(-0.96039564) q[1];
sx q[1];
rz(-1.7405225) q[1];
sx q[1];
rz(2.3849934) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3698759) q[0];
sx q[0];
rz(-1.4780095) q[0];
sx q[0];
rz(1.4183527) q[0];
rz(-pi) q[1];
x q[1];
rz(0.43811626) q[2];
sx q[2];
rz(-2.6926059) q[2];
sx q[2];
rz(-0.24030906) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3488256) q[1];
sx q[1];
rz(-1.7276689) q[1];
sx q[1];
rz(-1.3385593) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29715085) q[3];
sx q[3];
rz(-0.20016709) q[3];
sx q[3];
rz(2.196892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1371655) q[2];
sx q[2];
rz(-0.1903154) q[2];
sx q[2];
rz(-2.9471617) q[2];
rz(0.91313177) q[3];
sx q[3];
rz(-1.7539932) q[3];
sx q[3];
rz(-2.156179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5450491) q[0];
sx q[0];
rz(-2.5248435) q[0];
sx q[0];
rz(3.074926) q[0];
rz(-2.8170259) q[1];
sx q[1];
rz(-1.6371744) q[1];
sx q[1];
rz(0.98888046) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37799997) q[0];
sx q[0];
rz(-2.6110296) q[0];
sx q[0];
rz(2.2647122) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0481911) q[2];
sx q[2];
rz(-1.5548692) q[2];
sx q[2];
rz(1.2707368) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.40438548) q[1];
sx q[1];
rz(-2.7215241) q[1];
sx q[1];
rz(-1.7350446) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51560651) q[3];
sx q[3];
rz(-1.139384) q[3];
sx q[3];
rz(-0.069375667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4618335) q[2];
sx q[2];
rz(-0.89670783) q[2];
sx q[2];
rz(-0.40763339) q[2];
rz(2.3729825) q[3];
sx q[3];
rz(-1.8278443) q[3];
sx q[3];
rz(-2.2176946) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3826564) q[0];
sx q[0];
rz(-1.8396682) q[0];
sx q[0];
rz(-0.60920238) q[0];
rz(-0.095104782) q[1];
sx q[1];
rz(-1.8895878) q[1];
sx q[1];
rz(-0.87337714) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3082335) q[0];
sx q[0];
rz(-1.5033659) q[0];
sx q[0];
rz(-1.7504577) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.74553211) q[2];
sx q[2];
rz(-2.9420256) q[2];
sx q[2];
rz(-1.0996639) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1904618) q[1];
sx q[1];
rz(-1.8776263) q[1];
sx q[1];
rz(-0.61913403) q[1];
rz(2.196225) q[3];
sx q[3];
rz(-0.67694596) q[3];
sx q[3];
rz(-2.7834746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0218899) q[2];
sx q[2];
rz(-0.78803524) q[2];
sx q[2];
rz(0.68230391) q[2];
rz(0.37426379) q[3];
sx q[3];
rz(-1.5687317) q[3];
sx q[3];
rz(-0.16690978) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4436214) q[0];
sx q[0];
rz(-1.3598096) q[0];
sx q[0];
rz(-0.95296729) q[0];
rz(-0.8264181) q[1];
sx q[1];
rz(-2.4024139) q[1];
sx q[1];
rz(-1.7451161) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29165927) q[0];
sx q[0];
rz(-0.33432654) q[0];
sx q[0];
rz(2.8773984) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2730607) q[2];
sx q[2];
rz(-0.93431384) q[2];
sx q[2];
rz(2.0993078) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4954088) q[1];
sx q[1];
rz(-0.42802654) q[1];
sx q[1];
rz(2.3133548) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8793328) q[3];
sx q[3];
rz(-1.6256623) q[3];
sx q[3];
rz(2.7454387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.46618) q[2];
sx q[2];
rz(-2.7853577) q[2];
sx q[2];
rz(0.15979016) q[2];
rz(-2.8397078) q[3];
sx q[3];
rz(-0.92697898) q[3];
sx q[3];
rz(-0.2872428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
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
rz(0.13327577) q[1];
sx q[1];
rz(-1.6242846) q[1];
sx q[1];
rz(-0.12856738) q[1];
rz(-1.8272022) q[2];
sx q[2];
rz(-0.64654965) q[2];
sx q[2];
rz(1.8241573) q[2];
rz(-1.521048) q[3];
sx q[3];
rz(-1.0377025) q[3];
sx q[3];
rz(2.514537) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
