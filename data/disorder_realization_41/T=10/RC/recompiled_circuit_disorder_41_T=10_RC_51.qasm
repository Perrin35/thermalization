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
rz(-0.2149166) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2673187) q[0];
sx q[0];
rz(-1.312717) q[0];
sx q[0];
rz(-1.8549071) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29834892) q[2];
sx q[2];
rz(-2.6539408) q[2];
sx q[2];
rz(-1.2717441) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3494527) q[1];
sx q[1];
rz(-0.98590241) q[1];
sx q[1];
rz(-2.5930415) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0825726) q[3];
sx q[3];
rz(-0.91592741) q[3];
sx q[3];
rz(2.0522576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.73137838) q[2];
sx q[2];
rz(-1.6822858) q[2];
sx q[2];
rz(0.56420502) q[2];
rz(1.365186) q[3];
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
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0441701) q[0];
sx q[0];
rz(-1.928227) q[0];
sx q[0];
rz(2.2136097) q[0];
rz(-1.1652975) q[1];
sx q[1];
rz(-1.5382643) q[1];
sx q[1];
rz(-0.89675084) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39820406) q[0];
sx q[0];
rz(-2.2687015) q[0];
sx q[0];
rz(-2.8845805) q[0];
x q[1];
rz(-2.7361761) q[2];
sx q[2];
rz(-2.9581796) q[2];
sx q[2];
rz(0.90454067) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.6078728) q[1];
sx q[1];
rz(-2.3405511) q[1];
sx q[1];
rz(-2.968722) q[1];
rz(-pi) q[2];
rz(0.81929042) q[3];
sx q[3];
rz(-1.8835526) q[3];
sx q[3];
rz(-2.6889338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.26560489) q[2];
sx q[2];
rz(-2.6066055) q[2];
sx q[2];
rz(1.0401475) q[2];
rz(1.4552207) q[3];
sx q[3];
rz(-1.8486332) q[3];
sx q[3];
rz(0.35475981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4002157) q[0];
sx q[0];
rz(-0.57755661) q[0];
sx q[0];
rz(-1.0282015) q[0];
rz(-1.0785412) q[1];
sx q[1];
rz(-0.56285793) q[1];
sx q[1];
rz(-2.7064586) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2745278) q[0];
sx q[0];
rz(-1.2194249) q[0];
sx q[0];
rz(-1.4562796) q[0];
rz(-pi) q[1];
rz(-0.25031309) q[2];
sx q[2];
rz(-1.7134943) q[2];
sx q[2];
rz(-0.38562361) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.55089009) q[1];
sx q[1];
rz(-1.0895551) q[1];
sx q[1];
rz(-0.18443702) q[1];
rz(1.9153254) q[3];
sx q[3];
rz(-0.83762729) q[3];
sx q[3];
rz(1.4602349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6083287) q[2];
sx q[2];
rz(-1.8412795) q[2];
sx q[2];
rz(2.8386774) q[2];
rz(1.3251925) q[3];
sx q[3];
rz(-1.9830827) q[3];
sx q[3];
rz(0.091025092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19642297) q[0];
sx q[0];
rz(-1.7315995) q[0];
sx q[0];
rz(2.2241425) q[0];
rz(-0.67287412) q[1];
sx q[1];
rz(-2.0560975) q[1];
sx q[1];
rz(2.8767169) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0850071) q[0];
sx q[0];
rz(-1.5406973) q[0];
sx q[0];
rz(0.65926512) q[0];
rz(2.8633966) q[2];
sx q[2];
rz(-1.858466) q[2];
sx q[2];
rz(-2.1759335) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.781573) q[1];
sx q[1];
rz(-0.85687602) q[1];
sx q[1];
rz(-0.54846958) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81569205) q[3];
sx q[3];
rz(-0.94592735) q[3];
sx q[3];
rz(0.16251414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.778487) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.251579) q[0];
sx q[0];
rz(-0.13680923) q[0];
sx q[0];
rz(-2.662861) q[0];
rz(-2.1084673) q[1];
sx q[1];
rz(-0.9712351) q[1];
sx q[1];
rz(-0.95265257) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1088516) q[0];
sx q[0];
rz(-1.6881588) q[0];
sx q[0];
rz(-2.5928241) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6636159) q[2];
sx q[2];
rz(-2.6307719) q[2];
sx q[2];
rz(-0.36885992) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.488727) q[1];
sx q[1];
rz(-1.1796724) q[1];
sx q[1];
rz(-2.8847242) q[1];
rz(-2.6635366) q[3];
sx q[3];
rz(-1.5336256) q[3];
sx q[3];
rz(0.93070785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4218563) q[2];
sx q[2];
rz(-0.36281261) q[2];
sx q[2];
rz(-0.53058132) q[2];
rz(1.7355708) q[3];
sx q[3];
rz(-1.1281745) q[3];
sx q[3];
rz(-0.83166844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.56753165) q[0];
sx q[0];
rz(-1.4968137) q[0];
sx q[0];
rz(-1.5166327) q[0];
rz(1.8364871) q[1];
sx q[1];
rz(-1.790698) q[1];
sx q[1];
rz(2.9690202) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5694002) q[0];
sx q[0];
rz(-1.5844371) q[0];
sx q[0];
rz(-1.1543247) q[0];
rz(-2.0270945) q[2];
sx q[2];
rz(-1.3367532) q[2];
sx q[2];
rz(-0.32548387) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8166568) q[1];
sx q[1];
rz(-1.4656193) q[1];
sx q[1];
rz(1.9865958) q[1];
rz(-pi) q[2];
rz(2.8302912) q[3];
sx q[3];
rz(-2.4875896) q[3];
sx q[3];
rz(-1.9037387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3036348) q[2];
sx q[2];
rz(-1.4540318) q[2];
sx q[2];
rz(-1.1266358) q[2];
rz(-0.78222328) q[3];
sx q[3];
rz(-1.2354847) q[3];
sx q[3];
rz(1.8036028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.85309) q[0];
sx q[0];
rz(-2.8350916) q[0];
sx q[0];
rz(0.66147584) q[0];
rz(-2.181197) q[1];
sx q[1];
rz(-1.7405225) q[1];
sx q[1];
rz(-2.3849934) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3698759) q[0];
sx q[0];
rz(-1.6635832) q[0];
sx q[0];
rz(1.7232399) q[0];
rz(-pi) q[1];
rz(-0.43811626) q[2];
sx q[2];
rz(-0.44898673) q[2];
sx q[2];
rz(-0.24030906) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7797864) q[1];
sx q[1];
rz(-2.8621319) q[1];
sx q[1];
rz(-0.96868412) q[1];
rz(-pi) q[2];
rz(-2.9499801) q[3];
sx q[3];
rz(-1.629047) q[3];
sx q[3];
rz(-2.2239457) q[3];
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
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5450491) q[0];
sx q[0];
rz(-0.61674917) q[0];
sx q[0];
rz(3.074926) q[0];
rz(-0.32456675) q[1];
sx q[1];
rz(-1.5044183) q[1];
sx q[1];
rz(0.98888046) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37799997) q[0];
sx q[0];
rz(-0.53056301) q[0];
sx q[0];
rz(0.87688045) q[0];
rz(-pi) q[1];
rz(2.9724389) q[2];
sx q[2];
rz(-0.094745853) q[2];
sx q[2];
rz(2.6731238) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1253742) q[1];
sx q[1];
rz(-1.5040633) q[1];
sx q[1];
rz(-1.9858422) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.057468) q[3];
sx q[3];
rz(-1.1063965) q[3];
sx q[3];
rz(-1.7341136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4618335) q[2];
sx q[2];
rz(-0.89670783) q[2];
sx q[2];
rz(-2.7339593) q[2];
rz(-2.3729825) q[3];
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
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3826564) q[0];
sx q[0];
rz(-1.3019245) q[0];
sx q[0];
rz(2.5323903) q[0];
rz(-0.095104782) q[1];
sx q[1];
rz(-1.2520049) q[1];
sx q[1];
rz(0.87337714) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3082335) q[0];
sx q[0];
rz(-1.6382268) q[0];
sx q[0];
rz(1.391135) q[0];
rz(-pi) q[1];
rz(-0.74553211) q[2];
sx q[2];
rz(-2.9420256) q[2];
sx q[2];
rz(-1.0996639) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9511309) q[1];
sx q[1];
rz(-1.2639664) q[1];
sx q[1];
rz(0.61913403) q[1];
rz(-2.196225) q[3];
sx q[3];
rz(-2.4646467) q[3];
sx q[3];
rz(0.35811801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1197027) q[2];
sx q[2];
rz(-2.3535574) q[2];
sx q[2];
rz(-0.68230391) q[2];
rz(0.37426379) q[3];
sx q[3];
rz(-1.5687317) q[3];
sx q[3];
rz(2.9746829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6122702) q[0];
sx q[0];
rz(-1.6565874) q[0];
sx q[0];
rz(0.32353185) q[0];
rz(-1.868532) q[2];
sx q[2];
rz(-0.93431384) q[2];
sx q[2];
rz(-1.0422848) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4954088) q[1];
sx q[1];
rz(-2.7135661) q[1];
sx q[1];
rz(2.3133548) q[1];
x q[2];
rz(2.8793328) q[3];
sx q[3];
rz(-1.6256623) q[3];
sx q[3];
rz(-2.7454387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.46618) q[2];
sx q[2];
rz(-0.35623494) q[2];
sx q[2];
rz(-2.9818025) q[2];
rz(2.8397078) q[3];
sx q[3];
rz(-0.92697898) q[3];
sx q[3];
rz(-2.8543499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.044534279) q[0];
sx q[0];
rz(-0.67561588) q[0];
sx q[0];
rz(-1.5560879) q[0];
rz(0.13327577) q[1];
sx q[1];
rz(-1.6242846) q[1];
sx q[1];
rz(-0.12856738) q[1];
rz(-2.2014387) q[2];
sx q[2];
rz(-1.4174145) q[2];
sx q[2];
rz(-2.6819475) q[2];
rz(-3.0575183) q[3];
sx q[3];
rz(-2.6064059) q[3];
sx q[3];
rz(2.4168766) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
