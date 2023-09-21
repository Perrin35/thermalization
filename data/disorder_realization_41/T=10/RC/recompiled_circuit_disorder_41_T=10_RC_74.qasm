OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.4364606) q[0];
sx q[0];
rz(-0.55186614) q[0];
sx q[0];
rz(-3.119757) q[0];
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
rz(2.7267933) q[0];
sx q[0];
rz(-0.3814632) q[0];
sx q[0];
rz(-2.3261855) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6724042) q[2];
sx q[2];
rz(-1.4326296) q[2];
sx q[2];
rz(-0.56433041) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3494527) q[1];
sx q[1];
rz(-0.98590241) q[1];
sx q[1];
rz(0.5485512) q[1];
x q[2];
rz(-2.5932556) q[3];
sx q[3];
rz(-0.79474802) q[3];
sx q[3];
rz(-0.37219513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.73137838) q[2];
sx q[2];
rz(-1.6822858) q[2];
sx q[2];
rz(2.5773876) q[2];
rz(-1.7764067) q[3];
sx q[3];
rz(-0.44962883) q[3];
sx q[3];
rz(1.2692497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0441701) q[0];
sx q[0];
rz(-1.2133657) q[0];
sx q[0];
rz(-0.92798293) q[0];
rz(1.9762951) q[1];
sx q[1];
rz(-1.5382643) q[1];
sx q[1];
rz(2.2448418) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7433886) q[0];
sx q[0];
rz(-0.87289116) q[0];
sx q[0];
rz(2.8845805) q[0];
x q[1];
rz(-0.40541655) q[2];
sx q[2];
rz(-0.18341309) q[2];
sx q[2];
rz(-2.237052) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7795418) q[1];
sx q[1];
rz(-0.78501399) q[1];
sx q[1];
rz(-1.3951468) q[1];
rz(2.3223022) q[3];
sx q[3];
rz(-1.8835526) q[3];
sx q[3];
rz(-0.4526588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8759878) q[2];
sx q[2];
rz(-2.6066055) q[2];
sx q[2];
rz(-1.0401475) q[2];
rz(1.6863719) q[3];
sx q[3];
rz(-1.8486332) q[3];
sx q[3];
rz(-0.35475981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(2.4002157) q[0];
sx q[0];
rz(-2.564036) q[0];
sx q[0];
rz(1.0282015) q[0];
rz(1.0785412) q[1];
sx q[1];
rz(-0.56285793) q[1];
sx q[1];
rz(-0.43513402) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4774287) q[0];
sx q[0];
rz(-1.6782883) q[0];
sx q[0];
rz(0.35350032) q[0];
x q[1];
rz(-1.7180195) q[2];
sx q[2];
rz(-1.3230811) q[2];
sx q[2];
rz(-1.9927646) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0355465) q[1];
sx q[1];
rz(-1.4075081) q[1];
sx q[1];
rz(-1.0825023) q[1];
rz(-0.76336236) q[3];
sx q[3];
rz(-1.3171139) q[3];
sx q[3];
rz(2.7953479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.53326398) q[2];
sx q[2];
rz(-1.3003131) q[2];
sx q[2];
rz(-2.8386774) q[2];
rz(-1.8164002) q[3];
sx q[3];
rz(-1.9830827) q[3];
sx q[3];
rz(-3.0505676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19642297) q[0];
sx q[0];
rz(-1.4099932) q[0];
sx q[0];
rz(-0.91745013) q[0];
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
rz(-0.53752758) q[0];
sx q[0];
rz(-2.2297105) q[0];
sx q[0];
rz(1.6088681) q[0];
rz(-pi) q[1];
rz(2.8633966) q[2];
sx q[2];
rz(-1.858466) q[2];
sx q[2];
rz(0.96565914) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5323822) q[1];
sx q[1];
rz(-0.86984837) q[1];
sx q[1];
rz(2.1125395) q[1];
rz(-pi) q[2];
rz(-0.81569205) q[3];
sx q[3];
rz(-0.94592735) q[3];
sx q[3];
rz(0.16251414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.778487) q[2];
sx q[2];
rz(-0.48831707) q[2];
sx q[2];
rz(-1.5765566) q[2];
rz(1.0270843) q[3];
sx q[3];
rz(-2.4001207) q[3];
sx q[3];
rz(1.1013793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.89001369) q[0];
sx q[0];
rz(-3.0047834) q[0];
sx q[0];
rz(2.662861) q[0];
rz(1.0331253) q[1];
sx q[1];
rz(-2.1703576) q[1];
sx q[1];
rz(0.95265257) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.032741) q[0];
sx q[0];
rz(-1.4534338) q[0];
sx q[0];
rz(-0.54876859) q[0];
rz(1.2383934) q[2];
sx q[2];
rz(-1.1754416) q[2];
sx q[2];
rz(-0.3619286) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6528656) q[1];
sx q[1];
rz(-1.9619202) q[1];
sx q[1];
rz(-2.8847242) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0609344) q[3];
sx q[3];
rz(-0.47938743) q[3];
sx q[3];
rz(-0.56848923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4218563) q[2];
sx q[2];
rz(-2.77878) q[2];
sx q[2];
rz(-2.6110113) q[2];
rz(1.7355708) q[3];
sx q[3];
rz(-1.1281745) q[3];
sx q[3];
rz(2.3099242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.574061) q[0];
sx q[0];
rz(-1.4968137) q[0];
sx q[0];
rz(-1.6249599) q[0];
rz(1.3051055) q[1];
sx q[1];
rz(-1.790698) q[1];
sx q[1];
rz(-2.9690202) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0046376) q[0];
sx q[0];
rz(-1.1543659) q[0];
sx q[0];
rz(0.01491551) q[0];
rz(2.0270945) q[2];
sx q[2];
rz(-1.8048394) q[2];
sx q[2];
rz(-0.32548387) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8166568) q[1];
sx q[1];
rz(-1.4656193) q[1];
sx q[1];
rz(-1.9865958) q[1];
rz(-pi) q[2];
rz(-1.3401838) q[3];
sx q[3];
rz(-2.1884544) q[3];
sx q[3];
rz(-2.2889083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3036348) q[2];
sx q[2];
rz(-1.4540318) q[2];
sx q[2];
rz(2.0149569) q[2];
rz(-0.78222328) q[3];
sx q[3];
rz(-1.9061079) q[3];
sx q[3];
rz(-1.8036028) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.85309) q[0];
sx q[0];
rz(-0.30650109) q[0];
sx q[0];
rz(-0.66147584) q[0];
rz(0.96039564) q[1];
sx q[1];
rz(-1.7405225) q[1];
sx q[1];
rz(0.75659928) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.341757) q[0];
sx q[0];
rz(-2.9633187) q[0];
sx q[0];
rz(1.0210277) q[0];
x q[1];
rz(-1.3691749) q[2];
sx q[2];
rz(-1.1668418) q[2];
sx q[2];
rz(0.23922761) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.25890299) q[1];
sx q[1];
rz(-1.3414624) q[1];
sx q[1];
rz(2.980466) q[1];
x q[2];
rz(1.6301304) q[3];
sx q[3];
rz(-1.3795128) q[3];
sx q[3];
rz(-0.64185601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1371655) q[2];
sx q[2];
rz(-2.9512773) q[2];
sx q[2];
rz(-0.19443092) q[2];
rz(-2.2284609) q[3];
sx q[3];
rz(-1.7539932) q[3];
sx q[3];
rz(0.98541361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5450491) q[0];
sx q[0];
rz(-2.5248435) q[0];
sx q[0];
rz(3.074926) q[0];
rz(-0.32456675) q[1];
sx q[1];
rz(-1.5044183) q[1];
sx q[1];
rz(0.98888046) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7522404) q[0];
sx q[0];
rz(-1.9703431) q[0];
sx q[0];
rz(-0.35895343) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9724389) q[2];
sx q[2];
rz(-0.094745853) q[2];
sx q[2];
rz(0.46846889) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1253742) q[1];
sx q[1];
rz(-1.6375293) q[1];
sx q[1];
rz(-1.1557505) q[1];
rz(-2.6259861) q[3];
sx q[3];
rz(-1.139384) q[3];
sx q[3];
rz(0.069375667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6797592) q[2];
sx q[2];
rz(-2.2448848) q[2];
sx q[2];
rz(-0.40763339) q[2];
rz(0.76861012) q[3];
sx q[3];
rz(-1.8278443) q[3];
sx q[3];
rz(-0.9238981) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
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
rz(-0.60920238) q[0];
rz(3.0464879) q[1];
sx q[1];
rz(-1.2520049) q[1];
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
rz(-1.7500449) q[0];
sx q[0];
rz(0.068530131) q[0];
x q[1];
rz(0.74553211) q[2];
sx q[2];
rz(-0.19956707) q[2];
sx q[2];
rz(-1.0996639) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1904618) q[1];
sx q[1];
rz(-1.2639664) q[1];
sx q[1];
rz(2.5224586) q[1];
x q[2];
rz(-0.99336266) q[3];
sx q[3];
rz(-1.9462898) q[3];
sx q[3];
rz(1.7253699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0218899) q[2];
sx q[2];
rz(-2.3535574) q[2];
sx q[2];
rz(-0.68230391) q[2];
rz(2.7673289) q[3];
sx q[3];
rz(-1.5687317) q[3];
sx q[3];
rz(-2.9746829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69797126) q[0];
sx q[0];
rz(-1.781783) q[0];
sx q[0];
rz(2.1886254) q[0];
rz(0.8264181) q[1];
sx q[1];
rz(-2.4024139) q[1];
sx q[1];
rz(-1.3964765) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6122702) q[0];
sx q[0];
rz(-1.4850052) q[0];
sx q[0];
rz(-2.8180608) q[0];
rz(-pi) q[1];
rz(-1.2730607) q[2];
sx q[2];
rz(-0.93431384) q[2];
sx q[2];
rz(-2.0993078) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4954088) q[1];
sx q[1];
rz(-2.7135661) q[1];
sx q[1];
rz(-0.82823786) q[1];
x q[2];
rz(-0.26225984) q[3];
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
rz(0.46618) q[2];
sx q[2];
rz(-0.35623494) q[2];
sx q[2];
rz(0.15979016) q[2];
rz(-2.8397078) q[3];
sx q[3];
rz(-0.92697898) q[3];
sx q[3];
rz(2.8543499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044534279) q[0];
sx q[0];
rz(-2.4659768) q[0];
sx q[0];
rz(1.5855047) q[0];
rz(0.13327577) q[1];
sx q[1];
rz(-1.6242846) q[1];
sx q[1];
rz(-0.12856738) q[1];
rz(1.8272022) q[2];
sx q[2];
rz(-2.495043) q[2];
sx q[2];
rz(-1.3174353) q[2];
rz(1.6205447) q[3];
sx q[3];
rz(-1.0377025) q[3];
sx q[3];
rz(2.514537) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
