OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0410864) q[0];
sx q[0];
rz(-0.20809986) q[0];
sx q[0];
rz(2.7383374) q[0];
rz(2.9085605) q[1];
sx q[1];
rz(-1.7014528) q[1];
sx q[1];
rz(-0.22388248) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3849523) q[0];
sx q[0];
rz(-2.1813574) q[0];
sx q[0];
rz(-1.8288307) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86647948) q[2];
sx q[2];
rz(-1.472855) q[2];
sx q[2];
rz(-0.78664727) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.63131166) q[1];
sx q[1];
rz(-1.3461539) q[1];
sx q[1];
rz(-0.90490492) q[1];
rz(0.68239642) q[3];
sx q[3];
rz(-0.93554893) q[3];
sx q[3];
rz(1.531383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4441121) q[2];
sx q[2];
rz(-2.8287502) q[2];
sx q[2];
rz(3.1057788) q[2];
rz(0.73016417) q[3];
sx q[3];
rz(-1.1532447) q[3];
sx q[3];
rz(-0.20619503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8908454) q[0];
sx q[0];
rz(-0.99531168) q[0];
sx q[0];
rz(-0.19749755) q[0];
rz(-2.6644871) q[1];
sx q[1];
rz(-0.66480607) q[1];
sx q[1];
rz(0.8055996) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8715848) q[0];
sx q[0];
rz(-2.2257937) q[0];
sx q[0];
rz(1.3359469) q[0];
rz(2.9764011) q[2];
sx q[2];
rz(-2.7562263) q[2];
sx q[2];
rz(-1.8866273) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2195314) q[1];
sx q[1];
rz(-2.1967255) q[1];
sx q[1];
rz(-0.74447592) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74975852) q[3];
sx q[3];
rz(-2.9364412) q[3];
sx q[3];
rz(1.4914025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3104559) q[2];
sx q[2];
rz(-1.6418991) q[2];
sx q[2];
rz(1.4031225) q[2];
rz(-0.63203114) q[3];
sx q[3];
rz(-2.8681614) q[3];
sx q[3];
rz(3.0145751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8359351) q[0];
sx q[0];
rz(-0.62863612) q[0];
sx q[0];
rz(2.054731) q[0];
rz(2.944259) q[1];
sx q[1];
rz(-1.5060164) q[1];
sx q[1];
rz(2.8089583) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64651782) q[0];
sx q[0];
rz(-1.7103638) q[0];
sx q[0];
rz(-2.0049176) q[0];
rz(-pi) q[1];
rz(-2.3057322) q[2];
sx q[2];
rz(-1.5243013) q[2];
sx q[2];
rz(2.495386) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3217993) q[1];
sx q[1];
rz(-1.7463356) q[1];
sx q[1];
rz(1.3513397) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8007711) q[3];
sx q[3];
rz(-0.85715961) q[3];
sx q[3];
rz(-0.35060397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3745594) q[2];
sx q[2];
rz(-1.5908396) q[2];
sx q[2];
rz(-1.4683051) q[2];
rz(3.1324006) q[3];
sx q[3];
rz(-2.6720948) q[3];
sx q[3];
rz(-0.79206842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62995768) q[0];
sx q[0];
rz(-0.63794962) q[0];
sx q[0];
rz(-1.6988423) q[0];
rz(-2.1171782) q[1];
sx q[1];
rz(-1.2316278) q[1];
sx q[1];
rz(2.3146497) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0684375) q[0];
sx q[0];
rz(-0.19864635) q[0];
sx q[0];
rz(-1.1852926) q[0];
rz(-0.154279) q[2];
sx q[2];
rz(-2.2529229) q[2];
sx q[2];
rz(-0.23274225) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6232439) q[1];
sx q[1];
rz(-2.001013) q[1];
sx q[1];
rz(1.3395549) q[1];
rz(-pi) q[2];
rz(-0.8884646) q[3];
sx q[3];
rz(-1.6649705) q[3];
sx q[3];
rz(1.7674131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3247165) q[2];
sx q[2];
rz(-0.78410316) q[2];
sx q[2];
rz(1.8640222) q[2];
rz(-0.68814021) q[3];
sx q[3];
rz(-1.0253996) q[3];
sx q[3];
rz(2.1791606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14030305) q[0];
sx q[0];
rz(-1.7004509) q[0];
sx q[0];
rz(0.21743123) q[0];
rz(0.28542074) q[1];
sx q[1];
rz(-2.5358584) q[1];
sx q[1];
rz(-0.4981471) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0764387) q[0];
sx q[0];
rz(-2.3992043) q[0];
sx q[0];
rz(-2.0641293) q[0];
x q[1];
rz(2.010538) q[2];
sx q[2];
rz(-2.5806107) q[2];
sx q[2];
rz(1.9064685) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5707055) q[1];
sx q[1];
rz(-1.0859617) q[1];
sx q[1];
rz(-0.2356694) q[1];
rz(-pi) q[2];
rz(-2.7015649) q[3];
sx q[3];
rz(-1.3413359) q[3];
sx q[3];
rz(-0.10073951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6238326) q[2];
sx q[2];
rz(-2.0443003) q[2];
sx q[2];
rz(-0.1499873) q[2];
rz(0.013966694) q[3];
sx q[3];
rz(-1.59168) q[3];
sx q[3];
rz(-0.51378957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.49195313) q[0];
sx q[0];
rz(-2.5998901) q[0];
sx q[0];
rz(0.81277043) q[0];
rz(1.4179519) q[1];
sx q[1];
rz(-1.6287454) q[1];
sx q[1];
rz(-1.0636122) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2246216) q[0];
sx q[0];
rz(-1.2422891) q[0];
sx q[0];
rz(-0.55060951) q[0];
x q[1];
rz(0.10165693) q[2];
sx q[2];
rz(-0.16451193) q[2];
sx q[2];
rz(-0.042334231) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9199976) q[1];
sx q[1];
rz(-2.1694899) q[1];
sx q[1];
rz(1.5569975) q[1];
rz(-pi) q[2];
rz(-2.205414) q[3];
sx q[3];
rz(-2.8293316) q[3];
sx q[3];
rz(-2.0978417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.873988) q[2];
sx q[2];
rz(-2.5817817) q[2];
sx q[2];
rz(1.7248636) q[2];
rz(0.060976107) q[3];
sx q[3];
rz(-2.2305326) q[3];
sx q[3];
rz(-1.6772259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7199719) q[0];
sx q[0];
rz(-0.82141972) q[0];
sx q[0];
rz(-3.022497) q[0];
rz(-2.6630317) q[1];
sx q[1];
rz(-2.3275972) q[1];
sx q[1];
rz(-2.4518769) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0499303) q[0];
sx q[0];
rz(-2.1436467) q[0];
sx q[0];
rz(0.86688231) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5096171) q[2];
sx q[2];
rz(-0.58210974) q[2];
sx q[2];
rz(1.7012973) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7642245) q[1];
sx q[1];
rz(-0.86478327) q[1];
sx q[1];
rz(-3.0824667) q[1];
rz(1.4004565) q[3];
sx q[3];
rz(-2.4217477) q[3];
sx q[3];
rz(0.094739044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0488284) q[2];
sx q[2];
rz(-3.0807639) q[2];
sx q[2];
rz(-2.0737341) q[2];
rz(0.16718665) q[3];
sx q[3];
rz(-1.7696295) q[3];
sx q[3];
rz(-3.0021477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6086513) q[0];
sx q[0];
rz(-0.81357384) q[0];
sx q[0];
rz(2.2331878) q[0];
rz(-0.99336973) q[1];
sx q[1];
rz(-2.9790331) q[1];
sx q[1];
rz(2.2672674) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68087606) q[0];
sx q[0];
rz(-1.5534288) q[0];
sx q[0];
rz(-0.13219035) q[0];
rz(-1.4972009) q[2];
sx q[2];
rz(-2.8831867) q[2];
sx q[2];
rz(-2.0143353) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5914375) q[1];
sx q[1];
rz(-2.9936325) q[1];
sx q[1];
rz(1.9846538) q[1];
rz(-2.1217974) q[3];
sx q[3];
rz(-1.2606818) q[3];
sx q[3];
rz(-1.9581025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.61451644) q[2];
sx q[2];
rz(-0.15915844) q[2];
sx q[2];
rz(-2.2873774) q[2];
rz(2.6863875) q[3];
sx q[3];
rz(-2.5762317) q[3];
sx q[3];
rz(0.2555041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.663317) q[0];
sx q[0];
rz(-1.9238967) q[0];
sx q[0];
rz(-1.664337) q[0];
rz(-1.5746337) q[1];
sx q[1];
rz(-2.4576371) q[1];
sx q[1];
rz(2.4050567) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7304717) q[0];
sx q[0];
rz(-2.6101042) q[0];
sx q[0];
rz(2.0474252) q[0];
x q[1];
rz(1.4669424) q[2];
sx q[2];
rz(-1.9496634) q[2];
sx q[2];
rz(-1.9571638) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.29727259) q[1];
sx q[1];
rz(-1.2579009) q[1];
sx q[1];
rz(-1.3144397) q[1];
rz(2.1677954) q[3];
sx q[3];
rz(-1.5415808) q[3];
sx q[3];
rz(-2.4643498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7443202) q[2];
sx q[2];
rz(-0.87817764) q[2];
sx q[2];
rz(1.5273904) q[2];
rz(0.39198908) q[3];
sx q[3];
rz(-1.1992998) q[3];
sx q[3];
rz(-3.0919302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.024260661) q[0];
sx q[0];
rz(-1.4735104) q[0];
sx q[0];
rz(-2.9397553) q[0];
rz(-0.2233389) q[1];
sx q[1];
rz(-0.52450648) q[1];
sx q[1];
rz(-0.39252678) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46014547) q[0];
sx q[0];
rz(-2.0325615) q[0];
sx q[0];
rz(-0.11790922) q[0];
x q[1];
rz(0.052458737) q[2];
sx q[2];
rz(-0.83485583) q[2];
sx q[2];
rz(-0.8142161) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1087677) q[1];
sx q[1];
rz(-0.6995753) q[1];
sx q[1];
rz(-2.2684069) q[1];
rz(-pi) q[2];
x q[2];
rz(0.48247561) q[3];
sx q[3];
rz(-1.2996718) q[3];
sx q[3];
rz(1.9737873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4203804) q[2];
sx q[2];
rz(-2.1537697) q[2];
sx q[2];
rz(0.58050275) q[2];
rz(2.765559) q[3];
sx q[3];
rz(-2.1707363) q[3];
sx q[3];
rz(1.6002801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15039438) q[0];
sx q[0];
rz(-1.5924441) q[0];
sx q[0];
rz(1.7572255) q[0];
rz(1.4906384) q[1];
sx q[1];
rz(-1.2532267) q[1];
sx q[1];
rz(0.86029235) q[1];
rz(-2.2223086) q[2];
sx q[2];
rz(-2.2625661) q[2];
sx q[2];
rz(1.7023466) q[2];
rz(-1.5483472) q[3];
sx q[3];
rz(-1.8825363) q[3];
sx q[3];
rz(-2.0288449) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
