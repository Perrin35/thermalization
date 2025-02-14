OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3833157) q[0];
sx q[0];
rz(-1.436469) q[0];
sx q[0];
rz(-2.1767148) q[0];
rz(-0.75196737) q[1];
sx q[1];
rz(5.8536018) q[1];
sx q[1];
rz(11.516973) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9193662) q[0];
sx q[0];
rz(-2.1458186) q[0];
sx q[0];
rz(2.7336043) q[0];
x q[1];
rz(-1.8951178) q[2];
sx q[2];
rz(-1.7056381) q[2];
sx q[2];
rz(0.91712778) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.2856372) q[1];
sx q[1];
rz(-2.1826943) q[1];
sx q[1];
rz(2.6735071) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14127381) q[3];
sx q[3];
rz(-2.5577099) q[3];
sx q[3];
rz(-0.48030765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.45399484) q[2];
sx q[2];
rz(-1.5291841) q[2];
sx q[2];
rz(-3.122984) q[2];
rz(-2.5971557) q[3];
sx q[3];
rz(-0.33014044) q[3];
sx q[3];
rz(-2.4936567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7740771) q[0];
sx q[0];
rz(-1.6870966) q[0];
sx q[0];
rz(0.53502214) q[0];
rz(-2.6016443) q[1];
sx q[1];
rz(-0.63553634) q[1];
sx q[1];
rz(-1.3998869) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2311418) q[0];
sx q[0];
rz(-0.60908356) q[0];
sx q[0];
rz(-0.85394359) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82615699) q[2];
sx q[2];
rz(-0.75010502) q[2];
sx q[2];
rz(0.78924417) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3185841) q[1];
sx q[1];
rz(-2.6806147) q[1];
sx q[1];
rz(0.58580841) q[1];
x q[2];
rz(-1.0008282) q[3];
sx q[3];
rz(-1.0335575) q[3];
sx q[3];
rz(1.1972103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0672368) q[2];
sx q[2];
rz(-1.3398193) q[2];
sx q[2];
rz(-0.92607099) q[2];
rz(-0.98008424) q[3];
sx q[3];
rz(-0.96433774) q[3];
sx q[3];
rz(2.4721691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7922908) q[0];
sx q[0];
rz(-1.8632357) q[0];
sx q[0];
rz(1.5128304) q[0];
rz(-2.8577562) q[1];
sx q[1];
rz(-2.2161039) q[1];
sx q[1];
rz(1.1121174) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6215206) q[0];
sx q[0];
rz(-3.1315098) q[0];
sx q[0];
rz(0.79147379) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6792504) q[2];
sx q[2];
rz(-1.6466093) q[2];
sx q[2];
rz(-0.52559847) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6472811) q[1];
sx q[1];
rz(-2.0485281) q[1];
sx q[1];
rz(1.7812438) q[1];
x q[2];
rz(-1.2073969) q[3];
sx q[3];
rz(-1.8135241) q[3];
sx q[3];
rz(-0.12871615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5230368) q[2];
sx q[2];
rz(-1.3448389) q[2];
sx q[2];
rz(-1.1019361) q[2];
rz(2.2526422) q[3];
sx q[3];
rz(-0.44822732) q[3];
sx q[3];
rz(-2.2811208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76982826) q[0];
sx q[0];
rz(-0.17426057) q[0];
sx q[0];
rz(0.68921047) q[0];
rz(-0.31461942) q[1];
sx q[1];
rz(-0.92051053) q[1];
sx q[1];
rz(1.7291732) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1293761) q[0];
sx q[0];
rz(-2.3123992) q[0];
sx q[0];
rz(2.0022805) q[0];
rz(-pi) q[1];
rz(-2.7233549) q[2];
sx q[2];
rz(-0.23242885) q[2];
sx q[2];
rz(1.4168036) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1656944) q[1];
sx q[1];
rz(-1.1912575) q[1];
sx q[1];
rz(-1.6326153) q[1];
rz(0.19417089) q[3];
sx q[3];
rz(-1.0946858) q[3];
sx q[3];
rz(0.74060696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7249001) q[2];
sx q[2];
rz(-2.576454) q[2];
sx q[2];
rz(0.55595428) q[2];
rz(1.2712831) q[3];
sx q[3];
rz(-1.8794182) q[3];
sx q[3];
rz(-0.2909734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3274662) q[0];
sx q[0];
rz(-3.0004582) q[0];
sx q[0];
rz(-2.5471174) q[0];
rz(0.56198436) q[1];
sx q[1];
rz(-0.85834208) q[1];
sx q[1];
rz(2.1655703) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41244477) q[0];
sx q[0];
rz(-1.2093822) q[0];
sx q[0];
rz(-2.0301719) q[0];
x q[1];
rz(1.541829) q[2];
sx q[2];
rz(-2.3246362) q[2];
sx q[2];
rz(3.0532233) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4415293) q[1];
sx q[1];
rz(-2.2231327) q[1];
sx q[1];
rz(1.9487914) q[1];
x q[2];
rz(-0.30140437) q[3];
sx q[3];
rz(-1.3632953) q[3];
sx q[3];
rz(1.0338155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.48656616) q[2];
sx q[2];
rz(-1.8883294) q[2];
sx q[2];
rz(-2.5308934) q[2];
rz(0.30580172) q[3];
sx q[3];
rz(-0.96547258) q[3];
sx q[3];
rz(1.3293728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.3139451) q[0];
sx q[0];
rz(-0.018445404) q[0];
sx q[0];
rz(-1.4661283) q[0];
rz(0.77955359) q[1];
sx q[1];
rz(-1.9773217) q[1];
sx q[1];
rz(2.4028042) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8429564) q[0];
sx q[0];
rz(-1.8514575) q[0];
sx q[0];
rz(2.929351) q[0];
rz(-pi) q[1];
rz(0.03348695) q[2];
sx q[2];
rz(-0.30695633) q[2];
sx q[2];
rz(-0.54881664) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.14911095) q[1];
sx q[1];
rz(-0.68435366) q[1];
sx q[1];
rz(2.1696287) q[1];
x q[2];
rz(-1.1973279) q[3];
sx q[3];
rz(-1.9091144) q[3];
sx q[3];
rz(-0.23972971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.61591992) q[2];
sx q[2];
rz(-1.9024666) q[2];
sx q[2];
rz(-2.2789148) q[2];
rz(-2.681813) q[3];
sx q[3];
rz(-1.516187) q[3];
sx q[3];
rz(1.1427243) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2343242) q[0];
sx q[0];
rz(-0.37721226) q[0];
sx q[0];
rz(2.1211076) q[0];
rz(-0.057295784) q[1];
sx q[1];
rz(-1.658193) q[1];
sx q[1];
rz(-2.3540156) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6660992) q[0];
sx q[0];
rz(-1.3331873) q[0];
sx q[0];
rz(2.4852024) q[0];
x q[1];
rz(1.4992981) q[2];
sx q[2];
rz(-1.3132957) q[2];
sx q[2];
rz(1.6704287) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5912093) q[1];
sx q[1];
rz(-2.0616643) q[1];
sx q[1];
rz(-0.20259133) q[1];
x q[2];
rz(-1.5191139) q[3];
sx q[3];
rz(-2.5530836) q[3];
sx q[3];
rz(1.1872329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.82137498) q[2];
sx q[2];
rz(-1.3221952) q[2];
sx q[2];
rz(-2.5569432) q[2];
rz(-2.4892877) q[3];
sx q[3];
rz(-1.2290596) q[3];
sx q[3];
rz(-1.339213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(-2.898191) q[0];
sx q[0];
rz(-1.3404055) q[0];
sx q[0];
rz(-0.32824326) q[0];
rz(-0.39168656) q[1];
sx q[1];
rz(-1.990254) q[1];
sx q[1];
rz(-1.6627056) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97349629) q[0];
sx q[0];
rz(-0.35051051) q[0];
sx q[0];
rz(-2.4309733) q[0];
x q[1];
rz(-0.52430341) q[2];
sx q[2];
rz(-1.9881762) q[2];
sx q[2];
rz(-2.5427172) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0285769) q[1];
sx q[1];
rz(-1.2927107) q[1];
sx q[1];
rz(-0.098829513) q[1];
rz(2.9313179) q[3];
sx q[3];
rz(-1.9326903) q[3];
sx q[3];
rz(-1.1392322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0216003) q[2];
sx q[2];
rz(-1.3730717) q[2];
sx q[2];
rz(2.3528986) q[2];
rz(0.24104077) q[3];
sx q[3];
rz(-2.2153416) q[3];
sx q[3];
rz(2.327976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4929844) q[0];
sx q[0];
rz(-2.6032175) q[0];
sx q[0];
rz(-0.96187821) q[0];
rz(-2.9073763) q[1];
sx q[1];
rz(-2.0030256) q[1];
sx q[1];
rz(-2.403517) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3320112) q[0];
sx q[0];
rz(-1.1763078) q[0];
sx q[0];
rz(1.3138735) q[0];
rz(2.5066824) q[2];
sx q[2];
rz(-1.2037841) q[2];
sx q[2];
rz(0.5762595) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2710423) q[1];
sx q[1];
rz(-1.1525453) q[1];
sx q[1];
rz(-2.8224432) q[1];
x q[2];
rz(1.8515737) q[3];
sx q[3];
rz(-0.92377907) q[3];
sx q[3];
rz(-0.10296497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6568079) q[2];
sx q[2];
rz(-1.022555) q[2];
sx q[2];
rz(1.8812995) q[2];
rz(1.3336746) q[3];
sx q[3];
rz(-2.1402054) q[3];
sx q[3];
rz(-0.26237747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8986847) q[0];
sx q[0];
rz(-0.81166357) q[0];
sx q[0];
rz(0.86724487) q[0];
rz(2.1278837) q[1];
sx q[1];
rz(-1.8127245) q[1];
sx q[1];
rz(1.0702466) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1742221) q[0];
sx q[0];
rz(-1.4486607) q[0];
sx q[0];
rz(1.4184667) q[0];
rz(0.70906822) q[2];
sx q[2];
rz(-0.96154172) q[2];
sx q[2];
rz(2.9035062) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.874843) q[1];
sx q[1];
rz(-1.0581746) q[1];
sx q[1];
rz(-2.0155409) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6813047) q[3];
sx q[3];
rz(-1.7095057) q[3];
sx q[3];
rz(-1.8294301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9762207) q[2];
sx q[2];
rz(-0.82087159) q[2];
sx q[2];
rz(-1.2316068) q[2];
rz(1.5008789) q[3];
sx q[3];
rz(-0.21017635) q[3];
sx q[3];
rz(0.73178449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3043542) q[0];
sx q[0];
rz(-1.9914892) q[0];
sx q[0];
rz(-1.7627841) q[0];
rz(-2.7267743) q[1];
sx q[1];
rz(-1.7638313) q[1];
sx q[1];
rz(1.5324963) q[1];
rz(1.6258705) q[2];
sx q[2];
rz(-1.1255506) q[2];
sx q[2];
rz(-0.19565565) q[2];
rz(2.1555156) q[3];
sx q[3];
rz(-2.232983) q[3];
sx q[3];
rz(1.9017526) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
