OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7497082) q[0];
sx q[0];
rz(-2.9449129) q[0];
sx q[0];
rz(-1.1893907) q[0];
rz(0.2285129) q[1];
sx q[1];
rz(-0.84140468) q[1];
sx q[1];
rz(-2.7639311) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084378622) q[0];
sx q[0];
rz(-2.9998261) q[0];
sx q[0];
rz(-2.111582) q[0];
x q[1];
rz(-1.471465) q[2];
sx q[2];
rz(-0.28495312) q[2];
sx q[2];
rz(-0.003665912) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9610112) q[1];
sx q[1];
rz(-2.2716224) q[1];
sx q[1];
rz(-0.58971528) q[1];
rz(-pi) q[2];
rz(-1.7114867) q[3];
sx q[3];
rz(-1.1555539) q[3];
sx q[3];
rz(-1.5822441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1315786) q[2];
sx q[2];
rz(-2.6476314) q[2];
sx q[2];
rz(0.67260355) q[2];
rz(-0.16942313) q[3];
sx q[3];
rz(-2.7526581) q[3];
sx q[3];
rz(1.3385564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.910903) q[0];
sx q[0];
rz(-0.90086532) q[0];
sx q[0];
rz(0.22856523) q[0];
rz(2.9810492) q[1];
sx q[1];
rz(-1.4385782) q[1];
sx q[1];
rz(0.28796089) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32591336) q[0];
sx q[0];
rz(-2.8796112) q[0];
sx q[0];
rz(2.322305) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38950133) q[2];
sx q[2];
rz(-2.1283538) q[2];
sx q[2];
rz(-2.8745289) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.49111734) q[1];
sx q[1];
rz(-0.75956356) q[1];
sx q[1];
rz(2.526545) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.037022) q[3];
sx q[3];
rz(-0.35211709) q[3];
sx q[3];
rz(-0.62192384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.59445375) q[2];
sx q[2];
rz(-1.2524266) q[2];
sx q[2];
rz(-1.4734369) q[2];
rz(2.2041221) q[3];
sx q[3];
rz(-2.6963186) q[3];
sx q[3];
rz(-0.37500769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32524747) q[0];
sx q[0];
rz(-0.49961093) q[0];
sx q[0];
rz(0.26741272) q[0];
rz(-1.7193517) q[1];
sx q[1];
rz(-1.1317252) q[1];
sx q[1];
rz(0.95169383) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0824453) q[0];
sx q[0];
rz(-1.7623616) q[0];
sx q[0];
rz(3.0491769) q[0];
x q[1];
rz(2.0864262) q[2];
sx q[2];
rz(-1.7679169) q[2];
sx q[2];
rz(-2.4899763) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0905076) q[1];
sx q[1];
rz(-1.2625492) q[1];
sx q[1];
rz(-0.80866637) q[1];
x q[2];
rz(2.2610407) q[3];
sx q[3];
rz(-2.3909702) q[3];
sx q[3];
rz(0.39117884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0114228) q[2];
sx q[2];
rz(-1.495785) q[2];
sx q[2];
rz(2.7974131) q[2];
rz(2.2551645) q[3];
sx q[3];
rz(-2.8635946) q[3];
sx q[3];
rz(-2.5755431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13720559) q[0];
sx q[0];
rz(-1.6442278) q[0];
sx q[0];
rz(-1.244506) q[0];
rz(-3.0124774) q[1];
sx q[1];
rz(-1.3543509) q[1];
sx q[1];
rz(2.7688162) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6308206) q[0];
sx q[0];
rz(-1.6850867) q[0];
sx q[0];
rz(-0.75829102) q[0];
x q[1];
rz(-1.1431085) q[2];
sx q[2];
rz(-2.3017985) q[2];
sx q[2];
rz(2.1172303) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8744295) q[1];
sx q[1];
rz(-2.0544555) q[1];
sx q[1];
rz(2.8664385) q[1];
rz(-pi) q[2];
rz(2.8126206) q[3];
sx q[3];
rz(-1.1165459) q[3];
sx q[3];
rz(-2.8153552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.63561511) q[2];
sx q[2];
rz(-1.4900692) q[2];
sx q[2];
rz(-0.90488952) q[2];
rz(-2.7010226) q[3];
sx q[3];
rz(-2.6456656) q[3];
sx q[3];
rz(-2.1499965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47914094) q[0];
sx q[0];
rz(-0.15563706) q[0];
sx q[0];
rz(1.2878081) q[0];
rz(-1.4783391) q[1];
sx q[1];
rz(-1.1454502) q[1];
sx q[1];
rz(0.53422654) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14311612) q[0];
sx q[0];
rz(-1.2623708) q[0];
sx q[0];
rz(0.31750676) q[0];
rz(-pi) q[1];
rz(-2.1557501) q[2];
sx q[2];
rz(-0.57136977) q[2];
sx q[2];
rz(-2.5845598) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.91671645) q[1];
sx q[1];
rz(-1.4955048) q[1];
sx q[1];
rz(3.0776943) q[1];
rz(-pi) q[2];
rz(-0.71877919) q[3];
sx q[3];
rz(-1.6703509) q[3];
sx q[3];
rz(2.0841141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5010895) q[2];
sx q[2];
rz(-1.9260294) q[2];
sx q[2];
rz(-0.14349288) q[2];
rz(1.3714553) q[3];
sx q[3];
rz(-1.2797132) q[3];
sx q[3];
rz(2.9218856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4390398) q[0];
sx q[0];
rz(-2.2655903) q[0];
sx q[0];
rz(-2.904536) q[0];
rz(-1.9006231) q[1];
sx q[1];
rz(-2.3126912) q[1];
sx q[1];
rz(1.7664849) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9532167) q[0];
sx q[0];
rz(-1.2326476) q[0];
sx q[0];
rz(1.2921278) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0639406) q[2];
sx q[2];
rz(-0.75040557) q[2];
sx q[2];
rz(2.1466308) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3875229) q[1];
sx q[1];
rz(-1.9807528) q[1];
sx q[1];
rz(-1.5974664) q[1];
rz(-2.6508413) q[3];
sx q[3];
rz(-1.1928344) q[3];
sx q[3];
rz(2.2839387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6087626) q[2];
sx q[2];
rz(-2.4217748) q[2];
sx q[2];
rz(2.9928845) q[2];
rz(0.016629774) q[3];
sx q[3];
rz(-2.7745268) q[3];
sx q[3];
rz(2.9490411) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6454813) q[0];
sx q[0];
rz(-2.8493024) q[0];
sx q[0];
rz(-0.87316978) q[0];
rz(-0.9219777) q[1];
sx q[1];
rz(-1.1154122) q[1];
sx q[1];
rz(-3.057664) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8221995) q[0];
sx q[0];
rz(-2.5376352) q[0];
sx q[0];
rz(-1.970406) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7883045) q[2];
sx q[2];
rz(-1.7751667) q[2];
sx q[2];
rz(0.7962966) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.88372701) q[1];
sx q[1];
rz(-0.95974937) q[1];
sx q[1];
rz(2.3962767) q[1];
rz(-pi) q[2];
x q[2];
rz(0.069025741) q[3];
sx q[3];
rz(-0.8419753) q[3];
sx q[3];
rz(2.98711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.38368791) q[2];
sx q[2];
rz(-2.7010475) q[2];
sx q[2];
rz(2.596358) q[2];
rz(0.43045726) q[3];
sx q[3];
rz(-2.0038219) q[3];
sx q[3];
rz(-2.8366413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51171821) q[0];
sx q[0];
rz(-3.1261303) q[0];
sx q[0];
rz(-1.9301201) q[0];
rz(0.97310549) q[1];
sx q[1];
rz(-2.548023) q[1];
sx q[1];
rz(-2.1957695) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57293939) q[0];
sx q[0];
rz(-1.6488254) q[0];
sx q[0];
rz(-0.22139876) q[0];
x q[1];
rz(0.85407599) q[2];
sx q[2];
rz(-1.0208566) q[2];
sx q[2];
rz(-2.9727109) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9981873) q[1];
sx q[1];
rz(-0.37969509) q[1];
sx q[1];
rz(-0.077296301) q[1];
rz(-pi) q[2];
rz(2.4754727) q[3];
sx q[3];
rz(-1.5770116) q[3];
sx q[3];
rz(2.4339649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2716081) q[2];
sx q[2];
rz(-1.4166069) q[2];
sx q[2];
rz(2.8611709) q[2];
rz(-2.5366606) q[3];
sx q[3];
rz(-1.0792024) q[3];
sx q[3];
rz(-2.159507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.1542926) q[0];
sx q[0];
rz(-2.2315318) q[0];
sx q[0];
rz(-0.61532414) q[0];
rz(0.92957169) q[1];
sx q[1];
rz(-2.2832182) q[1];
sx q[1];
rz(-2.5659134) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070411365) q[0];
sx q[0];
rz(-2.1343263) q[0];
sx q[0];
rz(1.0067183) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3747146) q[2];
sx q[2];
rz(-1.6819281) q[2];
sx q[2];
rz(2.9642504) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.8262144) q[1];
sx q[1];
rz(-1.678102) q[1];
sx q[1];
rz(3.1158434) q[1];
rz(-pi) q[2];
rz(-1.0300893) q[3];
sx q[3];
rz(-1.4095777) q[3];
sx q[3];
rz(-2.2459523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3118887) q[2];
sx q[2];
rz(-0.76278937) q[2];
sx q[2];
rz(0.021961948) q[2];
rz(2.9711376) q[3];
sx q[3];
rz(-2.1178092) q[3];
sx q[3];
rz(-0.26486614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72572529) q[0];
sx q[0];
rz(-0.9265582) q[0];
sx q[0];
rz(-0.6341933) q[0];
rz(0.028907396) q[1];
sx q[1];
rz(-0.78866619) q[1];
sx q[1];
rz(0.96910563) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.076981912) q[0];
sx q[0];
rz(-1.4990028) q[0];
sx q[0];
rz(-3.0781156) q[0];
x q[1];
rz(-1.3766039) q[2];
sx q[2];
rz(-2.3079434) q[2];
sx q[2];
rz(-0.091723524) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9173968) q[1];
sx q[1];
rz(-0.90985137) q[1];
sx q[1];
rz(-1.3437273) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7006847) q[3];
sx q[3];
rz(-1.8280067) q[3];
sx q[3];
rz(-0.5514901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6897631) q[2];
sx q[2];
rz(-1.4844866) q[2];
sx q[2];
rz(0.36007145) q[2];
rz(-1.6137971) q[3];
sx q[3];
rz(-2.4751622) q[3];
sx q[3];
rz(0.56267363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2071028) q[0];
sx q[0];
rz(-1.5705382) q[0];
sx q[0];
rz(-1.6194153) q[0];
rz(3.0974401) q[1];
sx q[1];
rz(-1.4587198) q[1];
sx q[1];
rz(-1.1062467) q[1];
rz(0.85514851) q[2];
sx q[2];
rz(-2.6519041) q[2];
sx q[2];
rz(2.3408163) q[2];
rz(-2.5526657) q[3];
sx q[3];
rz(-0.52994655) q[3];
sx q[3];
rz(0.51701057) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];