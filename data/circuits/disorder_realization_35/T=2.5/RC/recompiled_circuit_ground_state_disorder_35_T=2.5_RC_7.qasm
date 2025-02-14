OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2090981) q[0];
sx q[0];
rz(-3.0734835) q[0];
sx q[0];
rz(-0.77686247) q[0];
rz(1.8344954) q[1];
sx q[1];
rz(-2.1721462) q[1];
sx q[1];
rz(-1.3219272) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37047077) q[0];
sx q[0];
rz(-1.7092472) q[0];
sx q[0];
rz(-2.2884503) q[0];
rz(-pi) q[1];
rz(-2.5547682) q[2];
sx q[2];
rz(-0.57673645) q[2];
sx q[2];
rz(-1.5962034) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.81404954) q[1];
sx q[1];
rz(-2.7316878) q[1];
sx q[1];
rz(-1.4973212) q[1];
x q[2];
rz(-2.9453927) q[3];
sx q[3];
rz(-1.2854738) q[3];
sx q[3];
rz(1.6183492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0778568) q[2];
sx q[2];
rz(-1.7813762) q[2];
sx q[2];
rz(-2.784101) q[2];
rz(0.56719559) q[3];
sx q[3];
rz(-1.4128069) q[3];
sx q[3];
rz(2.7020057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1883063) q[0];
sx q[0];
rz(-2.8585241) q[0];
sx q[0];
rz(1.8081007) q[0];
rz(-2.7045344) q[1];
sx q[1];
rz(-0.60085618) q[1];
sx q[1];
rz(2.3016047) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6463722) q[0];
sx q[0];
rz(-2.6073808) q[0];
sx q[0];
rz(-1.87508) q[0];
rz(0.2510906) q[2];
sx q[2];
rz(-0.60949627) q[2];
sx q[2];
rz(1.9025546) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4639791) q[1];
sx q[1];
rz(-2.2871454) q[1];
sx q[1];
rz(-2.3542627) q[1];
x q[2];
rz(2.1206109) q[3];
sx q[3];
rz(-1.6040168) q[3];
sx q[3];
rz(2.0149505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7175032) q[2];
sx q[2];
rz(-1.4428029) q[2];
sx q[2];
rz(-1.7355512) q[2];
rz(0.16215912) q[3];
sx q[3];
rz(-1.5608965) q[3];
sx q[3];
rz(2.6774008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22299396) q[0];
sx q[0];
rz(-3.0610237) q[0];
sx q[0];
rz(2.719847) q[0];
rz(-0.84284198) q[1];
sx q[1];
rz(-1.3858567) q[1];
sx q[1];
rz(-2.8657894) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11711794) q[0];
sx q[0];
rz(-1.8016651) q[0];
sx q[0];
rz(-0.33322115) q[0];
rz(-pi) q[1];
rz(-2.7245443) q[2];
sx q[2];
rz(-2.6337503) q[2];
sx q[2];
rz(-1.5323973) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.35728243) q[1];
sx q[1];
rz(-1.6228592) q[1];
sx q[1];
rz(2.8207458) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6869322) q[3];
sx q[3];
rz(-0.7051055) q[3];
sx q[3];
rz(3.1381315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.80321035) q[2];
sx q[2];
rz(-2.3440177) q[2];
sx q[2];
rz(-0.69022834) q[2];
rz(-2.7817182) q[3];
sx q[3];
rz(-1.7706324) q[3];
sx q[3];
rz(2.7580822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91367078) q[0];
sx q[0];
rz(-2.0991195) q[0];
sx q[0];
rz(-2.2699455) q[0];
rz(-2.7323885) q[1];
sx q[1];
rz(-0.49574655) q[1];
sx q[1];
rz(-2.8292378) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6808514) q[0];
sx q[0];
rz(-0.4609403) q[0];
sx q[0];
rz(-0.26682202) q[0];
rz(1.4710452) q[2];
sx q[2];
rz(-1.3653697) q[2];
sx q[2];
rz(-0.86459407) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.40939366) q[1];
sx q[1];
rz(-2.7229561) q[1];
sx q[1];
rz(-2.6217209) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4931025) q[3];
sx q[3];
rz(-1.7269508) q[3];
sx q[3];
rz(0.79532571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.006134) q[2];
sx q[2];
rz(-1.4486855) q[2];
sx q[2];
rz(0.82497605) q[2];
rz(-1.4247591) q[3];
sx q[3];
rz(-0.32130876) q[3];
sx q[3];
rz(-2.7062374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7216126) q[0];
sx q[0];
rz(-1.551832) q[0];
sx q[0];
rz(-2.2671674) q[0];
rz(-2.8127316) q[1];
sx q[1];
rz(-1.3855653) q[1];
sx q[1];
rz(2.0430476) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1060396) q[0];
sx q[0];
rz(-1.6811724) q[0];
sx q[0];
rz(-2.7644185) q[0];
rz(-pi) q[1];
rz(2.3057865) q[2];
sx q[2];
rz(-1.5440189) q[2];
sx q[2];
rz(-1.4876897) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9096127) q[1];
sx q[1];
rz(-1.8706338) q[1];
sx q[1];
rz(1.058387) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9727398) q[3];
sx q[3];
rz(-0.62905772) q[3];
sx q[3];
rz(2.9161616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.141779) q[2];
sx q[2];
rz(-1.9580611) q[2];
sx q[2];
rz(0.67406526) q[2];
rz(1.6541803) q[3];
sx q[3];
rz(-1.6964648) q[3];
sx q[3];
rz(-2.8580247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0312626) q[0];
sx q[0];
rz(-1.0227579) q[0];
sx q[0];
rz(-1.7559825) q[0];
rz(-1.1194057) q[1];
sx q[1];
rz(-1.4675354) q[1];
sx q[1];
rz(-1.3853692) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.114217) q[0];
sx q[0];
rz(-1.3189335) q[0];
sx q[0];
rz(-2.2172539) q[0];
x q[1];
rz(0.49853168) q[2];
sx q[2];
rz(-1.4953002) q[2];
sx q[2];
rz(-0.99270644) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.54810235) q[1];
sx q[1];
rz(-1.2498651) q[1];
sx q[1];
rz(-1.3671419) q[1];
x q[2];
rz(-2.9117111) q[3];
sx q[3];
rz(-0.15860367) q[3];
sx q[3];
rz(1.5029217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0783656) q[2];
sx q[2];
rz(-2.207022) q[2];
sx q[2];
rz(-0.37609491) q[2];
rz(-1.1345351) q[3];
sx q[3];
rz(-2.7781656) q[3];
sx q[3];
rz(-2.334972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.019808708) q[0];
sx q[0];
rz(-0.60096318) q[0];
sx q[0];
rz(3.0390749) q[0];
rz(0.58492297) q[1];
sx q[1];
rz(-0.94766098) q[1];
sx q[1];
rz(-1.1835416) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4485156) q[0];
sx q[0];
rz(-1.7127556) q[0];
sx q[0];
rz(-1.4407115) q[0];
rz(-pi) q[1];
rz(0.66254424) q[2];
sx q[2];
rz(-1.8383611) q[2];
sx q[2];
rz(-0.91987687) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.21321024) q[1];
sx q[1];
rz(-2.2021535) q[1];
sx q[1];
rz(0.31698314) q[1];
rz(2.6449219) q[3];
sx q[3];
rz(-0.52049812) q[3];
sx q[3];
rz(-0.48540533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3465603) q[2];
sx q[2];
rz(-2.2117895) q[2];
sx q[2];
rz(2.9417876) q[2];
rz(-2.0810769) q[3];
sx q[3];
rz(-2.6001055) q[3];
sx q[3];
rz(3.0919891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26315966) q[0];
sx q[0];
rz(-1.4819772) q[0];
sx q[0];
rz(0.31295452) q[0];
rz(-2.1921659) q[1];
sx q[1];
rz(-1.7890309) q[1];
sx q[1];
rz(1.688028) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87438238) q[0];
sx q[0];
rz(-1.6613164) q[0];
sx q[0];
rz(2.9541914) q[0];
x q[1];
rz(3.0309148) q[2];
sx q[2];
rz(-0.62244906) q[2];
sx q[2];
rz(1.3712997) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6387299) q[1];
sx q[1];
rz(-2.4045378) q[1];
sx q[1];
rz(0.33691892) q[1];
x q[2];
rz(1.9096776) q[3];
sx q[3];
rz(-1.1383071) q[3];
sx q[3];
rz(-0.56604715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7877385) q[2];
sx q[2];
rz(-2.165803) q[2];
sx q[2];
rz(1.806951) q[2];
rz(1.8169962) q[3];
sx q[3];
rz(-0.66671222) q[3];
sx q[3];
rz(-1.9273531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1052833) q[0];
sx q[0];
rz(-2.5258625) q[0];
sx q[0];
rz(2.4435254) q[0];
rz(2.2659194) q[1];
sx q[1];
rz(-2.0798637) q[1];
sx q[1];
rz(-0.36044136) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35931906) q[0];
sx q[0];
rz(-0.8018612) q[0];
sx q[0];
rz(-2.5129621) q[0];
rz(-pi) q[1];
rz(0.38273229) q[2];
sx q[2];
rz(-1.5110516) q[2];
sx q[2];
rz(-1.1381696) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1680582) q[1];
sx q[1];
rz(-2.413001) q[1];
sx q[1];
rz(0.97961564) q[1];
x q[2];
rz(0.54887532) q[3];
sx q[3];
rz(-1.4272808) q[3];
sx q[3];
rz(1.9560069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2299819) q[2];
sx q[2];
rz(-1.7155827) q[2];
sx q[2];
rz(-2.3243375) q[2];
rz(-2.3115555) q[3];
sx q[3];
rz(-2.6006112) q[3];
sx q[3];
rz(-2.9221007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63477409) q[0];
sx q[0];
rz(-2.2932678) q[0];
sx q[0];
rz(-0.032055227) q[0];
rz(-1.3196779) q[1];
sx q[1];
rz(-1.3751043) q[1];
sx q[1];
rz(0.65418902) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7830156) q[0];
sx q[0];
rz(-0.23046432) q[0];
sx q[0];
rz(-2.5925267) q[0];
rz(-0.78829576) q[2];
sx q[2];
rz(-1.9018146) q[2];
sx q[2];
rz(0.85110474) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4465877) q[1];
sx q[1];
rz(-0.69952974) q[1];
sx q[1];
rz(-2.3746731) q[1];
x q[2];
rz(-1.1919674) q[3];
sx q[3];
rz(-2.4740268) q[3];
sx q[3];
rz(2.6258056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0111982) q[2];
sx q[2];
rz(-1.0426714) q[2];
sx q[2];
rz(-2.1350258) q[2];
rz(-2.448163) q[3];
sx q[3];
rz(-1.8208241) q[3];
sx q[3];
rz(-1.2815732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7158159) q[0];
sx q[0];
rz(-0.67831138) q[0];
sx q[0];
rz(-0.074445733) q[0];
rz(-2.2609932) q[1];
sx q[1];
rz(-2.0279299) q[1];
sx q[1];
rz(0.50150064) q[1];
rz(0.975561) q[2];
sx q[2];
rz(-0.79339334) q[2];
sx q[2];
rz(-1.7421772) q[2];
rz(1.952259) q[3];
sx q[3];
rz(-2.2305924) q[3];
sx q[3];
rz(2.3371405) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
