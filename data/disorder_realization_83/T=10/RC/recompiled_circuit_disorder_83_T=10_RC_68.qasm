OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5873592) q[0];
sx q[0];
rz(-0.98266196) q[0];
sx q[0];
rz(2.3798556) q[0];
rz(0.76454437) q[1];
sx q[1];
rz(-2.0643056) q[1];
sx q[1];
rz(0.74365562) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3418158) q[0];
sx q[0];
rz(-1.8108978) q[0];
sx q[0];
rz(-3.0738897) q[0];
rz(0.6526297) q[2];
sx q[2];
rz(-1.4748117) q[2];
sx q[2];
rz(1.3053615) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.57063369) q[1];
sx q[1];
rz(-1.1753517) q[1];
sx q[1];
rz(2.8407437) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7573962) q[3];
sx q[3];
rz(-0.36986923) q[3];
sx q[3];
rz(-0.44934011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4347697) q[2];
sx q[2];
rz(-2.7396024) q[2];
sx q[2];
rz(-0.10786954) q[2];
rz(-0.14262959) q[3];
sx q[3];
rz(-1.7248036) q[3];
sx q[3];
rz(-0.67255783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07664872) q[0];
sx q[0];
rz(-2.3363484) q[0];
sx q[0];
rz(-2.9108677) q[0];
rz(-1.8143066) q[1];
sx q[1];
rz(-2.4704411) q[1];
sx q[1];
rz(-3.1006295) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.165867) q[0];
sx q[0];
rz(-1.8164993) q[0];
sx q[0];
rz(0.16733549) q[0];
rz(-pi) q[1];
rz(-0.83611739) q[2];
sx q[2];
rz(-0.40514075) q[2];
sx q[2];
rz(0.9678313) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4119271) q[1];
sx q[1];
rz(-1.1150868) q[1];
sx q[1];
rz(0.63974849) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9713692) q[3];
sx q[3];
rz(-1.922732) q[3];
sx q[3];
rz(0.34917253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9791947) q[2];
sx q[2];
rz(-1.4731864) q[2];
sx q[2];
rz(0.56817788) q[2];
rz(-0.5125106) q[3];
sx q[3];
rz(-2.5458702) q[3];
sx q[3];
rz(-0.56186831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7115962) q[0];
sx q[0];
rz(-1.6226409) q[0];
sx q[0];
rz(-2.6913753) q[0];
rz(-1.846116) q[1];
sx q[1];
rz(-1.1643012) q[1];
sx q[1];
rz(-0.67726642) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8836425) q[0];
sx q[0];
rz(-1.5188688) q[0];
sx q[0];
rz(-0.16846637) q[0];
rz(-pi) q[1];
rz(0.1707503) q[2];
sx q[2];
rz(-3.1351334) q[2];
sx q[2];
rz(3.1250172) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8202159) q[1];
sx q[1];
rz(-1.7019338) q[1];
sx q[1];
rz(-1.3908435) q[1];
x q[2];
rz(-0.47795313) q[3];
sx q[3];
rz(-1.9216929) q[3];
sx q[3];
rz(1.705845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.21851097) q[2];
sx q[2];
rz(-1.1542902) q[2];
sx q[2];
rz(-0.046860524) q[2];
rz(0.81165195) q[3];
sx q[3];
rz(-2.4620158) q[3];
sx q[3];
rz(2.9714382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99455225) q[0];
sx q[0];
rz(-3.0299598) q[0];
sx q[0];
rz(-3.0901093) q[0];
rz(-0.4908081) q[1];
sx q[1];
rz(-0.97508109) q[1];
sx q[1];
rz(-1.9690008) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20576142) q[0];
sx q[0];
rz(-1.7044221) q[0];
sx q[0];
rz(-0.66098102) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5355859) q[2];
sx q[2];
rz(-1.5835985) q[2];
sx q[2];
rz(-2.9738604) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6231411) q[1];
sx q[1];
rz(-1.271283) q[1];
sx q[1];
rz(-1.9407942) q[1];
rz(-0.074313642) q[3];
sx q[3];
rz(-2.001611) q[3];
sx q[3];
rz(2.5527692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2525758) q[2];
sx q[2];
rz(-0.93095195) q[2];
sx q[2];
rz(-2.6546997) q[2];
rz(1.1934818) q[3];
sx q[3];
rz(-0.32961696) q[3];
sx q[3];
rz(0.025432767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.837773) q[0];
sx q[0];
rz(-2.1210414) q[0];
sx q[0];
rz(1.3254962) q[0];
rz(0.59108132) q[1];
sx q[1];
rz(-2.4609844) q[1];
sx q[1];
rz(-2.4868734) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39010534) q[0];
sx q[0];
rz(-0.57919466) q[0];
sx q[0];
rz(0.67436995) q[0];
rz(-2.7227313) q[2];
sx q[2];
rz(-1.7423811) q[2];
sx q[2];
rz(2.3562743) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3760738) q[1];
sx q[1];
rz(-0.79080938) q[1];
sx q[1];
rz(0.54733025) q[1];
rz(-pi) q[2];
rz(0.22722865) q[3];
sx q[3];
rz(-0.39441808) q[3];
sx q[3];
rz(-1.3487181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4513662) q[2];
sx q[2];
rz(-1.2008685) q[2];
sx q[2];
rz(0.5711242) q[2];
rz(-0.58978224) q[3];
sx q[3];
rz(-2.6815806) q[3];
sx q[3];
rz(-3.0736198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1968483) q[0];
sx q[0];
rz(-1.9510883) q[0];
sx q[0];
rz(-0.29904547) q[0];
rz(1.8213182) q[1];
sx q[1];
rz(-2.8838005) q[1];
sx q[1];
rz(1.4978131) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70443557) q[0];
sx q[0];
rz(-0.91059443) q[0];
sx q[0];
rz(3.0755088) q[0];
rz(-pi) q[1];
rz(0.80070337) q[2];
sx q[2];
rz(-2.3869037) q[2];
sx q[2];
rz(0.18037361) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0648246) q[1];
sx q[1];
rz(-1.2829363) q[1];
sx q[1];
rz(2.0516146) q[1];
rz(-pi) q[2];
rz(-2.8619814) q[3];
sx q[3];
rz(-1.8108484) q[3];
sx q[3];
rz(-0.9595426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6254639) q[2];
sx q[2];
rz(-1.4279782) q[2];
sx q[2];
rz(-3.1141172) q[2];
rz(0.52250683) q[3];
sx q[3];
rz(-0.80111879) q[3];
sx q[3];
rz(2.3251422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7291173) q[0];
sx q[0];
rz(-1.1183879) q[0];
sx q[0];
rz(0.11418848) q[0];
rz(-2.1633637) q[1];
sx q[1];
rz(-0.54662919) q[1];
sx q[1];
rz(2.3506929) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0349802) q[0];
sx q[0];
rz(-2.5133586) q[0];
sx q[0];
rz(-1.4094704) q[0];
rz(-pi) q[1];
rz(0.74366624) q[2];
sx q[2];
rz(-1.3197834) q[2];
sx q[2];
rz(2.2199092) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0485718) q[1];
sx q[1];
rz(-1.2067817) q[1];
sx q[1];
rz(-0.084333468) q[1];
rz(-0.88364717) q[3];
sx q[3];
rz(-1.8568608) q[3];
sx q[3];
rz(1.4253695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.87166446) q[2];
sx q[2];
rz(-2.0285138) q[2];
sx q[2];
rz(0.86501914) q[2];
rz(-2.4690752) q[3];
sx q[3];
rz(-1.1285684) q[3];
sx q[3];
rz(3.045936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(1.6487811) q[0];
sx q[0];
rz(-2.6219941) q[0];
sx q[0];
rz(0.46736026) q[0];
rz(-2.6043747) q[1];
sx q[1];
rz(-0.98058128) q[1];
sx q[1];
rz(0.25407243) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2652889) q[0];
sx q[0];
rz(-0.2067925) q[0];
sx q[0];
rz(-0.32949038) q[0];
rz(-pi) q[1];
x q[1];
rz(1.569283) q[2];
sx q[2];
rz(-1.5655893) q[2];
sx q[2];
rz(-1.4323915) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4916363) q[1];
sx q[1];
rz(-1.2600113) q[1];
sx q[1];
rz(-1.4375356) q[1];
rz(-3.0244175) q[3];
sx q[3];
rz(-2.5391038) q[3];
sx q[3];
rz(2.6206827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6415928) q[2];
sx q[2];
rz(-2.3708512) q[2];
sx q[2];
rz(0.33561486) q[2];
rz(2.8149758) q[3];
sx q[3];
rz(-0.89544046) q[3];
sx q[3];
rz(-1.7232822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0582054) q[0];
sx q[0];
rz(-2.931262) q[0];
sx q[0];
rz(2.0876419) q[0];
rz(0.15696934) q[1];
sx q[1];
rz(-1.7228246) q[1];
sx q[1];
rz(-2.1597247) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40030038) q[0];
sx q[0];
rz(-1.9809082) q[0];
sx q[0];
rz(-1.4492695) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0018714) q[2];
sx q[2];
rz(-1.47942) q[2];
sx q[2];
rz(0.56359529) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.30584221) q[1];
sx q[1];
rz(-1.7562477) q[1];
sx q[1];
rz(-3.0830543) q[1];
rz(-0.710886) q[3];
sx q[3];
rz(-1.3804187) q[3];
sx q[3];
rz(-0.058660942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92423576) q[2];
sx q[2];
rz(-1.057426) q[2];
sx q[2];
rz(2.1441933) q[2];
rz(-0.50619566) q[3];
sx q[3];
rz(-0.96118569) q[3];
sx q[3];
rz(-3.1072646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85703325) q[0];
sx q[0];
rz(-0.71145809) q[0];
sx q[0];
rz(0.51666623) q[0];
rz(2.754028) q[1];
sx q[1];
rz(-2.0811847) q[1];
sx q[1];
rz(-2.8732252) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70290138) q[0];
sx q[0];
rz(-0.99150204) q[0];
sx q[0];
rz(0.86433522) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1860113) q[2];
sx q[2];
rz(-0.68253839) q[2];
sx q[2];
rz(-0.62894422) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8908773) q[1];
sx q[1];
rz(-2.4836575) q[1];
sx q[1];
rz(-2.6943745) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9645421) q[3];
sx q[3];
rz(-2.5265794) q[3];
sx q[3];
rz(2.7494489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9364075) q[2];
sx q[2];
rz(-0.78339094) q[2];
sx q[2];
rz(0.56383413) q[2];
rz(1.9514203) q[3];
sx q[3];
rz(-2.1824013) q[3];
sx q[3];
rz(-2.4998375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.664809) q[0];
sx q[0];
rz(-1.8905147) q[0];
sx q[0];
rz(2.074194) q[0];
rz(-1.339636) q[1];
sx q[1];
rz(-1.7032774) q[1];
sx q[1];
rz(1.3443321) q[1];
rz(0.50189353) q[2];
sx q[2];
rz(-2.5645651) q[2];
sx q[2];
rz(-1.2679451) q[2];
rz(2.8292538) q[3];
sx q[3];
rz(-1.3795508) q[3];
sx q[3];
rz(3.0637904) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];