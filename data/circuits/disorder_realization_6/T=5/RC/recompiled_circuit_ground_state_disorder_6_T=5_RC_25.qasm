OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5175944) q[0];
sx q[0];
rz(-2.6156293) q[0];
sx q[0];
rz(2.9232803) q[0];
rz(-1.6649618) q[1];
sx q[1];
rz(-0.47774878) q[1];
sx q[1];
rz(-2.5876317) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71611878) q[0];
sx q[0];
rz(-2.0332575) q[0];
sx q[0];
rz(1.0008971) q[0];
rz(0.10138114) q[2];
sx q[2];
rz(-0.53115618) q[2];
sx q[2];
rz(0.44975933) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7377555) q[1];
sx q[1];
rz(-2.3970876) q[1];
sx q[1];
rz(1.2584524) q[1];
rz(-1.1470818) q[3];
sx q[3];
rz(-2.0689575) q[3];
sx q[3];
rz(-0.26396449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7654984) q[2];
sx q[2];
rz(-2.8480397) q[2];
sx q[2];
rz(1.8739088) q[2];
rz(-2.3119161) q[3];
sx q[3];
rz(-1.5209578) q[3];
sx q[3];
rz(-2.7177496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0881398) q[0];
sx q[0];
rz(-1.8064878) q[0];
sx q[0];
rz(1.394519) q[0];
rz(0.8949737) q[1];
sx q[1];
rz(-1.0286237) q[1];
sx q[1];
rz(-1.5825533) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5916067) q[0];
sx q[0];
rz(-0.9261407) q[0];
sx q[0];
rz(-2.2695973) q[0];
rz(-1.1505914) q[2];
sx q[2];
rz(-1.24345) q[2];
sx q[2];
rz(-0.98132747) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.30732111) q[1];
sx q[1];
rz(-0.91047344) q[1];
sx q[1];
rz(1.1701533) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8361735) q[3];
sx q[3];
rz(-1.9807743) q[3];
sx q[3];
rz(-0.78193203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6643657) q[2];
sx q[2];
rz(-1.5218647) q[2];
sx q[2];
rz(-3.1108372) q[2];
rz(-0.45298806) q[3];
sx q[3];
rz(-2.9121297) q[3];
sx q[3];
rz(2.9636813) q[3];
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
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7396616) q[0];
sx q[0];
rz(-3.0406096) q[0];
sx q[0];
rz(2.3133551) q[0];
rz(3.0939057) q[1];
sx q[1];
rz(-2.2779155) q[1];
sx q[1];
rz(-1.2275068) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9102893) q[0];
sx q[0];
rz(-2.4132015) q[0];
sx q[0];
rz(-0.69407082) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1647878) q[2];
sx q[2];
rz(-1.3719402) q[2];
sx q[2];
rz(0.59954294) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3674254) q[1];
sx q[1];
rz(-1.8313421) q[1];
sx q[1];
rz(-1.114964) q[1];
rz(-pi) q[2];
rz(2.2994386) q[3];
sx q[3];
rz(-1.531736) q[3];
sx q[3];
rz(2.0741163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.092827169) q[2];
sx q[2];
rz(-2.2266677) q[2];
sx q[2];
rz(1.1009334) q[2];
rz(0.01384211) q[3];
sx q[3];
rz(-1.775454) q[3];
sx q[3];
rz(2.3069416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58532995) q[0];
sx q[0];
rz(-1.4034554) q[0];
sx q[0];
rz(0.334326) q[0];
rz(-2.4090134) q[1];
sx q[1];
rz(-2.1926447) q[1];
sx q[1];
rz(-3.0308731) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80522711) q[0];
sx q[0];
rz(-2.3799161) q[0];
sx q[0];
rz(0.021999981) q[0];
rz(0.51038536) q[2];
sx q[2];
rz(-2.9523627) q[2];
sx q[2];
rz(-0.58923474) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9690937) q[1];
sx q[1];
rz(-2.2374472) q[1];
sx q[1];
rz(0.73760017) q[1];
rz(-pi) q[2];
x q[2];
rz(0.48583123) q[3];
sx q[3];
rz(-2.3049092) q[3];
sx q[3];
rz(-2.776718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.683814) q[2];
sx q[2];
rz(-2.8595698) q[2];
sx q[2];
rz(-1.4078183) q[2];
rz(-1.9862566) q[3];
sx q[3];
rz(-1.1563533) q[3];
sx q[3];
rz(-0.19481625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4487576) q[0];
sx q[0];
rz(-0.91788569) q[0];
sx q[0];
rz(0.99739972) q[0];
rz(1.5265441) q[1];
sx q[1];
rz(-2.5020182) q[1];
sx q[1];
rz(1.4195199) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90516312) q[0];
sx q[0];
rz(-1.4078724) q[0];
sx q[0];
rz(1.7927367) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9827955) q[2];
sx q[2];
rz(-2.3613075) q[2];
sx q[2];
rz(-0.71035383) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0816272) q[1];
sx q[1];
rz(-0.94725376) q[1];
sx q[1];
rz(-2.7452031) q[1];
x q[2];
rz(-1.6035242) q[3];
sx q[3];
rz(-1.447289) q[3];
sx q[3];
rz(1.1203958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2898966) q[2];
sx q[2];
rz(-0.09446129) q[2];
sx q[2];
rz(2.8186901) q[2];
rz(2.0276535) q[3];
sx q[3];
rz(-2.0506004) q[3];
sx q[3];
rz(2.7285301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36393976) q[0];
sx q[0];
rz(-0.33828619) q[0];
sx q[0];
rz(0.80663484) q[0];
rz(0.58397645) q[1];
sx q[1];
rz(-1.1176502) q[1];
sx q[1];
rz(1.1899828) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7936913) q[0];
sx q[0];
rz(-2.7712203) q[0];
sx q[0];
rz(1.9166975) q[0];
rz(-0.34752589) q[2];
sx q[2];
rz(-0.56171562) q[2];
sx q[2];
rz(2.7507741) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2417702) q[1];
sx q[1];
rz(-0.17647753) q[1];
sx q[1];
rz(1.2736257) q[1];
x q[2];
rz(0.1889008) q[3];
sx q[3];
rz(-2.402225) q[3];
sx q[3];
rz(2.1738571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.40818647) q[2];
sx q[2];
rz(-1.6383645) q[2];
sx q[2];
rz(0.1850941) q[2];
rz(-1.5271651) q[3];
sx q[3];
rz(-1.7388758) q[3];
sx q[3];
rz(-2.4534498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9516893) q[0];
sx q[0];
rz(-0.95174319) q[0];
sx q[0];
rz(2.9290747) q[0];
rz(2.3566133) q[1];
sx q[1];
rz(-1.4947944) q[1];
sx q[1];
rz(0.77883887) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76701421) q[0];
sx q[0];
rz(-0.85142577) q[0];
sx q[0];
rz(-2.3423751) q[0];
rz(2.6109527) q[2];
sx q[2];
rz(-0.65766774) q[2];
sx q[2];
rz(2.4462157) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6952371) q[1];
sx q[1];
rz(-1.8255705) q[1];
sx q[1];
rz(-2.9181515) q[1];
x q[2];
rz(-0.8605729) q[3];
sx q[3];
rz(-1.1906151) q[3];
sx q[3];
rz(0.036503867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.51398858) q[2];
sx q[2];
rz(-2.6476761) q[2];
sx q[2];
rz(0.40840515) q[2];
rz(0.70185316) q[3];
sx q[3];
rz(-2.2097094) q[3];
sx q[3];
rz(1.2592038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34847611) q[0];
sx q[0];
rz(-1.9261253) q[0];
sx q[0];
rz(3.075573) q[0];
rz(-1.5090212) q[1];
sx q[1];
rz(-2.0965818) q[1];
sx q[1];
rz(-2.1836233) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1692266) q[0];
sx q[0];
rz(-1.1145076) q[0];
sx q[0];
rz(-1.2537987) q[0];
rz(0.81664576) q[2];
sx q[2];
rz(-1.1743675) q[2];
sx q[2];
rz(0.81260896) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9878475) q[1];
sx q[1];
rz(-1.2586795) q[1];
sx q[1];
rz(0.69674833) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.417172) q[3];
sx q[3];
rz(-1.6699176) q[3];
sx q[3];
rz(-2.2528668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8336739) q[2];
sx q[2];
rz(-1.5235528) q[2];
sx q[2];
rz(-1.7306805) q[2];
rz(2.446567) q[3];
sx q[3];
rz(-1.4796939) q[3];
sx q[3];
rz(0.01785774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0628292) q[0];
sx q[0];
rz(-1.48209) q[0];
sx q[0];
rz(0.98989809) q[0];
rz(-0.46317378) q[1];
sx q[1];
rz(-1.4521234) q[1];
sx q[1];
rz(1.90082) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2808229) q[0];
sx q[0];
rz(-1.6284571) q[0];
sx q[0];
rz(-1.2755544) q[0];
rz(0.2026338) q[2];
sx q[2];
rz(-1.1362193) q[2];
sx q[2];
rz(0.045534924) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.74490101) q[1];
sx q[1];
rz(-0.38977515) q[1];
sx q[1];
rz(0.74662965) q[1];
x q[2];
rz(2.1677446) q[3];
sx q[3];
rz(-0.71248369) q[3];
sx q[3];
rz(-1.9402011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3339633) q[2];
sx q[2];
rz(-1.7619851) q[2];
sx q[2];
rz(-2.4228952) q[2];
rz(1.9888318) q[3];
sx q[3];
rz(-1.4216239) q[3];
sx q[3];
rz(-2.1271472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5929247) q[0];
sx q[0];
rz(-2.7935226) q[0];
sx q[0];
rz(0.80192178) q[0];
rz(2.0536664) q[1];
sx q[1];
rz(-1.3366924) q[1];
sx q[1];
rz(0.71802872) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1205667) q[0];
sx q[0];
rz(-0.88827288) q[0];
sx q[0];
rz(1.1521856) q[0];
x q[1];
rz(-2.3921591) q[2];
sx q[2];
rz(-0.66808703) q[2];
sx q[2];
rz(-1.6783448) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7321328) q[1];
sx q[1];
rz(-1.902305) q[1];
sx q[1];
rz(-3.0863239) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78022782) q[3];
sx q[3];
rz(-2.2426668) q[3];
sx q[3];
rz(1.5742009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3512909) q[2];
sx q[2];
rz(-0.21778926) q[2];
sx q[2];
rz(0.51266074) q[2];
rz(0.74448186) q[3];
sx q[3];
rz(-1.0267886) q[3];
sx q[3];
rz(2.3775533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22513334) q[0];
sx q[0];
rz(-1.8712578) q[0];
sx q[0];
rz(1.901392) q[0];
rz(-2.4254639) q[1];
sx q[1];
rz(-2.5934673) q[1];
sx q[1];
rz(0.60636884) q[1];
rz(2.0229522) q[2];
sx q[2];
rz(-0.24941988) q[2];
sx q[2];
rz(3.1385251) q[2];
rz(-1.8567139) q[3];
sx q[3];
rz(-0.34964041) q[3];
sx q[3];
rz(1.6007363) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
