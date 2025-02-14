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
rz(2.1428406) q[0];
sx q[0];
rz(-2.3698896) q[0];
sx q[0];
rz(-1.9089215) q[0];
rz(-1.9219037) q[1];
sx q[1];
rz(-0.86644679) q[1];
sx q[1];
rz(-2.7629857) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.632376) q[0];
sx q[0];
rz(-1.6086744) q[0];
sx q[0];
rz(1.8011455) q[0];
x q[1];
rz(-0.98388715) q[2];
sx q[2];
rz(-0.68418938) q[2];
sx q[2];
rz(0.049287576) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6191062) q[1];
sx q[1];
rz(-2.2402856) q[1];
sx q[1];
rz(1.8471902) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5903108) q[3];
sx q[3];
rz(-2.4935185) q[3];
sx q[3];
rz(0.66873811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0014701) q[2];
sx q[2];
rz(-1.1727611) q[2];
sx q[2];
rz(-2.6918461) q[2];
rz(0.55284119) q[3];
sx q[3];
rz(-2.5832085) q[3];
sx q[3];
rz(2.9302178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0387892) q[0];
sx q[0];
rz(-1.9961822) q[0];
sx q[0];
rz(0.75622028) q[0];
rz(2.4839632) q[1];
sx q[1];
rz(-2.2537474) q[1];
sx q[1];
rz(2.5530946) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3257826) q[0];
sx q[0];
rz(-1.9414158) q[0];
sx q[0];
rz(-2.4354706) q[0];
rz(-0.21958828) q[2];
sx q[2];
rz(-0.83634085) q[2];
sx q[2];
rz(1.7569923) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5246702) q[1];
sx q[1];
rz(-1.2782885) q[1];
sx q[1];
rz(-3.133417) q[1];
rz(-1.9733834) q[3];
sx q[3];
rz(-0.63791554) q[3];
sx q[3];
rz(1.0653374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6077891) q[2];
sx q[2];
rz(-1.6705931) q[2];
sx q[2];
rz(-0.076347366) q[2];
rz(-2.7401183) q[3];
sx q[3];
rz(-2.6186826) q[3];
sx q[3];
rz(1.128986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56979316) q[0];
sx q[0];
rz(-0.58572584) q[0];
sx q[0];
rz(0.3366003) q[0];
rz(1.0353237) q[1];
sx q[1];
rz(-2.2575049) q[1];
sx q[1];
rz(-0.050315637) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6099358) q[0];
sx q[0];
rz(-2.3667524) q[0];
sx q[0];
rz(-2.4365303) q[0];
rz(-0.29842768) q[2];
sx q[2];
rz(-0.61594916) q[2];
sx q[2];
rz(2.9263934) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.13408537) q[1];
sx q[1];
rz(-1.3614628) q[1];
sx q[1];
rz(0.34761859) q[1];
rz(-pi) q[2];
rz(-0.21304275) q[3];
sx q[3];
rz(-2.1318949) q[3];
sx q[3];
rz(1.1651426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7670373) q[2];
sx q[2];
rz(-0.87696806) q[2];
sx q[2];
rz(2.1602737) q[2];
rz(-2.0902925) q[3];
sx q[3];
rz(-1.4562573) q[3];
sx q[3];
rz(-0.013896996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1023079) q[0];
sx q[0];
rz(-1.6599052) q[0];
sx q[0];
rz(2.1639977) q[0];
rz(-2.5915937) q[1];
sx q[1];
rz(-2.1606052) q[1];
sx q[1];
rz(-2.1994793) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4956168) q[0];
sx q[0];
rz(-1.2496523) q[0];
sx q[0];
rz(-2.7944588) q[0];
rz(1.6245234) q[2];
sx q[2];
rz(-1.2135398) q[2];
sx q[2];
rz(-0.7815278) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7779322) q[1];
sx q[1];
rz(-2.1673551) q[1];
sx q[1];
rz(2.9371757) q[1];
rz(2.3650424) q[3];
sx q[3];
rz(-2.4403493) q[3];
sx q[3];
rz(1.6691735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.00086870988) q[2];
sx q[2];
rz(-2.7694747) q[2];
sx q[2];
rz(-1.5127399) q[2];
rz(0.77423972) q[3];
sx q[3];
rz(-2.1858678) q[3];
sx q[3];
rz(0.16404185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
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
rz(2.0416439) q[0];
sx q[0];
rz(-0.14687563) q[0];
sx q[0];
rz(-2.3351093) q[0];
rz(2.3355314) q[1];
sx q[1];
rz(-1.9652003) q[1];
sx q[1];
rz(0.75256601) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4680639) q[0];
sx q[0];
rz(-1.4554592) q[0];
sx q[0];
rz(-2.2542984) q[0];
rz(0.089293496) q[2];
sx q[2];
rz(-0.39555031) q[2];
sx q[2];
rz(1.7035005) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8991685) q[1];
sx q[1];
rz(-0.43546477) q[1];
sx q[1];
rz(-0.240761) q[1];
rz(1.7349929) q[3];
sx q[3];
rz(-1.1661652) q[3];
sx q[3];
rz(-1.5238289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.26087424) q[2];
sx q[2];
rz(-1.7267092) q[2];
sx q[2];
rz(-0.073337642) q[2];
rz(-0.65555769) q[3];
sx q[3];
rz(-0.95584241) q[3];
sx q[3];
rz(-2.5478794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23621479) q[0];
sx q[0];
rz(-1.0991993) q[0];
sx q[0];
rz(2.4429876) q[0];
rz(1.7504494) q[1];
sx q[1];
rz(-0.51934424) q[1];
sx q[1];
rz(-3.0379675) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2267644) q[0];
sx q[0];
rz(-1.6092759) q[0];
sx q[0];
rz(1.5917042) q[0];
x q[1];
rz(3.1049635) q[2];
sx q[2];
rz(-2.0994791) q[2];
sx q[2];
rz(1.5122459) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.99370804) q[1];
sx q[1];
rz(-1.4698878) q[1];
sx q[1];
rz(1.7624965) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7781939) q[3];
sx q[3];
rz(-1.816664) q[3];
sx q[3];
rz(1.2896001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.708272) q[2];
sx q[2];
rz(-0.89858133) q[2];
sx q[2];
rz(-1.9898604) q[2];
rz(0.084130675) q[3];
sx q[3];
rz(-0.79456544) q[3];
sx q[3];
rz(0.35552037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56650913) q[0];
sx q[0];
rz(-0.36324781) q[0];
sx q[0];
rz(1.890924) q[0];
rz(1.4051215) q[1];
sx q[1];
rz(-2.5650918) q[1];
sx q[1];
rz(1.0692495) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.273928) q[0];
sx q[0];
rz(-1.3153512) q[0];
sx q[0];
rz(-2.5219265) q[0];
rz(3.1282449) q[2];
sx q[2];
rz(-1.0927534) q[2];
sx q[2];
rz(1.0814217) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1508559) q[1];
sx q[1];
rz(-1.6526319) q[1];
sx q[1];
rz(0.86151716) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0642483) q[3];
sx q[3];
rz(-1.3367904) q[3];
sx q[3];
rz(1.4547494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8413267) q[2];
sx q[2];
rz(-1.0043251) q[2];
sx q[2];
rz(-1.1379854) q[2];
rz(2.581572) q[3];
sx q[3];
rz(-1.0617826) q[3];
sx q[3];
rz(1.2573857) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1244125) q[0];
sx q[0];
rz(-0.37473285) q[0];
sx q[0];
rz(2.4770233) q[0];
rz(0.38732227) q[1];
sx q[1];
rz(-2.1699984) q[1];
sx q[1];
rz(-1.9688781) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28139505) q[0];
sx q[0];
rz(-2.144964) q[0];
sx q[0];
rz(-2.6414394) q[0];
rz(-pi) q[1];
rz(1.7936158) q[2];
sx q[2];
rz(-2.0231721) q[2];
sx q[2];
rz(-1.4063175) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.203158) q[1];
sx q[1];
rz(-0.76701346) q[1];
sx q[1];
rz(2.2545283) q[1];
rz(-pi) q[2];
rz(-2.6083192) q[3];
sx q[3];
rz(-1.82996) q[3];
sx q[3];
rz(-2.5854098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.49501255) q[2];
sx q[2];
rz(-2.6436372) q[2];
sx q[2];
rz(-0.23769561) q[2];
rz(-2.2869535) q[3];
sx q[3];
rz(-1.6228638) q[3];
sx q[3];
rz(0.91116992) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89331996) q[0];
sx q[0];
rz(-1.8017636) q[0];
sx q[0];
rz(-2.7082537) q[0];
rz(-1.8386819) q[1];
sx q[1];
rz(-1.7830667) q[1];
sx q[1];
rz(3.0824331) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24640326) q[0];
sx q[0];
rz(-1.8120017) q[0];
sx q[0];
rz(3.1144322) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3174345) q[2];
sx q[2];
rz(-2.4640016) q[2];
sx q[2];
rz(-0.80582011) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8856992) q[1];
sx q[1];
rz(-2.7045193) q[1];
sx q[1];
rz(-2.423362) q[1];
x q[2];
rz(1.9800277) q[3];
sx q[3];
rz(-0.61823119) q[3];
sx q[3];
rz(2.359584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2032418) q[2];
sx q[2];
rz(-1.7277191) q[2];
sx q[2];
rz(2.4436277) q[2];
rz(-0.62567726) q[3];
sx q[3];
rz(-2.8901849) q[3];
sx q[3];
rz(1.3231369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.043592602) q[0];
sx q[0];
rz(-0.55857825) q[0];
sx q[0];
rz(-0.25398764) q[0];
rz(-0.99098539) q[1];
sx q[1];
rz(-1.5373983) q[1];
sx q[1];
rz(-3.1414247) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22826057) q[0];
sx q[0];
rz(-2.0096632) q[0];
sx q[0];
rz(-1.0717546) q[0];
x q[1];
rz(1.9591879) q[2];
sx q[2];
rz(-3.0157308) q[2];
sx q[2];
rz(2.2110155) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.907577) q[1];
sx q[1];
rz(-1.325851) q[1];
sx q[1];
rz(-2.8744389) q[1];
rz(1.1669257) q[3];
sx q[3];
rz(-1.7332826) q[3];
sx q[3];
rz(-1.6222909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.99074546) q[2];
sx q[2];
rz(-1.7898229) q[2];
sx q[2];
rz(3.0899437) q[2];
rz(2.1150186) q[3];
sx q[3];
rz(-2.5995422) q[3];
sx q[3];
rz(2.3718204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6386665) q[0];
sx q[0];
rz(-0.8994871) q[0];
sx q[0];
rz(0.58216397) q[0];
rz(2.7405986) q[1];
sx q[1];
rz(-0.6656701) q[1];
sx q[1];
rz(2.3013339) q[1];
rz(-0.51908334) q[2];
sx q[2];
rz(-1.097403) q[2];
sx q[2];
rz(2.335142) q[2];
rz(2.544315) q[3];
sx q[3];
rz(-0.18571449) q[3];
sx q[3];
rz(-1.0635536) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
