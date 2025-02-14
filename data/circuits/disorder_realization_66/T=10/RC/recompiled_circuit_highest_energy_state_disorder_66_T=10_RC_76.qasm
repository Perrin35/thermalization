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
rz(-0.99875206) q[0];
sx q[0];
rz(2.3698896) q[0];
sx q[0];
rz(10.657449) q[0];
rz(-1.9219037) q[1];
sx q[1];
rz(-0.86644679) q[1];
sx q[1];
rz(0.37860695) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0526992) q[0];
sx q[0];
rz(-1.3406154) q[0];
sx q[0];
rz(-3.102688) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1673752) q[2];
sx q[2];
rz(-1.2132036) q[2];
sx q[2];
rz(1.9973988) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.12593658) q[1];
sx q[1];
rz(-1.7864461) q[1];
sx q[1];
rz(-2.4531948) q[1];
rz(-pi) q[2];
rz(-0.572851) q[3];
sx q[3];
rz(-1.2490954) q[3];
sx q[3];
rz(2.6953146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0014701) q[2];
sx q[2];
rz(-1.1727611) q[2];
sx q[2];
rz(0.44974652) q[2];
rz(-2.5887515) q[3];
sx q[3];
rz(-2.5832085) q[3];
sx q[3];
rz(-0.21137485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1028035) q[0];
sx q[0];
rz(-1.9961822) q[0];
sx q[0];
rz(0.75622028) q[0];
rz(-2.4839632) q[1];
sx q[1];
rz(-2.2537474) q[1];
sx q[1];
rz(0.58849803) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.156608) q[0];
sx q[0];
rz(-2.3591908) q[0];
sx q[0];
rz(-0.53956145) q[0];
rz(-pi) q[1];
rz(-0.21958828) q[2];
sx q[2];
rz(-0.83634085) q[2];
sx q[2];
rz(-1.3846004) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.6169225) q[1];
sx q[1];
rz(-1.8633042) q[1];
sx q[1];
rz(0.0081756552) q[1];
x q[2];
rz(-2.858925) q[3];
sx q[3];
rz(-2.1506607) q[3];
sx q[3];
rz(-1.5527703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.53380352) q[2];
sx q[2];
rz(-1.4709996) q[2];
sx q[2];
rz(0.076347366) q[2];
rz(-0.40147436) q[3];
sx q[3];
rz(-0.52291003) q[3];
sx q[3];
rz(-2.0126066) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5717995) q[0];
sx q[0];
rz(-0.58572584) q[0];
sx q[0];
rz(-2.8049923) q[0];
rz(-1.0353237) q[1];
sx q[1];
rz(-2.2575049) q[1];
sx q[1];
rz(0.050315637) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49281989) q[0];
sx q[0];
rz(-2.0413715) q[0];
sx q[0];
rz(2.5008766) q[0];
rz(2.5467954) q[2];
sx q[2];
rz(-1.7414879) q[2];
sx q[2];
rz(1.5399982) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.13408537) q[1];
sx q[1];
rz(-1.3614628) q[1];
sx q[1];
rz(-2.7939741) q[1];
rz(-pi) q[2];
x q[2];
rz(0.21304275) q[3];
sx q[3];
rz(-1.0096978) q[3];
sx q[3];
rz(1.1651426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3745554) q[2];
sx q[2];
rz(-2.2646246) q[2];
sx q[2];
rz(0.98131895) q[2];
rz(-2.0902925) q[3];
sx q[3];
rz(-1.6853354) q[3];
sx q[3];
rz(0.013896996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0392847) q[0];
sx q[0];
rz(-1.6599052) q[0];
sx q[0];
rz(-0.977595) q[0];
rz(-2.5915937) q[1];
sx q[1];
rz(-2.1606052) q[1];
sx q[1];
rz(-2.1994793) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4993663) q[0];
sx q[0];
rz(-2.6731579) q[0];
sx q[0];
rz(-0.77418615) q[0];
rz(-pi) q[1];
rz(-2.7838628) q[2];
sx q[2];
rz(-1.5204645) q[2];
sx q[2];
rz(2.3711287) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0504024) q[1];
sx q[1];
rz(-1.7395294) q[1];
sx q[1];
rz(0.96442142) q[1];
rz(-pi) q[2];
rz(-2.5994456) q[3];
sx q[3];
rz(-2.039969) q[3];
sx q[3];
rz(0.5455324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1407239) q[2];
sx q[2];
rz(-0.372118) q[2];
sx q[2];
rz(-1.6288527) q[2];
rz(2.3673529) q[3];
sx q[3];
rz(-0.95572487) q[3];
sx q[3];
rz(-2.9775508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0999488) q[0];
sx q[0];
rz(-2.994717) q[0];
sx q[0];
rz(0.80648333) q[0];
rz(0.80606127) q[1];
sx q[1];
rz(-1.9652003) q[1];
sx q[1];
rz(-0.75256601) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8984932) q[0];
sx q[0];
rz(-0.69161915) q[0];
sx q[0];
rz(-1.3893632) q[0];
x q[1];
rz(-0.089293496) q[2];
sx q[2];
rz(-0.39555031) q[2];
sx q[2];
rz(1.4380921) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0322587) q[1];
sx q[1];
rz(-1.671549) q[1];
sx q[1];
rz(-0.42437683) q[1];
rz(-pi) q[2];
rz(2.7320382) q[3];
sx q[3];
rz(-1.7216276) q[3];
sx q[3];
rz(0.018168381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.26087424) q[2];
sx q[2];
rz(-1.7267092) q[2];
sx q[2];
rz(-0.073337642) q[2];
rz(-2.486035) q[3];
sx q[3];
rz(-2.1857502) q[3];
sx q[3];
rz(-2.5478794) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23621479) q[0];
sx q[0];
rz(-2.0423934) q[0];
sx q[0];
rz(-0.69860506) q[0];
rz(1.7504494) q[1];
sx q[1];
rz(-0.51934424) q[1];
sx q[1];
rz(0.10362518) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2267644) q[0];
sx q[0];
rz(-1.6092759) q[0];
sx q[0];
rz(1.5917042) q[0];
rz(-1.5081872) q[2];
sx q[2];
rz(-0.52982989) q[2];
sx q[2];
rz(1.7018715) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1478846) q[1];
sx q[1];
rz(-1.4698878) q[1];
sx q[1];
rz(-1.3790961) q[1];
rz(1.7781939) q[3];
sx q[3];
rz(-1.3249287) q[3];
sx q[3];
rz(1.2896001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.43332064) q[2];
sx q[2];
rz(-2.2430113) q[2];
sx q[2];
rz(1.9898604) q[2];
rz(-3.057462) q[3];
sx q[3];
rz(-0.79456544) q[3];
sx q[3];
rz(-2.7860723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5750835) q[0];
sx q[0];
rz(-2.7783448) q[0];
sx q[0];
rz(-1.890924) q[0];
rz(-1.7364712) q[1];
sx q[1];
rz(-0.57650081) q[1];
sx q[1];
rz(-1.0692495) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.273928) q[0];
sx q[0];
rz(-1.3153512) q[0];
sx q[0];
rz(-2.5219265) q[0];
x q[1];
rz(0.013347722) q[2];
sx q[2];
rz(-1.0927534) q[2];
sx q[2];
rz(-1.0814217) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9907367) q[1];
sx q[1];
rz(-1.6526319) q[1];
sx q[1];
rz(2.2800755) q[1];
x q[2];
rz(-1.8054784) q[3];
sx q[3];
rz(-1.6460287) q[3];
sx q[3];
rz(-3.0075775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.30026597) q[2];
sx q[2];
rz(-2.1372676) q[2];
sx q[2];
rz(-2.0036073) q[2];
rz(0.56002069) q[3];
sx q[3];
rz(-2.0798101) q[3];
sx q[3];
rz(-1.8842069) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1244125) q[0];
sx q[0];
rz(-2.7668598) q[0];
sx q[0];
rz(0.66456932) q[0];
rz(-0.38732227) q[1];
sx q[1];
rz(-2.1699984) q[1];
sx q[1];
rz(1.9688781) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8601976) q[0];
sx q[0];
rz(-2.144964) q[0];
sx q[0];
rz(0.50015323) q[0];
x q[1];
rz(0.42675777) q[2];
sx q[2];
rz(-0.50084693) q[2];
sx q[2];
rz(-0.92808135) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9786611) q[1];
sx q[1];
rz(-2.0245981) q[1];
sx q[1];
rz(-2.2125208) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6083192) q[3];
sx q[3];
rz(-1.3116326) q[3];
sx q[3];
rz(-0.55618286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6465801) q[2];
sx q[2];
rz(-2.6436372) q[2];
sx q[2];
rz(-0.23769561) q[2];
rz(-0.85463917) q[3];
sx q[3];
rz(-1.6228638) q[3];
sx q[3];
rz(2.2304227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.2482727) q[0];
sx q[0];
rz(-1.3398291) q[0];
sx q[0];
rz(-2.7082537) q[0];
rz(-1.8386819) q[1];
sx q[1];
rz(-1.358526) q[1];
sx q[1];
rz(0.059159577) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.823689) q[0];
sx q[0];
rz(-1.5971703) q[0];
sx q[0];
rz(-1.8120873) q[0];
x q[1];
rz(-2.9425609) q[2];
sx q[2];
rz(-2.2229554) q[2];
sx q[2];
rz(0.48495822) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6530919) q[1];
sx q[1];
rz(-1.8951804) q[1];
sx q[1];
rz(1.2725194) q[1];
rz(2.1489359) q[3];
sx q[3];
rz(-1.8035144) q[3];
sx q[3];
rz(-0.44905277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9383508) q[2];
sx q[2];
rz(-1.4138736) q[2];
sx q[2];
rz(-0.69796491) q[2];
rz(0.62567726) q[3];
sx q[3];
rz(-2.8901849) q[3];
sx q[3];
rz(-1.3231369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.043592602) q[0];
sx q[0];
rz(-0.55857825) q[0];
sx q[0];
rz(0.25398764) q[0];
rz(0.99098539) q[1];
sx q[1];
rz(-1.6041944) q[1];
sx q[1];
rz(-3.1414247) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5714665) q[0];
sx q[0];
rz(-1.1227221) q[0];
sx q[0];
rz(-0.49093935) q[0];
rz(1.1824047) q[2];
sx q[2];
rz(-3.0157308) q[2];
sx q[2];
rz(-2.2110155) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.40305576) q[1];
sx q[1];
rz(-1.8297894) q[1];
sx q[1];
rz(1.8243755) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17642085) q[3];
sx q[3];
rz(-1.9690445) q[3];
sx q[3];
rz(-0.12052025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1508472) q[2];
sx q[2];
rz(-1.3517697) q[2];
sx q[2];
rz(-3.0899437) q[2];
rz(-1.026574) q[3];
sx q[3];
rz(-2.5995422) q[3];
sx q[3];
rz(-0.76977229) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50292618) q[0];
sx q[0];
rz(-0.8994871) q[0];
sx q[0];
rz(0.58216397) q[0];
rz(0.40099405) q[1];
sx q[1];
rz(-2.4759226) q[1];
sx q[1];
rz(-0.84025875) q[1];
rz(2.3401716) q[2];
sx q[2];
rz(-0.68757551) q[2];
sx q[2];
rz(0.09085169) q[2];
rz(-0.15411986) q[3];
sx q[3];
rz(-1.4667635) q[3];
sx q[3];
rz(1.0964805) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
