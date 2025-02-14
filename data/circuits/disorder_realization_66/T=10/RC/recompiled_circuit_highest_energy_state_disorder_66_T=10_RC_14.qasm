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
rz(-0.771703) q[0];
sx q[0];
rz(-1.2326711) q[0];
rz(-1.9219037) q[1];
sx q[1];
rz(-0.86644679) q[1];
sx q[1];
rz(-2.7629857) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2216778) q[0];
sx q[0];
rz(-0.23338813) q[0];
sx q[0];
rz(-1.4063157) q[0];
rz(2.7173434) q[2];
sx q[2];
rz(-1.0165809) q[2];
sx q[2];
rz(-0.65997906) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5224865) q[1];
sx q[1];
rz(-0.90130708) q[1];
sx q[1];
rz(1.8471902) q[1];
rz(-pi) q[2];
rz(0.55128184) q[3];
sx q[3];
rz(-2.4935185) q[3];
sx q[3];
rz(-2.4728545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0014701) q[2];
sx q[2];
rz(-1.9688316) q[2];
sx q[2];
rz(0.44974652) q[2];
rz(0.55284119) q[3];
sx q[3];
rz(-0.55838412) q[3];
sx q[3];
rz(0.21137485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0387892) q[0];
sx q[0];
rz(-1.9961822) q[0];
sx q[0];
rz(0.75622028) q[0];
rz(-0.65762949) q[1];
sx q[1];
rz(-0.88784528) q[1];
sx q[1];
rz(0.58849803) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3257826) q[0];
sx q[0];
rz(-1.2001769) q[0];
sx q[0];
rz(-2.4354706) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82423721) q[2];
sx q[2];
rz(-1.7331799) q[2];
sx q[2];
rz(2.8069212) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.6169225) q[1];
sx q[1];
rz(-1.2782885) q[1];
sx q[1];
rz(-0.0081756552) q[1];
rz(-pi) q[2];
rz(-2.858925) q[3];
sx q[3];
rz(-0.99093197) q[3];
sx q[3];
rz(1.5527703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.53380352) q[2];
sx q[2];
rz(-1.4709996) q[2];
sx q[2];
rz(3.0652453) q[2];
rz(0.40147436) q[3];
sx q[3];
rz(-0.52291003) q[3];
sx q[3];
rz(-1.128986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5717995) q[0];
sx q[0];
rz(-2.5558668) q[0];
sx q[0];
rz(0.3366003) q[0];
rz(1.0353237) q[1];
sx q[1];
rz(-0.88408771) q[1];
sx q[1];
rz(0.050315637) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6099358) q[0];
sx q[0];
rz(-2.3667524) q[0];
sx q[0];
rz(-0.70506238) q[0];
x q[1];
rz(2.843165) q[2];
sx q[2];
rz(-2.5256435) q[2];
sx q[2];
rz(0.21519923) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0075073) q[1];
sx q[1];
rz(-1.7801298) q[1];
sx q[1];
rz(-2.7939741) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8953241) q[3];
sx q[3];
rz(-0.59610329) q[3];
sx q[3];
rz(-2.3625797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.1023079) q[0];
sx q[0];
rz(-1.4816875) q[0];
sx q[0];
rz(2.1639977) q[0];
rz(-2.5915937) q[1];
sx q[1];
rz(-2.1606052) q[1];
sx q[1];
rz(0.94211334) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4956168) q[0];
sx q[0];
rz(-1.8919404) q[0];
sx q[0];
rz(-0.34713388) q[0];
rz(-pi) q[1];
rz(-0.35772985) q[2];
sx q[2];
rz(-1.6211281) q[2];
sx q[2];
rz(2.3711287) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7779322) q[1];
sx q[1];
rz(-2.1673551) q[1];
sx q[1];
rz(0.20441695) q[1];
rz(2.5994456) q[3];
sx q[3];
rz(-2.039969) q[3];
sx q[3];
rz(-0.5455324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.00086870988) q[2];
sx q[2];
rz(-2.7694747) q[2];
sx q[2];
rz(1.5127399) q[2];
rz(0.77423972) q[3];
sx q[3];
rz(-0.95572487) q[3];
sx q[3];
rz(2.9775508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.0999488) q[0];
sx q[0];
rz(-2.994717) q[0];
sx q[0];
rz(0.80648333) q[0];
rz(-2.3355314) q[1];
sx q[1];
rz(-1.9652003) q[1];
sx q[1];
rz(2.3890266) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8984932) q[0];
sx q[0];
rz(-0.69161915) q[0];
sx q[0];
rz(1.7522295) q[0];
rz(0.3941329) q[2];
sx q[2];
rz(-1.536429) q[2];
sx q[2];
rz(0.050273035) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.10933398) q[1];
sx q[1];
rz(-1.671549) q[1];
sx q[1];
rz(2.7172158) q[1];
x q[2];
rz(-2.7769776) q[3];
sx q[3];
rz(-0.43496736) q[3];
sx q[3];
rz(-1.2193958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.26087424) q[2];
sx q[2];
rz(-1.7267092) q[2];
sx q[2];
rz(3.068255) q[2];
rz(0.65555769) q[3];
sx q[3];
rz(-0.95584241) q[3];
sx q[3];
rz(2.5478794) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23621479) q[0];
sx q[0];
rz(-1.0991993) q[0];
sx q[0];
rz(0.69860506) q[0];
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
rz(0.34322748) q[0];
sx q[0];
rz(-1.5916887) q[0];
sx q[0];
rz(-3.1031046) q[0];
rz(-pi) q[1];
rz(-0.03662911) q[2];
sx q[2];
rz(-2.0994791) q[2];
sx q[2];
rz(-1.6293467) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.99370804) q[1];
sx q[1];
rz(-1.6717049) q[1];
sx q[1];
rz(-1.3790961) q[1];
rz(2.4544434) q[3];
sx q[3];
rz(-0.32029974) q[3];
sx q[3];
rz(-2.564858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.43332064) q[2];
sx q[2];
rz(-0.89858133) q[2];
sx q[2];
rz(1.9898604) q[2];
rz(0.084130675) q[3];
sx q[3];
rz(-2.3470272) q[3];
sx q[3];
rz(2.7860723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5750835) q[0];
sx q[0];
rz(-2.7783448) q[0];
sx q[0];
rz(-1.2506686) q[0];
rz(-1.7364712) q[1];
sx q[1];
rz(-0.57650081) q[1];
sx q[1];
rz(-1.0692495) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0980643) q[0];
sx q[0];
rz(-2.4777924) q[0];
sx q[0];
rz(0.4225771) q[0];
x q[1];
rz(1.5450412) q[2];
sx q[2];
rz(-0.47821486) q[2];
sx q[2];
rz(-1.1104294) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7916966) q[1];
sx q[1];
rz(-0.86438599) q[1];
sx q[1];
rz(0.10766761) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.077344374) q[3];
sx q[3];
rz(-1.3367904) q[3];
sx q[3];
rz(-1.4547494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8413267) q[2];
sx q[2];
rz(-1.0043251) q[2];
sx q[2];
rz(2.0036073) q[2];
rz(0.56002069) q[3];
sx q[3];
rz(-2.0798101) q[3];
sx q[3];
rz(-1.8842069) q[3];
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
x q[2];
rz(-pi/2) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5636503) q[0];
sx q[0];
rz(-1.9852173) q[0];
sx q[0];
rz(-2.2060288) q[0];
rz(2.6792999) q[2];
sx q[2];
rz(-1.7708808) q[2];
sx q[2];
rz(-2.878396) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.091306578) q[1];
sx q[1];
rz(-1.0027409) q[1];
sx q[1];
rz(0.54692366) q[1];
rz(1.8694626) q[3];
sx q[3];
rz(-2.0844679) q[3];
sx q[3];
rz(0.86446867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6465801) q[2];
sx q[2];
rz(-0.4979555) q[2];
sx q[2];
rz(2.903897) q[2];
rz(2.2869535) q[3];
sx q[3];
rz(-1.5187289) q[3];
sx q[3];
rz(-2.2304227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89331996) q[0];
sx q[0];
rz(-1.8017636) q[0];
sx q[0];
rz(2.7082537) q[0];
rz(1.3029107) q[1];
sx q[1];
rz(-1.7830667) q[1];
sx q[1];
rz(-0.059159577) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35964701) q[0];
sx q[0];
rz(-0.24270013) q[0];
sx q[0];
rz(1.6807483) q[0];
rz(0.19903175) q[2];
sx q[2];
rz(-2.2229554) q[2];
sx q[2];
rz(-2.6566344) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2558935) q[1];
sx q[1];
rz(-2.7045193) q[1];
sx q[1];
rz(0.71823068) q[1];
rz(-2.8657984) q[3];
sx q[3];
rz(-1.0101476) q[3];
sx q[3];
rz(-1.8704853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9383508) q[2];
sx q[2];
rz(-1.7277191) q[2];
sx q[2];
rz(-2.4436277) q[2];
rz(0.62567726) q[3];
sx q[3];
rz(-0.2514078) q[3];
sx q[3];
rz(-1.8184557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0980001) q[0];
sx q[0];
rz(-0.55857825) q[0];
sx q[0];
rz(-0.25398764) q[0];
rz(-0.99098539) q[1];
sx q[1];
rz(-1.5373983) q[1];
sx q[1];
rz(-3.1414247) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9133321) q[0];
sx q[0];
rz(-1.1319295) q[0];
sx q[0];
rz(1.0717546) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0937121) q[2];
sx q[2];
rz(-1.6872395) q[2];
sx q[2];
rz(-2.6021985) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7533299) q[1];
sx q[1];
rz(-2.7811353) q[1];
sx q[1];
rz(0.75812104) q[1];
x q[2];
rz(1.9746669) q[3];
sx q[3];
rz(-1.40831) q[3];
sx q[3];
rz(1.5193017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1508472) q[2];
sx q[2];
rz(-1.7898229) q[2];
sx q[2];
rz(-0.051648971) q[2];
rz(2.1150186) q[3];
sx q[3];
rz(-0.54205042) q[3];
sx q[3];
rz(-2.3718204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6386665) q[0];
sx q[0];
rz(-2.2421056) q[0];
sx q[0];
rz(-2.5594287) q[0];
rz(0.40099405) q[1];
sx q[1];
rz(-2.4759226) q[1];
sx q[1];
rz(-0.84025875) q[1];
rz(2.1037965) q[2];
sx q[2];
rz(-1.1134951) q[2];
sx q[2];
rz(-2.122428) q[2];
rz(-1.6760679) q[3];
sx q[3];
rz(-1.7240763) q[3];
sx q[3];
rz(-0.45818466) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
