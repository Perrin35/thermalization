OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.90280688) q[0];
sx q[0];
rz(-1.1385318) q[0];
sx q[0];
rz(-2.0432552) q[0];
rz(1.1605473) q[1];
sx q[1];
rz(2.0055973) q[1];
sx q[1];
rz(10.027167) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1576958) q[0];
sx q[0];
rz(-2.0838782) q[0];
sx q[0];
rz(2.0747572) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6260398) q[2];
sx q[2];
rz(-0.85795244) q[2];
sx q[2];
rz(1.6875658) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8849395) q[1];
sx q[1];
rz(-2.6511484) q[1];
sx q[1];
rz(-1.8153193) q[1];
rz(-pi) q[2];
rz(3.1040339) q[3];
sx q[3];
rz(-1.9989487) q[3];
sx q[3];
rz(-1.1809837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7263553) q[2];
sx q[2];
rz(-1.0865612) q[2];
sx q[2];
rz(1.6437257) q[2];
rz(-3.030153) q[3];
sx q[3];
rz(-1.2846416) q[3];
sx q[3];
rz(2.1845412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9621256) q[0];
sx q[0];
rz(-2.2158556) q[0];
sx q[0];
rz(-2.977648) q[0];
rz(-0.97490772) q[1];
sx q[1];
rz(-0.68760741) q[1];
sx q[1];
rz(1.499739) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1346171) q[0];
sx q[0];
rz(-2.3762759) q[0];
sx q[0];
rz(1.2045443) q[0];
rz(-pi) q[1];
rz(-0.87514211) q[2];
sx q[2];
rz(-0.89584699) q[2];
sx q[2];
rz(0.067719134) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.30399366) q[1];
sx q[1];
rz(-1.6483867) q[1];
sx q[1];
rz(-2.4291888) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1466931) q[3];
sx q[3];
rz(-0.38299503) q[3];
sx q[3];
rz(0.3037616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6770596) q[2];
sx q[2];
rz(-0.50731069) q[2];
sx q[2];
rz(0.23059174) q[2];
rz(0.683189) q[3];
sx q[3];
rz(-2.4468827) q[3];
sx q[3];
rz(2.3901239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7815919) q[0];
sx q[0];
rz(-1.8311904) q[0];
sx q[0];
rz(0.42136425) q[0];
rz(1.6257809) q[1];
sx q[1];
rz(-2.6457364) q[1];
sx q[1];
rz(-2.1741507) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66146278) q[0];
sx q[0];
rz(-2.7217743) q[0];
sx q[0];
rz(-0.90645193) q[0];
rz(-pi) q[1];
rz(0.17379667) q[2];
sx q[2];
rz(-1.1090384) q[2];
sx q[2];
rz(0.51848903) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.9989717) q[1];
sx q[1];
rz(-1.2983067) q[1];
sx q[1];
rz(0.15010826) q[1];
rz(-pi) q[2];
rz(1.2382965) q[3];
sx q[3];
rz(-1.0345113) q[3];
sx q[3];
rz(-0.09418776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7385873) q[2];
sx q[2];
rz(-1.0120729) q[2];
sx q[2];
rz(-2.9834413) q[2];
rz(2.3824597) q[3];
sx q[3];
rz(-1.4594376) q[3];
sx q[3];
rz(-0.37091836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89449969) q[0];
sx q[0];
rz(-0.58421725) q[0];
sx q[0];
rz(-1.9547113) q[0];
rz(3.0184556) q[1];
sx q[1];
rz(-1.5712527) q[1];
sx q[1];
rz(-1.311696) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7191221) q[0];
sx q[0];
rz(-1.4234237) q[0];
sx q[0];
rz(3.0585994) q[0];
rz(1.4940631) q[2];
sx q[2];
rz(-1.4686246) q[2];
sx q[2];
rz(0.42021423) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32320106) q[1];
sx q[1];
rz(-1.9889327) q[1];
sx q[1];
rz(2.0050603) q[1];
rz(-pi) q[2];
rz(1.4315286) q[3];
sx q[3];
rz(-0.0041120681) q[3];
sx q[3];
rz(-2.508956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9436403) q[2];
sx q[2];
rz(-1.8215048) q[2];
sx q[2];
rz(0.8521592) q[2];
rz(-0.85396829) q[3];
sx q[3];
rz(-1.9210509) q[3];
sx q[3];
rz(1.0796237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78344408) q[0];
sx q[0];
rz(-1.3842979) q[0];
sx q[0];
rz(-3.070991) q[0];
rz(-2.103503) q[1];
sx q[1];
rz(-1.9990653) q[1];
sx q[1];
rz(-0.80931726) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22915086) q[0];
sx q[0];
rz(-0.99780592) q[0];
sx q[0];
rz(-1.6567162) q[0];
x q[1];
rz(0.74639928) q[2];
sx q[2];
rz(-2.0944907) q[2];
sx q[2];
rz(-0.46800287) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.63974158) q[1];
sx q[1];
rz(-1.0305522) q[1];
sx q[1];
rz(-1.2667202) q[1];
x q[2];
rz(-2.9633459) q[3];
sx q[3];
rz(-2.2615454) q[3];
sx q[3];
rz(2.8028298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.91850963) q[2];
sx q[2];
rz(-2.5727936) q[2];
sx q[2];
rz(1.9522379) q[2];
rz(0.14063028) q[3];
sx q[3];
rz(-0.46969241) q[3];
sx q[3];
rz(-0.3895337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73914948) q[0];
sx q[0];
rz(-0.91366714) q[0];
sx q[0];
rz(-1.6074578) q[0];
rz(2.2588579) q[1];
sx q[1];
rz(-2.2969756) q[1];
sx q[1];
rz(2.5260063) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.405571) q[0];
sx q[0];
rz(-1.478749) q[0];
sx q[0];
rz(2.7476601) q[0];
x q[1];
rz(0.65290218) q[2];
sx q[2];
rz(-0.98207475) q[2];
sx q[2];
rz(-1.9984286) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2193377) q[1];
sx q[1];
rz(-0.23592792) q[1];
sx q[1];
rz(-2.1990945) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1911931) q[3];
sx q[3];
rz(-1.5801893) q[3];
sx q[3];
rz(-0.91419807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.1659871) q[2];
sx q[2];
rz(-2.1773982) q[2];
sx q[2];
rz(-1.1791505) q[2];
rz(-2.9278582) q[3];
sx q[3];
rz(-1.7515747) q[3];
sx q[3];
rz(-0.24553044) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8971276) q[0];
sx q[0];
rz(-1.9944958) q[0];
sx q[0];
rz(3.097528) q[0];
rz(-2.7524718) q[1];
sx q[1];
rz(-2.6521284) q[1];
sx q[1];
rz(-2.3997831) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76352966) q[0];
sx q[0];
rz(-0.46398315) q[0];
sx q[0];
rz(0.65391175) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9474232) q[2];
sx q[2];
rz(-2.3583226) q[2];
sx q[2];
rz(-0.3585836) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.85838137) q[1];
sx q[1];
rz(-1.2492531) q[1];
sx q[1];
rz(1.47846) q[1];
x q[2];
rz(2.9030475) q[3];
sx q[3];
rz(-2.425634) q[3];
sx q[3];
rz(-1.3357287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6645708) q[2];
sx q[2];
rz(-1.7844113) q[2];
sx q[2];
rz(-0.50194293) q[2];
rz(-1.7160412) q[3];
sx q[3];
rz(-2.612096) q[3];
sx q[3];
rz(0.79290032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3821756) q[0];
sx q[0];
rz(-1.4074396) q[0];
sx q[0];
rz(-2.7446246) q[0];
rz(1.601903) q[1];
sx q[1];
rz(-2.868728) q[1];
sx q[1];
rz(2.4698965) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7714193) q[0];
sx q[0];
rz(-1.1259997) q[0];
sx q[0];
rz(-3.0014787) q[0];
rz(-2.0136859) q[2];
sx q[2];
rz(-1.9077565) q[2];
sx q[2];
rz(0.86293399) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7377598) q[1];
sx q[1];
rz(-1.2602401) q[1];
sx q[1];
rz(2.4871189) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.83064744) q[3];
sx q[3];
rz(-1.1201356) q[3];
sx q[3];
rz(-0.78307952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4063065) q[2];
sx q[2];
rz(-0.51484171) q[2];
sx q[2];
rz(-1.7740645) q[2];
rz(-0.96274084) q[3];
sx q[3];
rz(-1.5509501) q[3];
sx q[3];
rz(-0.65403691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.3108567) q[0];
sx q[0];
rz(-1.5565358) q[0];
sx q[0];
rz(-2.8209525) q[0];
rz(2.2036208) q[1];
sx q[1];
rz(-1.9357977) q[1];
sx q[1];
rz(-1.4042312) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36286255) q[0];
sx q[0];
rz(-1.6741279) q[0];
sx q[0];
rz(-0.10424239) q[0];
rz(-pi) q[1];
rz(-1.5455568) q[2];
sx q[2];
rz(-1.2938616) q[2];
sx q[2];
rz(-1.0482839) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.117489) q[1];
sx q[1];
rz(-1.5507332) q[1];
sx q[1];
rz(-1.5320918) q[1];
x q[2];
rz(-2.8912084) q[3];
sx q[3];
rz(-2.3634849) q[3];
sx q[3];
rz(2.476256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2001026) q[2];
sx q[2];
rz(-0.76277554) q[2];
sx q[2];
rz(0.74083677) q[2];
rz(-1.0624933) q[3];
sx q[3];
rz(-1.1081568) q[3];
sx q[3];
rz(1.197193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9366539) q[0];
sx q[0];
rz(-2.7118201) q[0];
sx q[0];
rz(0.75576654) q[0];
rz(1.8273805) q[1];
sx q[1];
rz(-1.5145489) q[1];
sx q[1];
rz(-0.72602138) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7814815) q[0];
sx q[0];
rz(-1.5364293) q[0];
sx q[0];
rz(-2.3590259) q[0];
rz(-pi) q[1];
rz(1.6976895) q[2];
sx q[2];
rz(-0.68918258) q[2];
sx q[2];
rz(0.5529595) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0997114) q[1];
sx q[1];
rz(-1.2234294) q[1];
sx q[1];
rz(0.16176407) q[1];
rz(1.8733379) q[3];
sx q[3];
rz(-2.6740251) q[3];
sx q[3];
rz(-2.9085161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8119592) q[2];
sx q[2];
rz(-2.6128431) q[2];
sx q[2];
rz(0.94919666) q[2];
rz(2.9694929) q[3];
sx q[3];
rz(-0.97730079) q[3];
sx q[3];
rz(2.7026091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028845499) q[0];
sx q[0];
rz(-1.5929359) q[0];
sx q[0];
rz(1.5449217) q[0];
rz(1.159809) q[1];
sx q[1];
rz(-1.8198967) q[1];
sx q[1];
rz(1.5130704) q[1];
rz(0.25111689) q[2];
sx q[2];
rz(-2.579183) q[2];
sx q[2];
rz(-0.05280799) q[2];
rz(-1.0817901) q[3];
sx q[3];
rz(-1.5086969) q[3];
sx q[3];
rz(2.8668565) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
