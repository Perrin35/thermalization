OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.43006858) q[0];
sx q[0];
rz(-3.0741337) q[0];
sx q[0];
rz(-0.67396069) q[0];
rz(2.824488) q[1];
sx q[1];
rz(-1.5082521) q[1];
sx q[1];
rz(-2.34692) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90097839) q[0];
sx q[0];
rz(-0.52871791) q[0];
sx q[0];
rz(2.3762977) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27172471) q[2];
sx q[2];
rz(-0.50479111) q[2];
sx q[2];
rz(1.7287849) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.049584576) q[1];
sx q[1];
rz(-1.2283748) q[1];
sx q[1];
rz(-0.53978668) q[1];
x q[2];
rz(-1.2875597) q[3];
sx q[3];
rz(-2.6784416) q[3];
sx q[3];
rz(1.2676261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.47544605) q[2];
sx q[2];
rz(-1.1300766) q[2];
sx q[2];
rz(1.8165992) q[2];
rz(0.31630668) q[3];
sx q[3];
rz(-0.91747147) q[3];
sx q[3];
rz(-0.90707183) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9544202) q[0];
sx q[0];
rz(-2.3311054) q[0];
sx q[0];
rz(-2.1818838) q[0];
rz(-1.6540487) q[1];
sx q[1];
rz(-2.4580749) q[1];
sx q[1];
rz(0.62746343) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42022959) q[0];
sx q[0];
rz(-0.010275928) q[0];
sx q[0];
rz(-2.4179439) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9299514) q[2];
sx q[2];
rz(-1.5492808) q[2];
sx q[2];
rz(-0.053552901) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5261425) q[1];
sx q[1];
rz(-2.8612988) q[1];
sx q[1];
rz(1.0326833) q[1];
x q[2];
rz(-0.012410951) q[3];
sx q[3];
rz(-1.59613) q[3];
sx q[3];
rz(-0.94143553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.54005694) q[2];
sx q[2];
rz(-1.7404218) q[2];
sx q[2];
rz(-0.97989782) q[2];
rz(-1.4510138) q[3];
sx q[3];
rz(-2.4308949) q[3];
sx q[3];
rz(-1.4427982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56101218) q[0];
sx q[0];
rz(-3.09364) q[0];
sx q[0];
rz(1.3797492) q[0];
rz(2.9648932) q[1];
sx q[1];
rz(-1.9243762) q[1];
sx q[1];
rz(1.3476936) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5603148) q[0];
sx q[0];
rz(-0.53239765) q[0];
sx q[0];
rz(0.48849948) q[0];
rz(-pi) q[1];
rz(-0.10920306) q[2];
sx q[2];
rz(-0.20143992) q[2];
sx q[2];
rz(2.2815454) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.99464097) q[1];
sx q[1];
rz(-0.1576345) q[1];
sx q[1];
rz(0.68871246) q[1];
x q[2];
rz(-1.755581) q[3];
sx q[3];
rz(-1.4818958) q[3];
sx q[3];
rz(0.88528663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1220864) q[2];
sx q[2];
rz(-2.2821189) q[2];
sx q[2];
rz(-1.9879831) q[2];
rz(0.38976088) q[3];
sx q[3];
rz(-0.4699769) q[3];
sx q[3];
rz(-0.78891689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47238123) q[0];
sx q[0];
rz(-0.62494576) q[0];
sx q[0];
rz(2.5770082) q[0];
rz(3.0584884) q[1];
sx q[1];
rz(-1.233498) q[1];
sx q[1];
rz(-2.9262503) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4190061) q[0];
sx q[0];
rz(-2.3424087) q[0];
sx q[0];
rz(-1.0976085) q[0];
rz(-pi) q[1];
rz(-2.1583546) q[2];
sx q[2];
rz(-1.431627) q[2];
sx q[2];
rz(0.61901865) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.54289651) q[1];
sx q[1];
rz(-2.9196432) q[1];
sx q[1];
rz(1.0735682) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68656355) q[3];
sx q[3];
rz(-2.4719704) q[3];
sx q[3];
rz(1.3470105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.6354436) q[2];
sx q[2];
rz(-2.052867) q[2];
sx q[2];
rz(2.5396458) q[2];
rz(0.69058949) q[3];
sx q[3];
rz(-0.75406051) q[3];
sx q[3];
rz(0.62600342) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2895806) q[0];
sx q[0];
rz(-0.99442011) q[0];
sx q[0];
rz(-2.1767298) q[0];
rz(-0.016618641) q[1];
sx q[1];
rz(-2.7162572) q[1];
sx q[1];
rz(-0.66666493) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4283838) q[0];
sx q[0];
rz(-0.41233006) q[0];
sx q[0];
rz(-1.8250188) q[0];
rz(-pi) q[1];
rz(0.35258106) q[2];
sx q[2];
rz(-1.0513628) q[2];
sx q[2];
rz(2.7378766) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1946698) q[1];
sx q[1];
rz(-1.2877081) q[1];
sx q[1];
rz(-2.4592295) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8673709) q[3];
sx q[3];
rz(-1.7843102) q[3];
sx q[3];
rz(-0.34542686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.91606402) q[2];
sx q[2];
rz(-1.3801882) q[2];
sx q[2];
rz(1.2529681) q[2];
rz(1.3145674) q[3];
sx q[3];
rz(-1.1809177) q[3];
sx q[3];
rz(-1.130828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0261633) q[0];
sx q[0];
rz(-0.9342497) q[0];
sx q[0];
rz(1.3622267) q[0];
rz(1.7419787) q[1];
sx q[1];
rz(-1.5796966) q[1];
sx q[1];
rz(1.6113575) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49028542) q[0];
sx q[0];
rz(-1.4978652) q[0];
sx q[0];
rz(-0.42071995) q[0];
rz(-pi) q[1];
rz(1.8958695) q[2];
sx q[2];
rz(-1.1522066) q[2];
sx q[2];
rz(-2.3649154) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.64165243) q[1];
sx q[1];
rz(-1.3556726) q[1];
sx q[1];
rz(-1.1670477) q[1];
x q[2];
rz(-1.3774032) q[3];
sx q[3];
rz(-2.8838727) q[3];
sx q[3];
rz(1.6267488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7851012) q[2];
sx q[2];
rz(-2.1337528) q[2];
sx q[2];
rz(0.029416857) q[2];
rz(-1.6507089) q[3];
sx q[3];
rz(-2.1865032) q[3];
sx q[3];
rz(-1.9194549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4239663) q[0];
sx q[0];
rz(-0.76537776) q[0];
sx q[0];
rz(0.18280612) q[0];
rz(-2.706066) q[1];
sx q[1];
rz(-0.85414779) q[1];
sx q[1];
rz(-2.6307154) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8683838) q[0];
sx q[0];
rz(-1.8330935) q[0];
sx q[0];
rz(-1.457085) q[0];
x q[1];
rz(0.80656959) q[2];
sx q[2];
rz(-0.91420805) q[2];
sx q[2];
rz(-0.29733959) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4527013) q[1];
sx q[1];
rz(-0.93033067) q[1];
sx q[1];
rz(-2.074261) q[1];
rz(-pi) q[2];
rz(-1.0484344) q[3];
sx q[3];
rz(-1.2921385) q[3];
sx q[3];
rz(0.1230965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9110979) q[2];
sx q[2];
rz(-2.348933) q[2];
sx q[2];
rz(2.5210209) q[2];
rz(0.42256045) q[3];
sx q[3];
rz(-1.3603323) q[3];
sx q[3];
rz(1.2703936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7997416) q[0];
sx q[0];
rz(-0.36062476) q[0];
sx q[0];
rz(-1.0890695) q[0];
rz(-1.0808806) q[1];
sx q[1];
rz(-1.8573449) q[1];
sx q[1];
rz(2.506315) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3978183) q[0];
sx q[0];
rz(-1.9282856) q[0];
sx q[0];
rz(0.1648358) q[0];
x q[1];
rz(-2.9674171) q[2];
sx q[2];
rz(-1.0798228) q[2];
sx q[2];
rz(-0.27506405) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2891846) q[1];
sx q[1];
rz(-1.8361366) q[1];
sx q[1];
rz(0.62383382) q[1];
rz(-pi) q[2];
rz(1.876272) q[3];
sx q[3];
rz(-2.5442122) q[3];
sx q[3];
rz(1.3437831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13715956) q[2];
sx q[2];
rz(-1.6798423) q[2];
sx q[2];
rz(0.00017246406) q[2];
rz(0.6862644) q[3];
sx q[3];
rz(-0.72824794) q[3];
sx q[3];
rz(0.85787684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2962608) q[0];
sx q[0];
rz(-2.2289472) q[0];
sx q[0];
rz(-2.4849179) q[0];
rz(0.36390057) q[1];
sx q[1];
rz(-2.6952422) q[1];
sx q[1];
rz(-2.231853) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37659392) q[0];
sx q[0];
rz(-1.24876) q[0];
sx q[0];
rz(-1.8565208) q[0];
rz(2.7934974) q[2];
sx q[2];
rz(-1.7762134) q[2];
sx q[2];
rz(1.5470488) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4917131) q[1];
sx q[1];
rz(-0.24398206) q[1];
sx q[1];
rz(0.53423832) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.25760381) q[3];
sx q[3];
rz(-1.8766878) q[3];
sx q[3];
rz(-2.5739939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.83071128) q[2];
sx q[2];
rz(-2.723366) q[2];
sx q[2];
rz(-2.6549784) q[2];
rz(-0.1085554) q[3];
sx q[3];
rz(-2.4311549) q[3];
sx q[3];
rz(1.9177115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9412823) q[0];
sx q[0];
rz(-2.0599984) q[0];
sx q[0];
rz(-2.6249028) q[0];
rz(1.5962881) q[1];
sx q[1];
rz(-1.0524806) q[1];
sx q[1];
rz(2.2470078) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45458083) q[0];
sx q[0];
rz(-2.5068388) q[0];
sx q[0];
rz(-1.4995585) q[0];
rz(1.9616227) q[2];
sx q[2];
rz(-0.50190364) q[2];
sx q[2];
rz(-2.6596136) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.48764187) q[1];
sx q[1];
rz(-0.36648053) q[1];
sx q[1];
rz(-1.3932863) q[1];
x q[2];
rz(0.92071269) q[3];
sx q[3];
rz(-1.4462785) q[3];
sx q[3];
rz(1.6684106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9937925) q[2];
sx q[2];
rz(-0.96269572) q[2];
sx q[2];
rz(-2.1195892) q[2];
rz(1.8036802) q[3];
sx q[3];
rz(-1.6513848) q[3];
sx q[3];
rz(-0.62906229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.548303) q[0];
sx q[0];
rz(-1.5240482) q[0];
sx q[0];
rz(-0.89292009) q[0];
rz(-2.100636) q[1];
sx q[1];
rz(-0.092408471) q[1];
sx q[1];
rz(1.668781) q[1];
rz(-0.54626089) q[2];
sx q[2];
rz(-2.1585474) q[2];
sx q[2];
rz(-1.3331158) q[2];
rz(1.7265504) q[3];
sx q[3];
rz(-2.8372436) q[3];
sx q[3];
rz(-2.5929034) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
