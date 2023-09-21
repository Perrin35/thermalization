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
rz(2.467632) q[0];
rz(-0.31710467) q[1];
sx q[1];
rz(-1.6333406) q[1];
sx q[1];
rz(-0.79467264) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90097839) q[0];
sx q[0];
rz(-0.52871791) q[0];
sx q[0];
rz(0.76529495) q[0];
rz(-pi) q[1];
rz(2.8698679) q[2];
sx q[2];
rz(-2.6368015) q[2];
sx q[2];
rz(1.4128078) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7197345) q[1];
sx q[1];
rz(-2.0761479) q[1];
sx q[1];
rz(1.964633) q[1];
x q[2];
rz(-1.854033) q[3];
sx q[3];
rz(-2.6784416) q[3];
sx q[3];
rz(1.8739665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6661466) q[2];
sx q[2];
rz(-1.1300766) q[2];
sx q[2];
rz(-1.3249935) q[2];
rz(-0.31630668) q[3];
sx q[3];
rz(-0.91747147) q[3];
sx q[3];
rz(-2.2345208) q[3];
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
rz(-1.1871724) q[0];
sx q[0];
rz(-2.3311054) q[0];
sx q[0];
rz(0.95970884) q[0];
rz(-1.487544) q[1];
sx q[1];
rz(-2.4580749) q[1];
sx q[1];
rz(2.5141292) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7146485) q[0];
sx q[0];
rz(-1.5639925) q[0];
sx q[0];
rz(-0.0077008458) q[0];
rz(-0.10208315) q[2];
sx q[2];
rz(-2.9288769) q[2];
sx q[2];
rz(1.724147) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6154502) q[1];
sx q[1];
rz(-2.8612988) q[1];
sx q[1];
rz(-1.0326833) q[1];
x q[2];
rz(-1.596132) q[3];
sx q[3];
rz(-1.5832033) q[3];
sx q[3];
rz(2.5125463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.54005694) q[2];
sx q[2];
rz(-1.7404218) q[2];
sx q[2];
rz(2.1616948) q[2];
rz(1.4510138) q[3];
sx q[3];
rz(-2.4308949) q[3];
sx q[3];
rz(1.4427982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56101218) q[0];
sx q[0];
rz(-3.09364) q[0];
sx q[0];
rz(-1.3797492) q[0];
rz(-0.17669949) q[1];
sx q[1];
rz(-1.2172164) q[1];
sx q[1];
rz(-1.3476936) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1130226) q[0];
sx q[0];
rz(-1.1060113) q[0];
sx q[0];
rz(1.3010498) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0323896) q[2];
sx q[2];
rz(-0.20143992) q[2];
sx q[2];
rz(-0.86004721) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2587535) q[1];
sx q[1];
rz(-1.6707318) q[1];
sx q[1];
rz(-0.12211166) q[1];
rz(-pi) q[2];
rz(1.1191145) q[3];
sx q[3];
rz(-2.936755) q[3];
sx q[3];
rz(1.1289489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0195062) q[2];
sx q[2];
rz(-2.2821189) q[2];
sx q[2];
rz(-1.9879831) q[2];
rz(2.7518318) q[3];
sx q[3];
rz(-2.6716158) q[3];
sx q[3];
rz(2.3526758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47238123) q[0];
sx q[0];
rz(-2.5166469) q[0];
sx q[0];
rz(-0.56458449) q[0];
rz(3.0584884) q[1];
sx q[1];
rz(-1.9080947) q[1];
sx q[1];
rz(-0.21534236) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1911083) q[0];
sx q[0];
rz(-1.9035625) q[0];
sx q[0];
rz(-2.3118408) q[0];
x q[1];
rz(-2.1583546) q[2];
sx q[2];
rz(-1.431627) q[2];
sx q[2];
rz(-2.522574) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1064062) q[1];
sx q[1];
rz(-1.7654997) q[1];
sx q[1];
rz(0.10722843) q[1];
rz(-pi) q[2];
rz(-2.0358884) q[3];
sx q[3];
rz(-2.0715189) q[3];
sx q[3];
rz(-0.53946686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5061491) q[2];
sx q[2];
rz(-2.052867) q[2];
sx q[2];
rz(2.5396458) q[2];
rz(2.4510032) q[3];
sx q[3];
rz(-2.3875321) q[3];
sx q[3];
rz(-2.5155892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2895806) q[0];
sx q[0];
rz(-2.1471725) q[0];
sx q[0];
rz(-0.96486282) q[0];
rz(3.124974) q[1];
sx q[1];
rz(-2.7162572) q[1];
sx q[1];
rz(2.4749277) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0913038) q[0];
sx q[0];
rz(-1.6717523) q[0];
sx q[0];
rz(1.1703277) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0275074) q[2];
sx q[2];
rz(-2.5230061) q[2];
sx q[2];
rz(-2.0999883) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.9469229) q[1];
sx q[1];
rz(-1.8538845) q[1];
sx q[1];
rz(-0.6823632) q[1];
rz(-pi) q[2];
rz(1.2742217) q[3];
sx q[3];
rz(-1.7843102) q[3];
sx q[3];
rz(0.34542686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2255286) q[2];
sx q[2];
rz(-1.7614044) q[2];
sx q[2];
rz(-1.8886245) q[2];
rz(1.8270252) q[3];
sx q[3];
rz(-1.1809177) q[3];
sx q[3];
rz(1.130828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0261633) q[0];
sx q[0];
rz(-0.9342497) q[0];
sx q[0];
rz(-1.7793659) q[0];
rz(-1.399614) q[1];
sx q[1];
rz(-1.5618961) q[1];
sx q[1];
rz(1.5302352) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49028542) q[0];
sx q[0];
rz(-1.4978652) q[0];
sx q[0];
rz(-2.7208727) q[0];
rz(-pi) q[1];
rz(-0.43892626) q[2];
sx q[2];
rz(-1.8669087) q[2];
sx q[2];
rz(-2.211328) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1215026) q[1];
sx q[1];
rz(-1.964718) q[1];
sx q[1];
rz(0.23328383) q[1];
x q[2];
rz(-1.3774032) q[3];
sx q[3];
rz(-2.8838727) q[3];
sx q[3];
rz(1.6267488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.35649148) q[2];
sx q[2];
rz(-2.1337528) q[2];
sx q[2];
rz(3.1121758) q[2];
rz(1.4908837) q[3];
sx q[3];
rz(-0.95508948) q[3];
sx q[3];
rz(1.9194549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(2.4239663) q[0];
sx q[0];
rz(-0.76537776) q[0];
sx q[0];
rz(-2.9587865) q[0];
rz(-0.43552661) q[1];
sx q[1];
rz(-2.2874449) q[1];
sx q[1];
rz(0.51087728) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2832527) q[0];
sx q[0];
rz(-2.8562299) q[0];
sx q[0];
rz(-0.39985379) q[0];
rz(-pi) q[1];
rz(-2.4099318) q[2];
sx q[2];
rz(-2.1795142) q[2];
sx q[2];
rz(-0.70639709) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.70799815) q[1];
sx q[1];
rz(-2.3494548) q[1];
sx q[1];
rz(-0.57452332) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0501409) q[3];
sx q[3];
rz(-0.58590349) q[3];
sx q[3];
rz(-1.0019765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9110979) q[2];
sx q[2];
rz(-0.7926597) q[2];
sx q[2];
rz(-2.5210209) q[2];
rz(-2.7190322) q[3];
sx q[3];
rz(-1.7812604) q[3];
sx q[3];
rz(-1.2703936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7997416) q[0];
sx q[0];
rz(-2.7809679) q[0];
sx q[0];
rz(-2.0525232) q[0];
rz(-1.0808806) q[1];
sx q[1];
rz(-1.8573449) q[1];
sx q[1];
rz(-0.63527766) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9104722) q[0];
sx q[0];
rz(-1.4164682) q[0];
sx q[0];
rz(1.2088103) q[0];
rz(-pi) q[1];
rz(-2.0681357) q[2];
sx q[2];
rz(-1.724223) q[2];
sx q[2];
rz(-1.2129601) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.63102555) q[1];
sx q[1];
rz(-0.67093611) q[1];
sx q[1];
rz(2.7061694) q[1];
x q[2];
rz(-0.99526309) q[3];
sx q[3];
rz(-1.7407773) q[3];
sx q[3];
rz(-3.1135524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0044331) q[2];
sx q[2];
rz(-1.6798423) q[2];
sx q[2];
rz(0.00017246406) q[2];
rz(0.6862644) q[3];
sx q[3];
rz(-2.4133447) q[3];
sx q[3];
rz(-0.85787684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2962608) q[0];
sx q[0];
rz(-2.2289472) q[0];
sx q[0];
rz(2.4849179) q[0];
rz(-0.36390057) q[1];
sx q[1];
rz(-2.6952422) q[1];
sx q[1];
rz(-0.90973967) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1014935) q[0];
sx q[0];
rz(-1.3001406) q[0];
sx q[0];
rz(0.33466848) q[0];
rz(-pi) q[1];
rz(-1.7889195) q[2];
sx q[2];
rz(-1.2303196) q[2];
sx q[2];
rz(3.043963) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4917131) q[1];
sx q[1];
rz(-2.8976106) q[1];
sx q[1];
rz(0.53423832) q[1];
rz(-pi) q[2];
rz(-2.8839888) q[3];
sx q[3];
rz(-1.2649049) q[3];
sx q[3];
rz(0.56759873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3108814) q[2];
sx q[2];
rz(-2.723366) q[2];
sx q[2];
rz(2.6549784) q[2];
rz(3.0330372) q[3];
sx q[3];
rz(-0.71043772) q[3];
sx q[3];
rz(-1.9177115) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20031032) q[0];
sx q[0];
rz(-1.0815942) q[0];
sx q[0];
rz(-2.6249028) q[0];
rz(1.5453045) q[1];
sx q[1];
rz(-1.0524806) q[1];
sx q[1];
rz(0.89458481) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45458083) q[0];
sx q[0];
rz(-2.5068388) q[0];
sx q[0];
rz(-1.6420341) q[0];
rz(1.1012494) q[2];
sx q[2];
rz(-1.38648) q[2];
sx q[2];
rz(1.4354401) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6539508) q[1];
sx q[1];
rz(-0.36648053) q[1];
sx q[1];
rz(1.7483064) q[1];
rz(-0.92071269) q[3];
sx q[3];
rz(-1.4462785) q[3];
sx q[3];
rz(-1.6684106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9937925) q[2];
sx q[2];
rz(-0.96269572) q[2];
sx q[2];
rz(2.1195892) q[2];
rz(1.3379124) q[3];
sx q[3];
rz(-1.6513848) q[3];
sx q[3];
rz(0.62906229) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.548303) q[0];
sx q[0];
rz(-1.5240482) q[0];
sx q[0];
rz(-0.89292009) q[0];
rz(2.100636) q[1];
sx q[1];
rz(-3.0491842) q[1];
sx q[1];
rz(-1.4728117) q[1];
rz(-2.5953318) q[2];
sx q[2];
rz(-0.98304521) q[2];
sx q[2];
rz(1.8084768) q[2];
rz(1.2699119) q[3];
sx q[3];
rz(-1.5242929) q[3];
sx q[3];
rz(2.268189) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];