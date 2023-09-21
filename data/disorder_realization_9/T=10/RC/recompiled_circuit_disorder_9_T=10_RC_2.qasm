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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2406143) q[0];
sx q[0];
rz(-2.6128747) q[0];
sx q[0];
rz(-2.3762977) q[0];
rz(0.48912666) q[2];
sx q[2];
rz(-1.4406275) q[2];
sx q[2];
rz(-0.39718539) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0104048) q[1];
sx q[1];
rz(-0.63001761) q[1];
sx q[1];
rz(-0.6063993) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2875597) q[3];
sx q[3];
rz(-2.6784416) q[3];
sx q[3];
rz(-1.2676261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6661466) q[2];
sx q[2];
rz(-2.011516) q[2];
sx q[2];
rz(1.3249935) q[2];
rz(-0.31630668) q[3];
sx q[3];
rz(-0.91747147) q[3];
sx q[3];
rz(-2.2345208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9544202) q[0];
sx q[0];
rz(-2.3311054) q[0];
sx q[0];
rz(2.1818838) q[0];
rz(-1.487544) q[1];
sx q[1];
rz(-2.4580749) q[1];
sx q[1];
rz(-0.62746343) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42022959) q[0];
sx q[0];
rz(-3.1313167) q[0];
sx q[0];
rz(-0.72364877) q[0];
rz(-pi) q[1];
rz(1.5487899) q[2];
sx q[2];
rz(-1.3592048) q[2];
sx q[2];
rz(-1.5218658) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1712449) q[1];
sx q[1];
rz(-1.3309609) q[1];
sx q[1];
rz(0.14648267) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1153568) q[3];
sx q[3];
rz(-0.028209837) q[3];
sx q[3];
rz(1.3970323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6015357) q[2];
sx q[2];
rz(-1.4011708) q[2];
sx q[2];
rz(-2.1616948) q[2];
rz(-1.6905789) q[3];
sx q[3];
rz(-0.7106978) q[3];
sx q[3];
rz(-1.4427982) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56101218) q[0];
sx q[0];
rz(-0.0479527) q[0];
sx q[0];
rz(-1.7618435) q[0];
rz(-2.9648932) q[1];
sx q[1];
rz(-1.2172164) q[1];
sx q[1];
rz(1.3476936) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7226669) q[0];
sx q[0];
rz(-1.330266) q[0];
sx q[0];
rz(-2.6618883) q[0];
rz(-pi) q[1];
rz(-2.9413207) q[2];
sx q[2];
rz(-1.592604) q[2];
sx q[2];
rz(0.81776103) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.29979953) q[1];
sx q[1];
rz(-1.449297) q[1];
sx q[1];
rz(1.4701162) q[1];
rz(-pi) q[2];
rz(-1.3860116) q[3];
sx q[3];
rz(-1.6596969) q[3];
sx q[3];
rz(0.88528663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1220864) q[2];
sx q[2];
rz(-0.85947376) q[2];
sx q[2];
rz(1.1536095) q[2];
rz(0.38976088) q[3];
sx q[3];
rz(-2.6716158) q[3];
sx q[3];
rz(-2.3526758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47238123) q[0];
sx q[0];
rz(-2.5166469) q[0];
sx q[0];
rz(0.56458449) q[0];
rz(-0.083104221) q[1];
sx q[1];
rz(-1.233498) q[1];
sx q[1];
rz(-2.9262503) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4190061) q[0];
sx q[0];
rz(-2.3424087) q[0];
sx q[0];
rz(1.0976085) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98323804) q[2];
sx q[2];
rz(-1.7099656) q[2];
sx q[2];
rz(0.61901865) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6268057) q[1];
sx q[1];
rz(-1.4656015) q[1];
sx q[1];
rz(-1.766596) q[1];
rz(-2.5921949) q[3];
sx q[3];
rz(-1.975172) q[3];
sx q[3];
rz(0.79493633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.6354436) q[2];
sx q[2];
rz(-1.0887257) q[2];
sx q[2];
rz(0.60194683) q[2];
rz(2.4510032) q[3];
sx q[3];
rz(-0.75406051) q[3];
sx q[3];
rz(-0.62600342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85201207) q[0];
sx q[0];
rz(-2.1471725) q[0];
sx q[0];
rz(-2.1767298) q[0];
rz(3.124974) q[1];
sx q[1];
rz(-0.42533541) q[1];
sx q[1];
rz(0.66666493) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0502888) q[0];
sx q[0];
rz(-1.4698403) q[0];
sx q[0];
rz(-1.1703277) q[0];
rz(-pi) q[1];
rz(1.0275074) q[2];
sx q[2];
rz(-2.5230061) q[2];
sx q[2];
rz(-2.0999883) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7409089) q[1];
sx q[1];
rz(-0.92036696) q[1];
sx q[1];
rz(-1.9294192) q[1];
x q[2];
rz(-1.8673709) q[3];
sx q[3];
rz(-1.3572825) q[3];
sx q[3];
rz(-0.34542686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.91606402) q[2];
sx q[2];
rz(-1.3801882) q[2];
sx q[2];
rz(-1.2529681) q[2];
rz(-1.8270252) q[3];
sx q[3];
rz(-1.9606749) q[3];
sx q[3];
rz(1.130828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0261633) q[0];
sx q[0];
rz(-0.9342497) q[0];
sx q[0];
rz(1.3622267) q[0];
rz(1.399614) q[1];
sx q[1];
rz(-1.5618961) q[1];
sx q[1];
rz(-1.5302352) q[1];
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
rz(0.62263454) q[2];
sx q[2];
rz(-0.52402516) q[2];
sx q[2];
rz(-0.084409075) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6758319) q[1];
sx q[1];
rz(-2.6869046) q[1];
sx q[1];
rz(2.0783706) q[1];
rz(-pi) q[2];
rz(-1.3176765) q[3];
sx q[3];
rz(-1.5217921) q[3];
sx q[3];
rz(-0.24310902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7851012) q[2];
sx q[2];
rz(-1.0078398) q[2];
sx q[2];
rz(0.029416857) q[2];
rz(-1.4908837) q[3];
sx q[3];
rz(-0.95508948) q[3];
sx q[3];
rz(-1.9194549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
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
rz(-2.2874449) q[1];
sx q[1];
rz(2.6307154) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2832527) q[0];
sx q[0];
rz(-2.8562299) q[0];
sx q[0];
rz(0.39985379) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80656959) q[2];
sx q[2];
rz(-0.91420805) q[2];
sx q[2];
rz(-0.29733959) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7055197) q[1];
sx q[1];
rz(-1.9680068) q[1];
sx q[1];
rz(-2.436609) q[1];
rz(-pi) q[2];
rz(1.0501409) q[3];
sx q[3];
rz(-2.5556892) q[3];
sx q[3];
rz(-1.0019765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2304948) q[2];
sx q[2];
rz(-2.348933) q[2];
sx q[2];
rz(2.5210209) q[2];
rz(-2.7190322) q[3];
sx q[3];
rz(-1.3603323) q[3];
sx q[3];
rz(-1.871199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34185103) q[0];
sx q[0];
rz(-2.7809679) q[0];
sx q[0];
rz(2.0525232) q[0];
rz(-2.0607121) q[1];
sx q[1];
rz(-1.2842478) q[1];
sx q[1];
rz(2.506315) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23112049) q[0];
sx q[0];
rz(-1.7251245) q[0];
sx q[0];
rz(-1.9327823) q[0];
rz(2.0681357) q[2];
sx q[2];
rz(-1.4173696) q[2];
sx q[2];
rz(1.9286326) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0465225) q[1];
sx q[1];
rz(-0.97192837) q[1];
sx q[1];
rz(-1.2477161) q[1];
rz(-pi) q[2];
rz(0.99526309) q[3];
sx q[3];
rz(-1.7407773) q[3];
sx q[3];
rz(3.1135524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.13715956) q[2];
sx q[2];
rz(-1.6798423) q[2];
sx q[2];
rz(-3.1414202) q[2];
rz(2.4553283) q[3];
sx q[3];
rz(-2.4133447) q[3];
sx q[3];
rz(0.85787684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2962608) q[0];
sx q[0];
rz(-2.2289472) q[0];
sx q[0];
rz(-2.4849179) q[0];
rz(-0.36390057) q[1];
sx q[1];
rz(-2.6952422) q[1];
sx q[1];
rz(2.231853) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1014935) q[0];
sx q[0];
rz(-1.3001406) q[0];
sx q[0];
rz(-0.33466848) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5932543) q[2];
sx q[2];
rz(-2.7395435) q[2];
sx q[2];
rz(-0.48834947) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.39975702) q[1];
sx q[1];
rz(-1.6941119) q[1];
sx q[1];
rz(0.21106212) q[1];
rz(-pi) q[2];
rz(-0.25760381) q[3];
sx q[3];
rz(-1.8766878) q[3];
sx q[3];
rz(-2.5739939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9412823) q[0];
sx q[0];
rz(-1.0815942) q[0];
sx q[0];
rz(2.6249028) q[0];
rz(1.5453045) q[1];
sx q[1];
rz(-1.0524806) q[1];
sx q[1];
rz(0.89458481) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45458083) q[0];
sx q[0];
rz(-2.5068388) q[0];
sx q[0];
rz(1.6420341) q[0];
x q[1];
rz(-1.1799699) q[2];
sx q[2];
rz(-0.50190364) q[2];
sx q[2];
rz(-2.6596136) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4641061) q[1];
sx q[1];
rz(-1.2103401) q[1];
sx q[1];
rz(3.0739215) q[1];
rz(1.3668725) q[3];
sx q[3];
rz(-0.66019928) q[3];
sx q[3];
rz(-0.064299671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.14780012) q[2];
sx q[2];
rz(-0.96269572) q[2];
sx q[2];
rz(-1.0220035) q[2];
rz(1.8036802) q[3];
sx q[3];
rz(-1.6513848) q[3];
sx q[3];
rz(-0.62906229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-1.548303) q[0];
sx q[0];
rz(-1.6175445) q[0];
sx q[0];
rz(2.2486726) q[0];
rz(-1.0409566) q[1];
sx q[1];
rz(-3.0491842) q[1];
sx q[1];
rz(-1.4728117) q[1];
rz(-0.54626089) q[2];
sx q[2];
rz(-2.1585474) q[2];
sx q[2];
rz(-1.3331158) q[2];
rz(-0.048687497) q[3];
sx q[3];
rz(-1.8713453) q[3];
sx q[3];
rz(0.71181675) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];