OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.5716612) q[0];
sx q[0];
rz(-0.067458955) q[0];
sx q[0];
rz(10.098739) q[0];
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
rz(0.062382467) q[0];
sx q[0];
rz(-1.9431207) q[0];
sx q[0];
rz(1.9553493) q[0];
x q[1];
rz(0.48912666) q[2];
sx q[2];
rz(-1.4406275) q[2];
sx q[2];
rz(-0.39718539) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0104048) q[1];
sx q[1];
rz(-0.63001761) q[1];
sx q[1];
rz(0.6063993) q[1];
rz(-pi) q[2];
x q[2];
rz(1.854033) q[3];
sx q[3];
rz(-2.6784416) q[3];
sx q[3];
rz(-1.8739665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6661466) q[2];
sx q[2];
rz(-2.011516) q[2];
sx q[2];
rz(-1.3249935) q[2];
rz(0.31630668) q[3];
sx q[3];
rz(-0.91747147) q[3];
sx q[3];
rz(-0.90707183) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9544202) q[0];
sx q[0];
rz(-0.81048727) q[0];
sx q[0];
rz(-0.95970884) q[0];
rz(-1.487544) q[1];
sx q[1];
rz(-0.68351775) q[1];
sx q[1];
rz(0.62746343) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42694416) q[0];
sx q[0];
rz(-1.5776002) q[0];
sx q[0];
rz(0.0077008458) q[0];
rz(-pi) q[1];
x q[1];
rz(0.10208315) q[2];
sx q[2];
rz(-2.9288769) q[2];
sx q[2];
rz(-1.724147) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6154502) q[1];
sx q[1];
rz(-0.28029385) q[1];
sx q[1];
rz(-1.0326833) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1291817) q[3];
sx q[3];
rz(-1.5454626) q[3];
sx q[3];
rz(2.2001571) q[3];
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
rz(-0.97989782) q[2];
rz(-1.4510138) q[3];
sx q[3];
rz(-2.4308949) q[3];
sx q[3];
rz(1.6987945) q[3];
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
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5805805) q[0];
sx q[0];
rz(-0.0479527) q[0];
sx q[0];
rz(1.7618435) q[0];
rz(0.17669949) q[1];
sx q[1];
rz(-1.9243762) q[1];
sx q[1];
rz(1.7938991) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028570024) q[0];
sx q[0];
rz(-2.0355814) q[0];
sx q[0];
rz(-1.3010498) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.548544) q[2];
sx q[2];
rz(-1.3705727) q[2];
sx q[2];
rz(2.3929838) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8828391) q[1];
sx q[1];
rz(-1.4708609) q[1];
sx q[1];
rz(-3.019481) q[1];
rz(-pi) q[2];
rz(2.0224781) q[3];
sx q[3];
rz(-0.20483769) q[3];
sx q[3];
rz(-2.0126437) q[3];
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
rz(-0.38976088) q[3];
sx q[3];
rz(-0.4699769) q[3];
sx q[3];
rz(-2.3526758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47238123) q[0];
sx q[0];
rz(-2.5166469) q[0];
sx q[0];
rz(-0.56458449) q[0];
rz(0.083104221) q[1];
sx q[1];
rz(-1.9080947) q[1];
sx q[1];
rz(-2.9262503) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1911083) q[0];
sx q[0];
rz(-1.2380301) q[0];
sx q[0];
rz(-2.3118408) q[0];
x q[1];
rz(-0.98323804) q[2];
sx q[2];
rz(-1.7099656) q[2];
sx q[2];
rz(-2.522574) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1064062) q[1];
sx q[1];
rz(-1.3760929) q[1];
sx q[1];
rz(-3.0343642) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68656355) q[3];
sx q[3];
rz(-2.4719704) q[3];
sx q[3];
rz(1.7945822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.6354436) q[2];
sx q[2];
rz(-1.0887257) q[2];
sx q[2];
rz(2.5396458) q[2];
rz(2.4510032) q[3];
sx q[3];
rz(-0.75406051) q[3];
sx q[3];
rz(-0.62600342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85201207) q[0];
sx q[0];
rz(-2.1471725) q[0];
sx q[0];
rz(-2.1767298) q[0];
rz(-0.016618641) q[1];
sx q[1];
rz(-0.42533541) q[1];
sx q[1];
rz(0.66666493) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43685164) q[0];
sx q[0];
rz(-1.1724823) q[0];
sx q[0];
rz(-0.10956357) q[0];
rz(1.0275074) q[2];
sx q[2];
rz(-0.61858656) q[2];
sx q[2];
rz(2.0999883) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.9469229) q[1];
sx q[1];
rz(-1.2877081) q[1];
sx q[1];
rz(-0.6823632) q[1];
x q[2];
rz(-0.93249647) q[3];
sx q[3];
rz(-2.7780048) q[3];
sx q[3];
rz(0.61908412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2255286) q[2];
sx q[2];
rz(-1.7614044) q[2];
sx q[2];
rz(-1.2529681) q[2];
rz(-1.8270252) q[3];
sx q[3];
rz(-1.1809177) q[3];
sx q[3];
rz(-1.130828) q[3];
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
rz(-1.1154293) q[0];
sx q[0];
rz(-0.9342497) q[0];
sx q[0];
rz(1.3622267) q[0];
rz(1.7419787) q[1];
sx q[1];
rz(-1.5796966) q[1];
sx q[1];
rz(1.6113575) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2419469) q[0];
sx q[0];
rz(-0.42662222) q[0];
sx q[0];
rz(-0.17701478) q[0];
rz(1.8958695) q[2];
sx q[2];
rz(-1.9893861) q[2];
sx q[2];
rz(-0.77667728) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1215026) q[1];
sx q[1];
rz(-1.964718) q[1];
sx q[1];
rz(-0.23328383) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8239162) q[3];
sx q[3];
rz(-1.6198006) q[3];
sx q[3];
rz(0.24310902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7851012) q[2];
sx q[2];
rz(-2.1337528) q[2];
sx q[2];
rz(-0.029416857) q[2];
rz(1.4908837) q[3];
sx q[3];
rz(-0.95508948) q[3];
sx q[3];
rz(1.9194549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4239663) q[0];
sx q[0];
rz(-2.3762149) q[0];
sx q[0];
rz(0.18280612) q[0];
rz(0.43552661) q[1];
sx q[1];
rz(-2.2874449) q[1];
sx q[1];
rz(-0.51087728) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2832527) q[0];
sx q[0];
rz(-0.28536277) q[0];
sx q[0];
rz(-0.39985379) q[0];
x q[1];
rz(2.3350231) q[2];
sx q[2];
rz(-0.91420805) q[2];
sx q[2];
rz(-2.8442531) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7055197) q[1];
sx q[1];
rz(-1.9680068) q[1];
sx q[1];
rz(-0.7049837) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0914518) q[3];
sx q[3];
rz(-0.58590349) q[3];
sx q[3];
rz(1.0019765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2304948) q[2];
sx q[2];
rz(-2.348933) q[2];
sx q[2];
rz(0.62057173) q[2];
rz(2.7190322) q[3];
sx q[3];
rz(-1.3603323) q[3];
sx q[3];
rz(1.871199) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7997416) q[0];
sx q[0];
rz(-0.36062476) q[0];
sx q[0];
rz(-2.0525232) q[0];
rz(-2.0607121) q[1];
sx q[1];
rz(-1.8573449) q[1];
sx q[1];
rz(-2.506315) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1875293) q[0];
sx q[0];
rz(-2.7494193) q[0];
sx q[0];
rz(1.1568882) q[0];
rz(-0.17417553) q[2];
sx q[2];
rz(-2.0617699) q[2];
sx q[2];
rz(2.8665286) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5105671) q[1];
sx q[1];
rz(-2.4706565) q[1];
sx q[1];
rz(0.4354233) q[1];
rz(1.876272) q[3];
sx q[3];
rz(-0.59738041) q[3];
sx q[3];
rz(-1.3437831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.13715956) q[2];
sx q[2];
rz(-1.4617504) q[2];
sx q[2];
rz(0.00017246406) q[2];
rz(-0.6862644) q[3];
sx q[3];
rz(-2.4133447) q[3];
sx q[3];
rz(0.85787684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8453318) q[0];
sx q[0];
rz(-2.2289472) q[0];
sx q[0];
rz(0.6566748) q[0];
rz(-0.36390057) q[1];
sx q[1];
rz(-0.44635043) q[1];
sx q[1];
rz(-2.231853) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0168403) q[0];
sx q[0];
rz(-2.7144055) q[0];
sx q[0];
rz(2.4401526) q[0];
rz(-pi) q[1];
rz(-1.3526731) q[2];
sx q[2];
rz(-1.2303196) q[2];
sx q[2];
rz(-3.043963) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7418356) q[1];
sx q[1];
rz(-1.6941119) q[1];
sx q[1];
rz(0.21106212) q[1];
x q[2];
rz(1.886456) q[3];
sx q[3];
rz(-1.8161895) q[3];
sx q[3];
rz(2.0592225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.83071128) q[2];
sx q[2];
rz(-2.723366) q[2];
sx q[2];
rz(-0.48661423) q[2];
rz(-0.1085554) q[3];
sx q[3];
rz(-0.71043772) q[3];
sx q[3];
rz(-1.9177115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9412823) q[0];
sx q[0];
rz(-2.0599984) q[0];
sx q[0];
rz(-0.51668984) q[0];
rz(1.5962881) q[1];
sx q[1];
rz(-1.0524806) q[1];
sx q[1];
rz(2.2470078) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5986225) q[0];
sx q[0];
rz(-0.93790903) q[0];
sx q[0];
rz(-0.052368725) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9355029) q[2];
sx q[2];
rz(-1.1098252) q[2];
sx q[2];
rz(-0.042630171) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.91720944) q[1];
sx q[1];
rz(-1.50748) q[1];
sx q[1];
rz(-1.2095832) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7747202) q[3];
sx q[3];
rz(-0.66019928) q[3];
sx q[3];
rz(-3.077293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.14780012) q[2];
sx q[2];
rz(-0.96269572) q[2];
sx q[2];
rz(-2.1195892) q[2];
rz(-1.8036802) q[3];
sx q[3];
rz(-1.6513848) q[3];
sx q[3];
rz(-2.5125304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.548303) q[0];
sx q[0];
rz(-1.6175445) q[0];
sx q[0];
rz(2.2486726) q[0];
rz(-2.100636) q[1];
sx q[1];
rz(-0.092408471) q[1];
sx q[1];
rz(1.668781) q[1];
rz(-0.9085761) q[2];
sx q[2];
rz(-2.3618345) q[2];
sx q[2];
rz(0.97710412) q[2];
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
