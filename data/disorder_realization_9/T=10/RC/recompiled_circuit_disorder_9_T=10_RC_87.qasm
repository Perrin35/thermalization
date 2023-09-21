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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.062382467) q[0];
sx q[0];
rz(-1.9431207) q[0];
sx q[0];
rz(-1.9553493) q[0];
rz(-2.652466) q[2];
sx q[2];
rz(-1.4406275) q[2];
sx q[2];
rz(2.7444073) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4218581) q[1];
sx q[1];
rz(-1.0654447) q[1];
sx q[1];
rz(-1.1769597) q[1];
rz(1.1236973) q[3];
sx q[3];
rz(-1.4456133) q[3];
sx q[3];
rz(3.0931635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.47544605) q[2];
sx q[2];
rz(-1.1300766) q[2];
sx q[2];
rz(-1.3249935) q[2];
rz(-2.825286) q[3];
sx q[3];
rz(-0.91747147) q[3];
sx q[3];
rz(-0.90707183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9544202) q[0];
sx q[0];
rz(-0.81048727) q[0];
sx q[0];
rz(2.1818838) q[0];
rz(-1.6540487) q[1];
sx q[1];
rz(-0.68351775) q[1];
sx q[1];
rz(2.5141292) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7146485) q[0];
sx q[0];
rz(-1.5639925) q[0];
sx q[0];
rz(0.0077008458) q[0];
rz(-pi) q[1];
rz(1.5487899) q[2];
sx q[2];
rz(-1.7823879) q[2];
sx q[2];
rz(1.5218658) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.97034772) q[1];
sx q[1];
rz(-1.8106318) q[1];
sx q[1];
rz(0.14648267) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.596132) q[3];
sx q[3];
rz(-1.5583894) q[3];
sx q[3];
rz(-2.5125463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.54005694) q[2];
sx q[2];
rz(-1.4011708) q[2];
sx q[2];
rz(-2.1616948) q[2];
rz(-1.4510138) q[3];
sx q[3];
rz(-2.4308949) q[3];
sx q[3];
rz(1.6987945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
rz(-2.5805805) q[0];
sx q[0];
rz(-0.0479527) q[0];
sx q[0];
rz(1.3797492) q[0];
rz(-0.17669949) q[1];
sx q[1];
rz(-1.9243762) q[1];
sx q[1];
rz(-1.7938991) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4189258) q[0];
sx q[0];
rz(-1.330266) q[0];
sx q[0];
rz(-2.6618883) q[0];
rz(1.548544) q[2];
sx q[2];
rz(-1.3705727) q[2];
sx q[2];
rz(0.74860886) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8417931) q[1];
sx q[1];
rz(-1.449297) q[1];
sx q[1];
rz(-1.4701162) q[1];
x q[2];
rz(-3.0511608) q[3];
sx q[3];
rz(-1.7548429) q[3];
sx q[3];
rz(-2.4726766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1220864) q[2];
sx q[2];
rz(-2.2821189) q[2];
sx q[2];
rz(-1.9879831) q[2];
rz(-0.38976088) q[3];
sx q[3];
rz(-2.6716158) q[3];
sx q[3];
rz(-0.78891689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47238123) q[0];
sx q[0];
rz(-0.62494576) q[0];
sx q[0];
rz(-0.56458449) q[0];
rz(-3.0584884) q[1];
sx q[1];
rz(-1.233498) q[1];
sx q[1];
rz(-0.21534236) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1911083) q[0];
sx q[0];
rz(-1.2380301) q[0];
sx q[0];
rz(-0.82975181) q[0];
rz(-0.16673659) q[2];
sx q[2];
rz(-2.1519289) q[2];
sx q[2];
rz(-0.8596479) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1064062) q[1];
sx q[1];
rz(-1.3760929) q[1];
sx q[1];
rz(-3.0343642) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4550291) q[3];
sx q[3];
rz(-0.6696223) q[3];
sx q[3];
rz(-1.7945822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.6354436) q[2];
sx q[2];
rz(-2.052867) q[2];
sx q[2];
rz(-0.60194683) q[2];
rz(0.69058949) q[3];
sx q[3];
rz(-0.75406051) q[3];
sx q[3];
rz(0.62600342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85201207) q[0];
sx q[0];
rz(-2.1471725) q[0];
sx q[0];
rz(2.1767298) q[0];
rz(0.016618641) q[1];
sx q[1];
rz(-0.42533541) q[1];
sx q[1];
rz(-0.66666493) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0913038) q[0];
sx q[0];
rz(-1.4698403) q[0];
sx q[0];
rz(1.9712649) q[0];
rz(-pi) q[1];
rz(-2.1140852) q[2];
sx q[2];
rz(-2.5230061) q[2];
sx q[2];
rz(1.0416043) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.9469229) q[1];
sx q[1];
rz(-1.8538845) q[1];
sx q[1];
rz(-2.4592295) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93249647) q[3];
sx q[3];
rz(-2.7780048) q[3];
sx q[3];
rz(-2.5225085) q[3];
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
rz(1.3145674) q[3];
sx q[3];
rz(-1.9606749) q[3];
sx q[3];
rz(-2.0107646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0261633) q[0];
sx q[0];
rz(-2.207343) q[0];
sx q[0];
rz(-1.7793659) q[0];
rz(1.7419787) q[1];
sx q[1];
rz(-1.5796966) q[1];
sx q[1];
rz(-1.5302352) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49028542) q[0];
sx q[0];
rz(-1.4978652) q[0];
sx q[0];
rz(-0.42071995) q[0];
rz(-pi) q[1];
x q[1];
rz(0.43892626) q[2];
sx q[2];
rz(-1.274684) q[2];
sx q[2];
rz(-2.211328) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.46576071) q[1];
sx q[1];
rz(-0.45468802) q[1];
sx q[1];
rz(2.0783706) q[1];
x q[2];
rz(1.3774032) q[3];
sx q[3];
rz(-2.8838727) q[3];
sx q[3];
rz(-1.6267488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7851012) q[2];
sx q[2];
rz(-2.1337528) q[2];
sx q[2];
rz(-3.1121758) q[2];
rz(-1.4908837) q[3];
sx q[3];
rz(-0.95508948) q[3];
sx q[3];
rz(-1.9194549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4239663) q[0];
sx q[0];
rz(-2.3762149) q[0];
sx q[0];
rz(-2.9587865) q[0];
rz(0.43552661) q[1];
sx q[1];
rz(-0.85414779) q[1];
sx q[1];
rz(0.51087728) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8683838) q[0];
sx q[0];
rz(-1.8330935) q[0];
sx q[0];
rz(1.6845076) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3235544) q[2];
sx q[2];
rz(-0.99070264) q[2];
sx q[2];
rz(1.8028508) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.70799815) q[1];
sx q[1];
rz(-2.3494548) q[1];
sx q[1];
rz(2.5670693) q[1];
rz(-pi) q[2];
rz(1.0501409) q[3];
sx q[3];
rz(-2.5556892) q[3];
sx q[3];
rz(-1.0019765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2304948) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34185103) q[0];
sx q[0];
rz(-0.36062476) q[0];
sx q[0];
rz(-2.0525232) q[0];
rz(-1.0808806) q[1];
sx q[1];
rz(-1.2842478) q[1];
sx q[1];
rz(0.63527766) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23112049) q[0];
sx q[0];
rz(-1.7251245) q[0];
sx q[0];
rz(-1.9327823) q[0];
x q[1];
rz(1.073457) q[2];
sx q[2];
rz(-1.4173696) q[2];
sx q[2];
rz(1.2129601) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.63102555) q[1];
sx q[1];
rz(-2.4706565) q[1];
sx q[1];
rz(-0.4354233) q[1];
rz(-pi) q[2];
rz(-1.2653207) q[3];
sx q[3];
rz(-2.5442122) q[3];
sx q[3];
rz(-1.7978096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.13715956) q[2];
sx q[2];
rz(-1.6798423) q[2];
sx q[2];
rz(0.00017246406) q[2];
rz(0.6862644) q[3];
sx q[3];
rz(-2.4133447) q[3];
sx q[3];
rz(2.2837158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2962608) q[0];
sx q[0];
rz(-2.2289472) q[0];
sx q[0];
rz(-0.6566748) q[0];
rz(0.36390057) q[1];
sx q[1];
rz(-2.6952422) q[1];
sx q[1];
rz(-2.231853) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1014935) q[0];
sx q[0];
rz(-1.3001406) q[0];
sx q[0];
rz(2.8069242) q[0];
rz(2.5932543) q[2];
sx q[2];
rz(-2.7395435) q[2];
sx q[2];
rz(-2.6532432) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.197387) q[1];
sx q[1];
rz(-1.3613609) q[1];
sx q[1];
rz(-1.6968813) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2551366) q[3];
sx q[3];
rz(-1.3254032) q[3];
sx q[3];
rz(1.0823702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.83071128) q[2];
sx q[2];
rz(-0.41822663) q[2];
sx q[2];
rz(0.48661423) q[2];
rz(-3.0330372) q[3];
sx q[3];
rz(-0.71043772) q[3];
sx q[3];
rz(-1.2238812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9412823) q[0];
sx q[0];
rz(-1.0815942) q[0];
sx q[0];
rz(0.51668984) q[0];
rz(-1.5962881) q[1];
sx q[1];
rz(-2.089112) q[1];
sx q[1];
rz(-0.89458481) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5986225) q[0];
sx q[0];
rz(-2.2036836) q[0];
sx q[0];
rz(3.0892239) q[0];
rz(-pi) q[1];
rz(-0.20608979) q[2];
sx q[2];
rz(-2.0317674) q[2];
sx q[2];
rz(3.0989625) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4641061) q[1];
sx q[1];
rz(-1.2103401) q[1];
sx q[1];
rz(-3.0739215) q[1];
x q[2];
rz(2.22088) q[3];
sx q[3];
rz(-1.6953141) q[3];
sx q[3];
rz(1.6684106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9937925) q[2];
sx q[2];
rz(-2.1788969) q[2];
sx q[2];
rz(-2.1195892) q[2];
rz(1.8036802) q[3];
sx q[3];
rz(-1.4902078) q[3];
sx q[3];
rz(-2.5125304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-0.90850496) q[2];
sx q[2];
rz(-2.0178595) q[2];
sx q[2];
rz(3.0541228) q[2];
rz(-1.8716807) q[3];
sx q[3];
rz(-1.5242929) q[3];
sx q[3];
rz(2.268189) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
