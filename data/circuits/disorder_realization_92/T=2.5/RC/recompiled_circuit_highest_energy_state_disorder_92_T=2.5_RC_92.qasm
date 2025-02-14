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
rz(0.25843698) q[0];
sx q[0];
rz(-0.28764495) q[0];
sx q[0];
rz(2.0339461) q[0];
rz(-1.8183174) q[1];
sx q[1];
rz(-0.32185093) q[1];
sx q[1];
rz(0.62825656) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0387255) q[0];
sx q[0];
rz(-2.9280781) q[0];
sx q[0];
rz(-0.98189129) q[0];
x q[1];
rz(1.1789178) q[2];
sx q[2];
rz(-0.88681839) q[2];
sx q[2];
rz(-0.98525233) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9440556) q[1];
sx q[1];
rz(-1.5345795) q[1];
sx q[1];
rz(0.22878583) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7006458) q[3];
sx q[3];
rz(-1.0797056) q[3];
sx q[3];
rz(0.90346149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.85311741) q[2];
sx q[2];
rz(-1.709781) q[2];
sx q[2];
rz(-0.34118578) q[2];
rz(0.8902542) q[3];
sx q[3];
rz(-1.2672) q[3];
sx q[3];
rz(-2.6064579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5219236) q[0];
sx q[0];
rz(-0.6837908) q[0];
sx q[0];
rz(2.4271915) q[0];
rz(1.5996541) q[1];
sx q[1];
rz(-2.6456412) q[1];
sx q[1];
rz(-3.0468859) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45825645) q[0];
sx q[0];
rz(-0.33412513) q[0];
sx q[0];
rz(-1.279356) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33659597) q[2];
sx q[2];
rz(-0.99592745) q[2];
sx q[2];
rz(2.6624555) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.75443711) q[1];
sx q[1];
rz(-1.2278986) q[1];
sx q[1];
rz(1.8217525) q[1];
rz(-pi) q[2];
rz(0.27908947) q[3];
sx q[3];
rz(-0.52625873) q[3];
sx q[3];
rz(-2.581439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.16227214) q[2];
sx q[2];
rz(-1.2459735) q[2];
sx q[2];
rz(-2.6675513) q[2];
rz(0.76215172) q[3];
sx q[3];
rz(-2.5049329) q[3];
sx q[3];
rz(0.4711841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.8629465) q[0];
sx q[0];
rz(-0.11676783) q[0];
sx q[0];
rz(-2.2886724) q[0];
rz(0.22311738) q[1];
sx q[1];
rz(-2.6631963) q[1];
sx q[1];
rz(2.6257637) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0313161) q[0];
sx q[0];
rz(-3.0678684) q[0];
sx q[0];
rz(2.323193) q[0];
rz(2.8858809) q[2];
sx q[2];
rz(-2.1923794) q[2];
sx q[2];
rz(-2.5598124) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8341537) q[1];
sx q[1];
rz(-0.98886988) q[1];
sx q[1];
rz(-0.80640275) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7985057) q[3];
sx q[3];
rz(-2.2394387) q[3];
sx q[3];
rz(2.9647788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8082661) q[2];
sx q[2];
rz(-2.1648679) q[2];
sx q[2];
rz(0.094956368) q[2];
rz(0.44148463) q[3];
sx q[3];
rz(-1.35651) q[3];
sx q[3];
rz(-1.1536185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3917291) q[0];
sx q[0];
rz(-0.57498217) q[0];
sx q[0];
rz(1.0561426) q[0];
rz(-2.1978343) q[1];
sx q[1];
rz(-2.8902003) q[1];
sx q[1];
rz(-1.6897078) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3849996) q[0];
sx q[0];
rz(-0.90314048) q[0];
sx q[0];
rz(-0.7636015) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9610434) q[2];
sx q[2];
rz(-2.0023291) q[2];
sx q[2];
rz(0.081588946) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.31409594) q[1];
sx q[1];
rz(-0.24330595) q[1];
sx q[1];
rz(-2.9929586) q[1];
rz(-pi) q[2];
x q[2];
rz(0.068030997) q[3];
sx q[3];
rz(-3.0306731) q[3];
sx q[3];
rz(2.7601931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4132495) q[2];
sx q[2];
rz(-2.3409833) q[2];
sx q[2];
rz(-2.441414) q[2];
rz(-0.87107301) q[3];
sx q[3];
rz(-2.0756523) q[3];
sx q[3];
rz(1.357249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(1.0102608) q[0];
sx q[0];
rz(-2.6483674) q[0];
sx q[0];
rz(-2.6569271) q[0];
rz(-0.87801814) q[1];
sx q[1];
rz(-1.9819825) q[1];
sx q[1];
rz(2.7427618) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53610401) q[0];
sx q[0];
rz(-1.4358504) q[0];
sx q[0];
rz(0.55079726) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2315751) q[2];
sx q[2];
rz(-0.89618635) q[2];
sx q[2];
rz(2.2774027) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.63020413) q[1];
sx q[1];
rz(-2.2126019) q[1];
sx q[1];
rz(-0.46496572) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5957803) q[3];
sx q[3];
rz(-2.1936676) q[3];
sx q[3];
rz(2.2198077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.47348076) q[2];
sx q[2];
rz(-1.6010189) q[2];
sx q[2];
rz(-0.87781805) q[2];
rz(0.38464883) q[3];
sx q[3];
rz(-0.84482241) q[3];
sx q[3];
rz(1.1788684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8574852) q[0];
sx q[0];
rz(-0.50898886) q[0];
sx q[0];
rz(0.2704764) q[0];
rz(-2.8349304) q[1];
sx q[1];
rz(-0.51391927) q[1];
sx q[1];
rz(-0.57672966) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8054415) q[0];
sx q[0];
rz(-2.9362943) q[0];
sx q[0];
rz(-2.0144256) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36251601) q[2];
sx q[2];
rz(-0.91372817) q[2];
sx q[2];
rz(2.9566979) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8851848) q[1];
sx q[1];
rz(-1.9659623) q[1];
sx q[1];
rz(-0.97104071) q[1];
rz(-pi) q[2];
rz(2.8068325) q[3];
sx q[3];
rz(-1.5840414) q[3];
sx q[3];
rz(-2.5760026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3660672) q[2];
sx q[2];
rz(-1.0696573) q[2];
sx q[2];
rz(-0.97037399) q[2];
rz(-0.35106418) q[3];
sx q[3];
rz(-0.23215663) q[3];
sx q[3];
rz(-2.8620913) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.905726) q[0];
sx q[0];
rz(-2.7577363) q[0];
sx q[0];
rz(0.50091499) q[0];
rz(2.4210335) q[1];
sx q[1];
rz(-0.84545285) q[1];
sx q[1];
rz(0.93773425) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2116579) q[0];
sx q[0];
rz(-1.6119566) q[0];
sx q[0];
rz(1.3672872) q[0];
rz(-pi) q[1];
rz(-2.4831122) q[2];
sx q[2];
rz(-2.5665433) q[2];
sx q[2];
rz(0.63840309) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8775505) q[1];
sx q[1];
rz(-0.94549798) q[1];
sx q[1];
rz(2.6912415) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.012091919) q[3];
sx q[3];
rz(-0.68688697) q[3];
sx q[3];
rz(1.6364607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.131989) q[2];
sx q[2];
rz(-0.14543532) q[2];
sx q[2];
rz(-2.2930938) q[2];
rz(-2.2961473) q[3];
sx q[3];
rz(-2.4852821) q[3];
sx q[3];
rz(-1.9158624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.41035143) q[0];
sx q[0];
rz(-2.1827965) q[0];
sx q[0];
rz(2.248726) q[0];
rz(-3.0800173) q[1];
sx q[1];
rz(-0.67191809) q[1];
sx q[1];
rz(0.16199131) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50349405) q[0];
sx q[0];
rz(-1.5862521) q[0];
sx q[0];
rz(1.6729805) q[0];
rz(2.4396877) q[2];
sx q[2];
rz(-1.4830952) q[2];
sx q[2];
rz(0.039393124) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.87016856) q[1];
sx q[1];
rz(-1.7080016) q[1];
sx q[1];
rz(2.8297564) q[1];
rz(-0.83655436) q[3];
sx q[3];
rz(-3.0020107) q[3];
sx q[3];
rz(-2.9849646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.091247678) q[2];
sx q[2];
rz(-2.7140736) q[2];
sx q[2];
rz(-2.4166935) q[2];
rz(-2.8643518) q[3];
sx q[3];
rz(-0.96479779) q[3];
sx q[3];
rz(3.0550756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8898833) q[0];
sx q[0];
rz(-2.4885663) q[0];
sx q[0];
rz(-0.025803331) q[0];
rz(1.7426527) q[1];
sx q[1];
rz(-1.502123) q[1];
sx q[1];
rz(0.10765156) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.066641971) q[0];
sx q[0];
rz(-2.5558976) q[0];
sx q[0];
rz(-1.8138712) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4716267) q[2];
sx q[2];
rz(-1.1749069) q[2];
sx q[2];
rz(0.63360032) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4835158) q[1];
sx q[1];
rz(-0.47464322) q[1];
sx q[1];
rz(-3.062611) q[1];
x q[2];
rz(-3.1237431) q[3];
sx q[3];
rz(-1.4571725) q[3];
sx q[3];
rz(-0.85434948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4699576) q[2];
sx q[2];
rz(-1.8822972) q[2];
sx q[2];
rz(-2.0476445) q[2];
rz(2.5661902) q[3];
sx q[3];
rz(-2.3016774) q[3];
sx q[3];
rz(-1.1168787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53330082) q[0];
sx q[0];
rz(-0.58654439) q[0];
sx q[0];
rz(0.71304148) q[0];
rz(-2.9167922) q[1];
sx q[1];
rz(-0.59200042) q[1];
sx q[1];
rz(1.0932659) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3521234) q[0];
sx q[0];
rz(-2.2575245) q[0];
sx q[0];
rz(2.5144469) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9399038) q[2];
sx q[2];
rz(-1.301577) q[2];
sx q[2];
rz(1.9112916) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.097979989) q[1];
sx q[1];
rz(-2.671125) q[1];
sx q[1];
rz(-0.2474252) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3087464) q[3];
sx q[3];
rz(-2.7152677) q[3];
sx q[3];
rz(0.48708068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1405868) q[2];
sx q[2];
rz(-2.4457377) q[2];
sx q[2];
rz(0.33983964) q[2];
rz(-3.0991683) q[3];
sx q[3];
rz(-1.8570615) q[3];
sx q[3];
rz(-0.62780082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.18702678) q[0];
sx q[0];
rz(-1.5196336) q[0];
sx q[0];
rz(-1.2484311) q[0];
rz(2.3595702) q[1];
sx q[1];
rz(-0.93688688) q[1];
sx q[1];
rz(-0.72601906) q[1];
rz(-0.72612957) q[2];
sx q[2];
rz(-0.48116044) q[2];
sx q[2];
rz(3.0143123) q[2];
rz(0.25039582) q[3];
sx q[3];
rz(-1.9848817) q[3];
sx q[3];
rz(3.1105697) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
