OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.84500399) q[0];
sx q[0];
rz(-0.72405976) q[0];
sx q[0];
rz(-1.6568503) q[0];
rz(4.3447189) q[1];
sx q[1];
rz(0.523518) q[1];
sx q[1];
rz(8.5365774) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0892544) q[0];
sx q[0];
rz(-1.6153781) q[0];
sx q[0];
rz(0.17104878) q[0];
rz(-pi) q[1];
rz(1.0339123) q[2];
sx q[2];
rz(-1.4281338) q[2];
sx q[2];
rz(-0.85096525) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8623212) q[1];
sx q[1];
rz(-0.16695484) q[1];
sx q[1];
rz(-3.1382986) q[1];
rz(0.69859759) q[3];
sx q[3];
rz(-1.5961831) q[3];
sx q[3];
rz(-0.7198402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.28329864) q[2];
sx q[2];
rz(-2.7259939) q[2];
sx q[2];
rz(1.1179914) q[2];
rz(-2.9962712) q[3];
sx q[3];
rz(-1.5829007) q[3];
sx q[3];
rz(-0.045923559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5685101) q[0];
sx q[0];
rz(-1.2643603) q[0];
sx q[0];
rz(-2.4867687) q[0];
rz(1.9251992) q[1];
sx q[1];
rz(-1.1647859) q[1];
sx q[1];
rz(0.12589802) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1032216) q[0];
sx q[0];
rz(-1.201259) q[0];
sx q[0];
rz(2.080337) q[0];
rz(-pi) q[1];
x q[1];
rz(0.89181487) q[2];
sx q[2];
rz(-1.0620772) q[2];
sx q[2];
rz(1.5430792) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.8996846) q[1];
sx q[1];
rz(-0.13026127) q[1];
sx q[1];
rz(1.3379407) q[1];
x q[2];
rz(1.5586073) q[3];
sx q[3];
rz(-0.72353957) q[3];
sx q[3];
rz(-0.42750588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0931603) q[2];
sx q[2];
rz(-1.9530714) q[2];
sx q[2];
rz(-2.753567) q[2];
rz(1.7175425) q[3];
sx q[3];
rz(-2.5035796) q[3];
sx q[3];
rz(0.58732906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5557142) q[0];
sx q[0];
rz(-2.3590187) q[0];
sx q[0];
rz(0.079332381) q[0];
rz(-0.084005984) q[1];
sx q[1];
rz(-2.3386798) q[1];
sx q[1];
rz(1.1598587) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9281611) q[0];
sx q[0];
rz(-2.0800989) q[0];
sx q[0];
rz(3.0470949) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2247505) q[2];
sx q[2];
rz(-2.4008304) q[2];
sx q[2];
rz(-3.0286718) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0821655) q[1];
sx q[1];
rz(-1.1929999) q[1];
sx q[1];
rz(-0.49281812) q[1];
rz(-pi) q[2];
rz(-0.056315259) q[3];
sx q[3];
rz(-1.0259797) q[3];
sx q[3];
rz(-2.44851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7362061) q[2];
sx q[2];
rz(-1.3233041) q[2];
sx q[2];
rz(-3.0333056) q[2];
rz(-2.4978499) q[3];
sx q[3];
rz(-2.0635922) q[3];
sx q[3];
rz(0.61557499) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.165034) q[0];
sx q[0];
rz(-1.3950011) q[0];
sx q[0];
rz(-1.4820341) q[0];
rz(2.4404793) q[1];
sx q[1];
rz(-1.2524403) q[1];
sx q[1];
rz(0.28465095) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68344224) q[0];
sx q[0];
rz(-1.0939286) q[0];
sx q[0];
rz(-3.0577858) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8585763) q[2];
sx q[2];
rz(-1.3735526) q[2];
sx q[2];
rz(1.6523199) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9173911) q[1];
sx q[1];
rz(-2.9926139) q[1];
sx q[1];
rz(-1.0148744) q[1];
x q[2];
rz(-2.1553667) q[3];
sx q[3];
rz(-0.75891906) q[3];
sx q[3];
rz(0.81134568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9099137) q[2];
sx q[2];
rz(-2.173013) q[2];
sx q[2];
rz(-2.5740734) q[2];
rz(0.41401687) q[3];
sx q[3];
rz(-1.1375789) q[3];
sx q[3];
rz(1.1119941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.652997) q[0];
sx q[0];
rz(-0.96019205) q[0];
sx q[0];
rz(1.8150785) q[0];
rz(-1.9891706) q[1];
sx q[1];
rz(-1.7633341) q[1];
sx q[1];
rz(-2.2036536) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65813488) q[0];
sx q[0];
rz(-1.5489274) q[0];
sx q[0];
rz(-1.4786426) q[0];
rz(-pi) q[1];
rz(2.6229834) q[2];
sx q[2];
rz(-1.0728288) q[2];
sx q[2];
rz(1.2155611) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2873958) q[1];
sx q[1];
rz(-1.5970267) q[1];
sx q[1];
rz(-1.5363974) q[1];
rz(-pi) q[2];
rz(1.7004847) q[3];
sx q[3];
rz(-1.9228336) q[3];
sx q[3];
rz(0.45613939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9662629) q[2];
sx q[2];
rz(-2.9523409) q[2];
sx q[2];
rz(-1.0413292) q[2];
rz(1.8390309) q[3];
sx q[3];
rz(-2.0077191) q[3];
sx q[3];
rz(-1.3180102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.702521) q[0];
sx q[0];
rz(-1.1477926) q[0];
sx q[0];
rz(0.55737108) q[0];
rz(2.5769261) q[1];
sx q[1];
rz(-0.70960418) q[1];
sx q[1];
rz(-2.5851137) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51150409) q[0];
sx q[0];
rz(-0.43404365) q[0];
sx q[0];
rz(-0.50891288) q[0];
rz(0.17369341) q[2];
sx q[2];
rz(-2.2787893) q[2];
sx q[2];
rz(-0.049335418) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6528875) q[1];
sx q[1];
rz(-0.31213752) q[1];
sx q[1];
rz(-2.1991792) q[1];
rz(1.6218833) q[3];
sx q[3];
rz(-0.79569492) q[3];
sx q[3];
rz(1.5339799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.53211987) q[2];
sx q[2];
rz(-0.76081053) q[2];
sx q[2];
rz(-1.2247941) q[2];
rz(1.5103643) q[3];
sx q[3];
rz(-1.8211726) q[3];
sx q[3];
rz(-1.640655) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0625793) q[0];
sx q[0];
rz(-1.2591079) q[0];
sx q[0];
rz(-0.042908948) q[0];
rz(0.91730109) q[1];
sx q[1];
rz(-0.62364548) q[1];
sx q[1];
rz(-2.6409805) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0450889) q[0];
sx q[0];
rz(-1.6052264) q[0];
sx q[0];
rz(1.7486497) q[0];
rz(-pi) q[1];
x q[1];
rz(2.826564) q[2];
sx q[2];
rz(-1.050596) q[2];
sx q[2];
rz(-2.5556285) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.73270479) q[1];
sx q[1];
rz(-1.171247) q[1];
sx q[1];
rz(-1.4399745) q[1];
x q[2];
rz(-2.72243) q[3];
sx q[3];
rz(-2.3401642) q[3];
sx q[3];
rz(-2.6475346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5381955) q[2];
sx q[2];
rz(-1.0791225) q[2];
sx q[2];
rz(2.4659757) q[2];
rz(2.7006941) q[3];
sx q[3];
rz(-1.4392122) q[3];
sx q[3];
rz(-0.52136695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.065141) q[0];
sx q[0];
rz(-2.8171709) q[0];
sx q[0];
rz(-2.0741529) q[0];
rz(-2.7087129) q[1];
sx q[1];
rz(-1.6378816) q[1];
sx q[1];
rz(1.4656461) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7536403) q[0];
sx q[0];
rz(-0.48250972) q[0];
sx q[0];
rz(-1.8637191) q[0];
rz(-pi) q[1];
x q[1];
rz(0.52532105) q[2];
sx q[2];
rz(-2.9267731) q[2];
sx q[2];
rz(-2.7810682) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9626179) q[1];
sx q[1];
rz(-1.0777506) q[1];
sx q[1];
rz(-0.38260539) q[1];
x q[2];
rz(0.96887178) q[3];
sx q[3];
rz(-1.0191917) q[3];
sx q[3];
rz(-1.5495891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.33891588) q[2];
sx q[2];
rz(-2.2854476) q[2];
sx q[2];
rz(-1.6652997) q[2];
rz(-0.87351292) q[3];
sx q[3];
rz(-2.2647808) q[3];
sx q[3];
rz(-0.4666369) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8840238) q[0];
sx q[0];
rz(-0.57618657) q[0];
sx q[0];
rz(-2.4043758) q[0];
rz(3.1228512) q[1];
sx q[1];
rz(-2.811921) q[1];
sx q[1];
rz(-2.2163056) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5718482) q[0];
sx q[0];
rz(-1.0751372) q[0];
sx q[0];
rz(0.75832383) q[0];
rz(-pi) q[1];
rz(2.972702) q[2];
sx q[2];
rz(-2.4714111) q[2];
sx q[2];
rz(0.18667135) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9097594) q[1];
sx q[1];
rz(-1.701046) q[1];
sx q[1];
rz(-0.18472437) q[1];
x q[2];
rz(0.41947375) q[3];
sx q[3];
rz(-0.89722108) q[3];
sx q[3];
rz(0.89299612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.76901889) q[2];
sx q[2];
rz(-1.6071268) q[2];
sx q[2];
rz(-1.0037237) q[2];
rz(-3.051493) q[3];
sx q[3];
rz(-3.1153479) q[3];
sx q[3];
rz(2.0285006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.039624778) q[0];
sx q[0];
rz(-1.5682546) q[0];
sx q[0];
rz(1.8819303) q[0];
rz(0.26578495) q[1];
sx q[1];
rz(-0.62756413) q[1];
sx q[1];
rz(-2.3840747) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.789061) q[0];
sx q[0];
rz(-1.7081982) q[0];
sx q[0];
rz(-1.765541) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1610356) q[2];
sx q[2];
rz(-1.7014116) q[2];
sx q[2];
rz(-0.57934258) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.65749189) q[1];
sx q[1];
rz(-1.7925646) q[1];
sx q[1];
rz(-2.0055254) q[1];
x q[2];
rz(1.5991454) q[3];
sx q[3];
rz(-2.2806014) q[3];
sx q[3];
rz(2.860644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.99047986) q[2];
sx q[2];
rz(-1.8258784) q[2];
sx q[2];
rz(-0.79375664) q[2];
rz(-2.5027067) q[3];
sx q[3];
rz(-0.34404889) q[3];
sx q[3];
rz(-1.0206153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4929852) q[0];
sx q[0];
rz(-2.6633371) q[0];
sx q[0];
rz(2.2289842) q[0];
rz(1.6408625) q[1];
sx q[1];
rz(-2.224557) q[1];
sx q[1];
rz(1.7932737) q[1];
rz(1.2407672) q[2];
sx q[2];
rz(-1.1541661) q[2];
sx q[2];
rz(3.0164568) q[2];
rz(0.32140857) q[3];
sx q[3];
rz(-1.0151498) q[3];
sx q[3];
rz(-2.2890454) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
