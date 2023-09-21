OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7735908) q[0];
sx q[0];
rz(-2.350783) q[0];
sx q[0];
rz(2.8074582) q[0];
rz(2.6842527) q[1];
sx q[1];
rz(-2.1973124) q[1];
sx q[1];
rz(1.9231208) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5028635) q[0];
sx q[0];
rz(-1.6067061) q[0];
sx q[0];
rz(-2.3907651) q[0];
x q[1];
rz(-0.081131447) q[2];
sx q[2];
rz(-0.46617026) q[2];
sx q[2];
rz(0.98035882) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6005046) q[1];
sx q[1];
rz(-2.0588015) q[1];
sx q[1];
rz(-2.0946676) q[1];
rz(1.7928042) q[3];
sx q[3];
rz(-1.3862002) q[3];
sx q[3];
rz(0.27664646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.618764) q[2];
sx q[2];
rz(-2.6601807) q[2];
sx q[2];
rz(-2.5640326) q[2];
rz(-1.9918359) q[3];
sx q[3];
rz(-1.3883608) q[3];
sx q[3];
rz(-2.4770588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98786551) q[0];
sx q[0];
rz(-0.58652121) q[0];
sx q[0];
rz(2.7541449) q[0];
rz(2.2024343) q[1];
sx q[1];
rz(-2.1444131) q[1];
sx q[1];
rz(-1.739025) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70799202) q[0];
sx q[0];
rz(-0.029768243) q[0];
sx q[0];
rz(3.004651) q[0];
rz(-1.3520794) q[2];
sx q[2];
rz(-1.0220851) q[2];
sx q[2];
rz(-0.95900853) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1837511) q[1];
sx q[1];
rz(-1.9722003) q[1];
sx q[1];
rz(-0.5387696) q[1];
rz(-pi) q[2];
rz(-2.5094633) q[3];
sx q[3];
rz(-1.4940133) q[3];
sx q[3];
rz(-1.1002822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.42276057) q[2];
sx q[2];
rz(-1.7451124) q[2];
sx q[2];
rz(-2.823901) q[2];
rz(2.9348532) q[3];
sx q[3];
rz(-2.5419149) q[3];
sx q[3];
rz(0.81682214) q[3];
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
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3690255) q[0];
sx q[0];
rz(-1.7018397) q[0];
sx q[0];
rz(1.4136219) q[0];
rz(0.47779045) q[1];
sx q[1];
rz(-1.7910035) q[1];
sx q[1];
rz(2.7405222) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.385752) q[0];
sx q[0];
rz(-2.7999561) q[0];
sx q[0];
rz(-1.8833478) q[0];
rz(-pi) q[1];
rz(-0.6443278) q[2];
sx q[2];
rz(-1.7228848) q[2];
sx q[2];
rz(2.6459141) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.33441015) q[1];
sx q[1];
rz(-1.5771598) q[1];
sx q[1];
rz(2.4114354) q[1];
rz(-1.1242261) q[3];
sx q[3];
rz(-0.43101573) q[3];
sx q[3];
rz(-2.9197846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5473189) q[2];
sx q[2];
rz(-1.6197562) q[2];
sx q[2];
rz(-0.55580124) q[2];
rz(2.1650971) q[3];
sx q[3];
rz(-0.54978168) q[3];
sx q[3];
rz(0.78021375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34898409) q[0];
sx q[0];
rz(-1.6225092) q[0];
sx q[0];
rz(1.697631) q[0];
rz(-1.6216888) q[1];
sx q[1];
rz(-0.65602055) q[1];
sx q[1];
rz(0.25340432) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69617535) q[0];
sx q[0];
rz(-1.6735958) q[0];
sx q[0];
rz(-2.1131383) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0565815) q[2];
sx q[2];
rz(-1.5554785) q[2];
sx q[2];
rz(-2.3682396) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6654012) q[1];
sx q[1];
rz(-0.63648495) q[1];
sx q[1];
rz(2.979216) q[1];
rz(-pi) q[2];
rz(-2.4921791) q[3];
sx q[3];
rz(-1.6522437) q[3];
sx q[3];
rz(-2.7922975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0306586) q[2];
sx q[2];
rz(-1.7548283) q[2];
sx q[2];
rz(-0.78732642) q[2];
rz(-2.2287255) q[3];
sx q[3];
rz(-0.74936167) q[3];
sx q[3];
rz(-2.1319938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.426429) q[0];
sx q[0];
rz(-0.63868317) q[0];
sx q[0];
rz(-0.062967904) q[0];
rz(-0.12403034) q[1];
sx q[1];
rz(-0.80563671) q[1];
sx q[1];
rz(-2.6834992) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23819085) q[0];
sx q[0];
rz(-1.5731249) q[0];
sx q[0];
rz(-2.6989614) q[0];
rz(0.96197084) q[2];
sx q[2];
rz(-1.1912727) q[2];
sx q[2];
rz(-2.6517207) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.071306989) q[1];
sx q[1];
rz(-2.3100393) q[1];
sx q[1];
rz(-2.9830095) q[1];
x q[2];
rz(1.7430274) q[3];
sx q[3];
rz(-0.50695626) q[3];
sx q[3];
rz(-2.7996922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1725585) q[2];
sx q[2];
rz(-0.92323747) q[2];
sx q[2];
rz(0.22949533) q[2];
rz(-0.0028006639) q[3];
sx q[3];
rz(-2.2715748) q[3];
sx q[3];
rz(-1.3389448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5795508) q[0];
sx q[0];
rz(-0.27844772) q[0];
sx q[0];
rz(-0.91947412) q[0];
rz(0.062285034) q[1];
sx q[1];
rz(-1.0039763) q[1];
sx q[1];
rz(-1.2671635) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7125268) q[0];
sx q[0];
rz(-0.8493087) q[0];
sx q[0];
rz(-2.6119786) q[0];
rz(0.74395545) q[2];
sx q[2];
rz(-1.4641054) q[2];
sx q[2];
rz(-0.26847408) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8333203) q[1];
sx q[1];
rz(-2.7566524) q[1];
sx q[1];
rz(-2.1182548) q[1];
x q[2];
rz(-0.24354981) q[3];
sx q[3];
rz(-1.1034414) q[3];
sx q[3];
rz(1.5817643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.9617812) q[2];
sx q[2];
rz(-0.87819019) q[2];
sx q[2];
rz(2.5578257) q[2];
rz(2.4328655) q[3];
sx q[3];
rz(-1.3112336) q[3];
sx q[3];
rz(3.1183929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07847438) q[0];
sx q[0];
rz(-2.6874976) q[0];
sx q[0];
rz(-2.069058) q[0];
rz(-2.5947) q[1];
sx q[1];
rz(-1.8976338) q[1];
sx q[1];
rz(-1.1118719) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5892964) q[0];
sx q[0];
rz(-2.6483905) q[0];
sx q[0];
rz(-1.3722377) q[0];
x q[1];
rz(1.215559) q[2];
sx q[2];
rz(-2.07395) q[2];
sx q[2];
rz(2.9687198) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3355616) q[1];
sx q[1];
rz(-0.30004382) q[1];
sx q[1];
rz(2.8801444) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92915793) q[3];
sx q[3];
rz(-2.540179) q[3];
sx q[3];
rz(1.8638368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.83773461) q[2];
sx q[2];
rz(-1.6656817) q[2];
sx q[2];
rz(-1.1676577) q[2];
rz(1.6052823) q[3];
sx q[3];
rz(-1.6746018) q[3];
sx q[3];
rz(-1.4060098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4090356) q[0];
sx q[0];
rz(-1.5696101) q[0];
sx q[0];
rz(2.4107966) q[0];
rz(-2.2413975) q[1];
sx q[1];
rz(-0.80454818) q[1];
sx q[1];
rz(0.75497595) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98159957) q[0];
sx q[0];
rz(-2.2845075) q[0];
sx q[0];
rz(0.1937565) q[0];
rz(-pi) q[1];
rz(2.412699) q[2];
sx q[2];
rz(-2.0260099) q[2];
sx q[2];
rz(-0.58065562) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6983812) q[1];
sx q[1];
rz(-1.6834007) q[1];
sx q[1];
rz(-2.5618782) q[1];
rz(-pi) q[2];
rz(2.158349) q[3];
sx q[3];
rz(-1.9907111) q[3];
sx q[3];
rz(1.1004694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3770611) q[2];
sx q[2];
rz(-1.3724644) q[2];
sx q[2];
rz(-2.5047452) q[2];
rz(-2.8751255) q[3];
sx q[3];
rz(-1.0576495) q[3];
sx q[3];
rz(1.586097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4144142) q[0];
sx q[0];
rz(-1.1295015) q[0];
sx q[0];
rz(-1.138858) q[0];
rz(-2.3873734) q[1];
sx q[1];
rz(-0.33640877) q[1];
sx q[1];
rz(0.019502217) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50197983) q[0];
sx q[0];
rz(-1.6507971) q[0];
sx q[0];
rz(-1.3638858) q[0];
rz(-2.2150061) q[2];
sx q[2];
rz(-2.8928061) q[2];
sx q[2];
rz(-0.34005806) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.94432482) q[1];
sx q[1];
rz(-1.3900583) q[1];
sx q[1];
rz(-2.1987869) q[1];
rz(-pi) q[2];
rz(2.9762514) q[3];
sx q[3];
rz(-1.8792361) q[3];
sx q[3];
rz(-0.43499085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1372244) q[2];
sx q[2];
rz(-1.4164111) q[2];
sx q[2];
rz(-2.2926889) q[2];
rz(2.7539339) q[3];
sx q[3];
rz(-2.013423) q[3];
sx q[3];
rz(1.5415812) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3432817) q[0];
sx q[0];
rz(-0.16769519) q[0];
sx q[0];
rz(2.6570901) q[0];
rz(-1.7548521) q[1];
sx q[1];
rz(-1.7157028) q[1];
sx q[1];
rz(-1.9932995) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5507817) q[0];
sx q[0];
rz(-1.5494487) q[0];
sx q[0];
rz(0.14069964) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63126385) q[2];
sx q[2];
rz(-1.3901276) q[2];
sx q[2];
rz(-1.0797015) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9705829) q[1];
sx q[1];
rz(-0.90921558) q[1];
sx q[1];
rz(-0.2172825) q[1];
x q[2];
rz(0.34928068) q[3];
sx q[3];
rz(-1.5216856) q[3];
sx q[3];
rz(1.1790566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.91223532) q[2];
sx q[2];
rz(-1.8476013) q[2];
sx q[2];
rz(0.36515507) q[2];
rz(-0.12864104) q[3];
sx q[3];
rz(-1.9059076) q[3];
sx q[3];
rz(0.45583367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1098332) q[0];
sx q[0];
rz(-2.3008627) q[0];
sx q[0];
rz(-1.536137) q[0];
rz(-0.96314349) q[1];
sx q[1];
rz(-1.2711202) q[1];
sx q[1];
rz(-1.0585379) q[1];
rz(-1.9671494) q[2];
sx q[2];
rz(-2.195993) q[2];
sx q[2];
rz(3.1396951) q[2];
rz(1.6297324) q[3];
sx q[3];
rz(-2.5806576) q[3];
sx q[3];
rz(2.9604119) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
