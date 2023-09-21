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
rz(1.4847423) q[0];
rz(1.2031263) q[1];
sx q[1];
rz(-0.523518) q[1];
sx q[1];
rz(2.2533921) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0523383) q[0];
sx q[0];
rz(-1.6153781) q[0];
sx q[0];
rz(0.17104878) q[0];
rz(-pi) q[1];
rz(1.2970096) q[2];
sx q[2];
rz(-2.5878776) q[2];
sx q[2];
rz(-0.48534976) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.29477316) q[1];
sx q[1];
rz(-1.5702489) q[1];
sx q[1];
rz(2.9746387) q[1];
rz(-pi) q[2];
x q[2];
rz(0.039460823) q[3];
sx q[3];
rz(-0.69898116) q[3];
sx q[3];
rz(2.3208502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.28329864) q[2];
sx q[2];
rz(-0.41559872) q[2];
sx q[2];
rz(1.1179914) q[2];
rz(-0.14532146) q[3];
sx q[3];
rz(-1.5829007) q[3];
sx q[3];
rz(-3.0956691) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5685101) q[0];
sx q[0];
rz(-1.8772323) q[0];
sx q[0];
rz(0.65482393) q[0];
rz(1.9251992) q[1];
sx q[1];
rz(-1.9768068) q[1];
sx q[1];
rz(-0.12589802) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1830131) q[0];
sx q[0];
rz(-0.61972451) q[0];
sx q[0];
rz(-2.241914) q[0];
rz(-pi) q[1];
rz(2.5198031) q[2];
sx q[2];
rz(-2.1513373) q[2];
sx q[2];
rz(0.34678005) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0071348) q[1];
sx q[1];
rz(-1.4440698) q[1];
sx q[1];
rz(0.030220672) q[1];
x q[2];
rz(-2.294299) q[3];
sx q[3];
rz(-1.5627268) q[3];
sx q[3];
rz(-1.9891667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.048432365) q[2];
sx q[2];
rz(-1.1885213) q[2];
sx q[2];
rz(-0.38802567) q[2];
rz(1.7175425) q[3];
sx q[3];
rz(-0.63801304) q[3];
sx q[3];
rz(2.5542636) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5858784) q[0];
sx q[0];
rz(-2.3590187) q[0];
sx q[0];
rz(-0.079332381) q[0];
rz(-3.0575867) q[1];
sx q[1];
rz(-0.80291286) q[1];
sx q[1];
rz(1.1598587) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9281611) q[0];
sx q[0];
rz(-2.0800989) q[0];
sx q[0];
rz(-3.0470949) q[0];
rz(-pi) q[1];
rz(0.91684219) q[2];
sx q[2];
rz(-0.7407623) q[2];
sx q[2];
rz(3.0286718) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0821655) q[1];
sx q[1];
rz(-1.1929999) q[1];
sx q[1];
rz(2.6487745) q[1];
rz(-pi) q[2];
x q[2];
rz(0.056315259) q[3];
sx q[3];
rz(-1.0259797) q[3];
sx q[3];
rz(-0.69308263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.40538654) q[2];
sx q[2];
rz(-1.3233041) q[2];
sx q[2];
rz(-0.1082871) q[2];
rz(-0.64374271) q[3];
sx q[3];
rz(-2.0635922) q[3];
sx q[3];
rz(2.5260177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.4404793) q[1];
sx q[1];
rz(-1.2524403) q[1];
sx q[1];
rz(2.8569417) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92589256) q[0];
sx q[0];
rz(-1.6452351) q[0];
sx q[0];
rz(-2.0490993) q[0];
rz(2.1842723) q[2];
sx q[2];
rz(-2.7942604) q[2];
sx q[2];
rz(-0.50328244) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.33672562) q[1];
sx q[1];
rz(-1.4443828) q[1];
sx q[1];
rz(-3.0625507) q[1];
rz(2.2399726) q[3];
sx q[3];
rz(-1.1812783) q[3];
sx q[3];
rz(0.31182409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9099137) q[2];
sx q[2];
rz(-2.173013) q[2];
sx q[2];
rz(0.56751928) q[2];
rz(-0.41401687) q[3];
sx q[3];
rz(-2.0040138) q[3];
sx q[3];
rz(1.1119941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.652997) q[0];
sx q[0];
rz(-2.1814006) q[0];
sx q[0];
rz(1.3265142) q[0];
rz(-1.9891706) q[1];
sx q[1];
rz(-1.7633341) q[1];
sx q[1];
rz(-2.2036536) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4612761) q[0];
sx q[0];
rz(-0.094705908) q[0];
sx q[0];
rz(-1.8041496) q[0];
rz(1.0114952) q[2];
sx q[2];
rz(-1.1202295) q[2];
sx q[2];
rz(0.089103854) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.282498) q[1];
sx q[1];
rz(-1.6051834) q[1];
sx q[1];
rz(-3.1153468) q[1];
rz(-pi) q[2];
rz(1.4411079) q[3];
sx q[3];
rz(-1.9228336) q[3];
sx q[3];
rz(-0.45613939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1753297) q[2];
sx q[2];
rz(-0.18925174) q[2];
sx q[2];
rz(1.0413292) q[2];
rz(-1.8390309) q[3];
sx q[3];
rz(-2.0077191) q[3];
sx q[3];
rz(1.3180102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43907169) q[0];
sx q[0];
rz(-1.1477926) q[0];
sx q[0];
rz(-0.55737108) q[0];
rz(0.56466651) q[1];
sx q[1];
rz(-2.4319885) q[1];
sx q[1];
rz(0.55647892) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0787449) q[0];
sx q[0];
rz(-1.1947462) q[0];
sx q[0];
rz(-1.7929121) q[0];
rz(-pi) q[1];
rz(-1.7700023) q[2];
sx q[2];
rz(-0.7253941) q[2];
sx q[2];
rz(2.9273916) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6528875) q[1];
sx q[1];
rz(-2.8294551) q[1];
sx q[1];
rz(-2.1991792) q[1];
rz(3.0895124) q[3];
sx q[3];
rz(-2.3651603) q[3];
sx q[3];
rz(1.6805502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.53211987) q[2];
sx q[2];
rz(-0.76081053) q[2];
sx q[2];
rz(1.9167985) q[2];
rz(1.6312284) q[3];
sx q[3];
rz(-1.8211726) q[3];
sx q[3];
rz(1.640655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-2.0790134) q[0];
sx q[0];
rz(-1.8824848) q[0];
sx q[0];
rz(-0.042908948) q[0];
rz(2.2242916) q[1];
sx q[1];
rz(-0.62364548) q[1];
sx q[1];
rz(-0.50061217) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0450889) q[0];
sx q[0];
rz(-1.6052264) q[0];
sx q[0];
rz(1.3929429) q[0];
rz(-pi) q[1];
rz(0.31502864) q[2];
sx q[2];
rz(-2.0909967) q[2];
sx q[2];
rz(-2.5556285) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7350406) q[1];
sx q[1];
rz(-2.7222689) q[1];
sx q[1];
rz(0.29962824) q[1];
rz(-1.9686437) q[3];
sx q[3];
rz(-2.2865191) q[3];
sx q[3];
rz(-2.0778823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5381955) q[2];
sx q[2];
rz(-2.0624702) q[2];
sx q[2];
rz(2.4659757) q[2];
rz(-0.44089857) q[3];
sx q[3];
rz(-1.4392122) q[3];
sx q[3];
rz(2.6202257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0764517) q[0];
sx q[0];
rz(-2.8171709) q[0];
sx q[0];
rz(-2.0741529) q[0];
rz(2.7087129) q[1];
sx q[1];
rz(-1.6378816) q[1];
sx q[1];
rz(-1.4656461) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4439125) q[0];
sx q[0];
rz(-1.4364103) q[0];
sx q[0];
rz(2.0356376) q[0];
rz(-pi) q[1];
rz(2.9550214) q[2];
sx q[2];
rz(-1.6779043) q[2];
sx q[2];
rz(0.69498108) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9379942) q[1];
sx q[1];
rz(-1.9059056) q[1];
sx q[1];
rz(2.0957698) q[1];
rz(-pi) q[2];
rz(-2.500324) q[3];
sx q[3];
rz(-2.0740168) q[3];
sx q[3];
rz(2.7748231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8026768) q[2];
sx q[2];
rz(-0.85614506) q[2];
sx q[2];
rz(-1.6652997) q[2];
rz(-2.2680797) q[3];
sx q[3];
rz(-2.2647808) q[3];
sx q[3];
rz(0.4666369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8840238) q[0];
sx q[0];
rz(-0.57618657) q[0];
sx q[0];
rz(2.4043758) q[0];
rz(-3.1228512) q[1];
sx q[1];
rz(-0.32967162) q[1];
sx q[1];
rz(0.92528701) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4243471) q[0];
sx q[0];
rz(-2.2204917) q[0];
sx q[0];
rz(-2.2109277) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66321744) q[2];
sx q[2];
rz(-1.6753917) q[2];
sx q[2];
rz(1.6246206) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9097594) q[1];
sx q[1];
rz(-1.701046) q[1];
sx q[1];
rz(0.18472437) q[1];
x q[2];
rz(2.0426644) q[3];
sx q[3];
rz(-2.3657551) q[3];
sx q[3];
rz(-2.8692506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3725738) q[2];
sx q[2];
rz(-1.5344658) q[2];
sx q[2];
rz(-2.1378689) q[2];
rz(-3.051493) q[3];
sx q[3];
rz(-0.026244791) q[3];
sx q[3];
rz(-2.0285006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1019679) q[0];
sx q[0];
rz(-1.573338) q[0];
sx q[0];
rz(1.2596624) q[0];
rz(-0.26578495) q[1];
sx q[1];
rz(-0.62756413) q[1];
sx q[1];
rz(-0.75751799) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3525317) q[0];
sx q[0];
rz(-1.7081982) q[0];
sx q[0];
rz(1.3760516) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1610356) q[2];
sx q[2];
rz(-1.440181) q[2];
sx q[2];
rz(0.57934258) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6706558) q[1];
sx q[1];
rz(-0.48479143) q[1];
sx q[1];
rz(2.0623341) q[1];
x q[2];
rz(3.1086139) q[3];
sx q[3];
rz(-2.4313201) q[3];
sx q[3];
rz(2.9041293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.99047986) q[2];
sx q[2];
rz(-1.8258784) q[2];
sx q[2];
rz(0.79375664) q[2];
rz(-2.5027067) q[3];
sx q[3];
rz(-0.34404889) q[3];
sx q[3];
rz(-1.0206153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4929852) q[0];
sx q[0];
rz(-2.6633371) q[0];
sx q[0];
rz(2.2289842) q[0];
rz(1.5007301) q[1];
sx q[1];
rz(-0.91703569) q[1];
sx q[1];
rz(-1.348319) q[1];
rz(0.43754229) q[2];
sx q[2];
rz(-1.2699288) q[2];
sx q[2];
rz(-1.5581836) q[2];
rz(1.100148) q[3];
sx q[3];
rz(-0.63334076) q[3];
sx q[3];
rz(0.28950194) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];