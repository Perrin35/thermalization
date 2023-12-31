OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.36800185) q[0];
sx q[0];
rz(-0.79080963) q[0];
sx q[0];
rz(0.33413449) q[0];
rz(2.6842527) q[1];
sx q[1];
rz(-2.1973124) q[1];
sx q[1];
rz(1.9231208) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5028635) q[0];
sx q[0];
rz(-1.6067061) q[0];
sx q[0];
rz(-0.7508276) q[0];
x q[1];
rz(-0.46484868) q[2];
sx q[2];
rz(-1.6072304) q[2];
sx q[2];
rz(-2.6236617) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4301181) q[1];
sx q[1];
rz(-2.4415486) q[1];
sx q[1];
rz(-0.75573604) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86745947) q[3];
sx q[3];
rz(-0.28775035) q[3];
sx q[3];
rz(-2.5301463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.618764) q[2];
sx q[2];
rz(-0.4814119) q[2];
sx q[2];
rz(-2.5640326) q[2];
rz(-1.1497568) q[3];
sx q[3];
rz(-1.7532319) q[3];
sx q[3];
rz(-2.4770588) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57099045) q[0];
sx q[0];
rz(-1.6002858) q[0];
sx q[0];
rz(-1.5748613) q[0];
x q[1];
rz(-2.8005373) q[2];
sx q[2];
rz(-2.555072) q[2];
sx q[2];
rz(-1.3618493) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.034033876) q[1];
sx q[1];
rz(-0.65980136) q[1];
sx q[1];
rz(-0.69114139) q[1];
x q[2];
rz(3.0121147) q[3];
sx q[3];
rz(-2.5054512) q[3];
sx q[3];
rz(0.57487088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.42276057) q[2];
sx q[2];
rz(-1.7451124) q[2];
sx q[2];
rz(-0.31769162) q[2];
rz(-0.20673949) q[3];
sx q[3];
rz(-0.59967774) q[3];
sx q[3];
rz(2.3247705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7725672) q[0];
sx q[0];
rz(-1.7018397) q[0];
sx q[0];
rz(1.7279708) q[0];
rz(-2.6638022) q[1];
sx q[1];
rz(-1.3505892) q[1];
sx q[1];
rz(-2.7405222) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-2.4972649) q[2];
sx q[2];
rz(-1.7228848) q[2];
sx q[2];
rz(-2.6459141) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9123147) q[1];
sx q[1];
rz(-2.4114128) q[1];
sx q[1];
rz(0.0095403949) q[1];
x q[2];
rz(2.9455455) q[3];
sx q[3];
rz(-1.184433) q[3];
sx q[3];
rz(2.8783609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.59427375) q[2];
sx q[2];
rz(-1.5218364) q[2];
sx q[2];
rz(-2.5857914) q[2];
rz(-0.9764955) q[3];
sx q[3];
rz(-0.54978168) q[3];
sx q[3];
rz(-2.3613789) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34898409) q[0];
sx q[0];
rz(-1.5190834) q[0];
sx q[0];
rz(-1.4439616) q[0];
rz(1.5199039) q[1];
sx q[1];
rz(-0.65602055) q[1];
sx q[1];
rz(-2.8881883) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2052106) q[0];
sx q[0];
rz(-1.0316327) q[0];
sx q[0];
rz(-0.11986952) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6035945) q[2];
sx q[2];
rz(-0.48600733) q[2];
sx q[2];
rz(-0.82644586) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6654012) q[1];
sx q[1];
rz(-2.5051077) q[1];
sx q[1];
rz(-2.979216) q[1];
rz(-2.4921791) q[3];
sx q[3];
rz(-1.489349) q[3];
sx q[3];
rz(-0.34929517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0306586) q[2];
sx q[2];
rz(-1.7548283) q[2];
sx q[2];
rz(-2.3542662) q[2];
rz(2.2287255) q[3];
sx q[3];
rz(-2.392231) q[3];
sx q[3];
rz(1.0095989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.71516365) q[0];
sx q[0];
rz(-0.63868317) q[0];
sx q[0];
rz(-3.0786247) q[0];
rz(-0.12403034) q[1];
sx q[1];
rz(-0.80563671) q[1];
sx q[1];
rz(-2.6834992) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3315017) q[0];
sx q[0];
rz(-2.0134263) q[0];
sx q[0];
rz(-1.5682194) q[0];
rz(-2.6890254) q[2];
sx q[2];
rz(-2.1308225) q[2];
sx q[2];
rz(-2.3134311) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.749436) q[1];
sx q[1];
rz(-1.4538308) q[1];
sx q[1];
rz(0.82526308) q[1];
rz(-pi) q[2];
rz(-3.0466988) q[3];
sx q[3];
rz(-1.0720383) q[3];
sx q[3];
rz(-0.14548485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1725585) q[2];
sx q[2];
rz(-2.2183552) q[2];
sx q[2];
rz(-2.9120973) q[2];
rz(-3.138792) q[3];
sx q[3];
rz(-0.87001785) q[3];
sx q[3];
rz(1.8026479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5795508) q[0];
sx q[0];
rz(-2.8631449) q[0];
sx q[0];
rz(0.91947412) q[0];
rz(3.0793076) q[1];
sx q[1];
rz(-1.0039763) q[1];
sx q[1];
rz(1.2671635) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4290659) q[0];
sx q[0];
rz(-2.292284) q[0];
sx q[0];
rz(2.6119786) q[0];
x q[1];
rz(-2.3976372) q[2];
sx q[2];
rz(-1.6774872) q[2];
sx q[2];
rz(0.26847408) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2515182) q[1];
sx q[1];
rz(-1.2444082) q[1];
sx q[1];
rz(-2.9337487) q[1];
rz(-pi) q[2];
rz(-1.1250161) q[3];
sx q[3];
rz(-0.52281724) q[3];
sx q[3];
rz(-1.0558053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1798114) q[2];
sx q[2];
rz(-2.2634025) q[2];
sx q[2];
rz(-2.5578257) q[2];
rz(-0.70872712) q[3];
sx q[3];
rz(-1.830359) q[3];
sx q[3];
rz(0.023199737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0631183) q[0];
sx q[0];
rz(-2.6874976) q[0];
sx q[0];
rz(1.0725347) q[0];
rz(-2.5947) q[1];
sx q[1];
rz(-1.2439589) q[1];
sx q[1];
rz(-2.0297208) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8060914) q[0];
sx q[0];
rz(-1.4772692) q[0];
sx q[0];
rz(2.0557687) q[0];
rz(-pi) q[1];
rz(-2.578031) q[2];
sx q[2];
rz(-2.5346018) q[2];
sx q[2];
rz(2.3129472) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0791196) q[1];
sx q[1];
rz(-1.8603431) q[1];
sx q[1];
rz(-1.6505961) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7518919) q[3];
sx q[3];
rz(-1.1003564) q[3];
sx q[3];
rz(-2.5999992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.303858) q[2];
sx q[2];
rz(-1.6656817) q[2];
sx q[2];
rz(-1.973935) q[2];
rz(-1.6052823) q[3];
sx q[3];
rz(-1.6746018) q[3];
sx q[3];
rz(-1.7355828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4090356) q[0];
sx q[0];
rz(-1.5696101) q[0];
sx q[0];
rz(-0.73079601) q[0];
rz(-0.90019512) q[1];
sx q[1];
rz(-0.80454818) q[1];
sx q[1];
rz(2.3866167) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1599931) q[0];
sx q[0];
rz(-0.85708517) q[0];
sx q[0];
rz(0.1937565) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.99004284) q[2];
sx q[2];
rz(-2.2120737) q[2];
sx q[2];
rz(-2.5255447) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0874487) q[1];
sx q[1];
rz(-2.1463697) q[1];
sx q[1];
rz(-1.7051484) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.98324361) q[3];
sx q[3];
rz(-1.1508815) q[3];
sx q[3];
rz(-1.1004694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.76453152) q[2];
sx q[2];
rz(-1.3724644) q[2];
sx q[2];
rz(-2.5047452) q[2];
rz(-0.26646715) q[3];
sx q[3];
rz(-2.0839432) q[3];
sx q[3];
rz(1.586097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72717845) q[0];
sx q[0];
rz(-1.1295015) q[0];
sx q[0];
rz(1.138858) q[0];
rz(-2.3873734) q[1];
sx q[1];
rz(-2.8051839) q[1];
sx q[1];
rz(3.1220904) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50197983) q[0];
sx q[0];
rz(-1.4907955) q[0];
sx q[0];
rz(1.3638858) q[0];
x q[1];
rz(-2.2150061) q[2];
sx q[2];
rz(-2.8928061) q[2];
sx q[2];
rz(2.8015346) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6448977) q[1];
sx q[1];
rz(-0.95458191) q[1];
sx q[1];
rz(0.22209054) q[1];
x q[2];
rz(-0.1653413) q[3];
sx q[3];
rz(-1.8792361) q[3];
sx q[3];
rz(2.7066018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1372244) q[2];
sx q[2];
rz(-1.7251816) q[2];
sx q[2];
rz(0.84890378) q[2];
rz(2.7539339) q[3];
sx q[3];
rz(-1.1281697) q[3];
sx q[3];
rz(-1.5415812) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7983109) q[0];
sx q[0];
rz(-0.16769519) q[0];
sx q[0];
rz(-0.48450255) q[0];
rz(-1.3867406) q[1];
sx q[1];
rz(-1.4258899) q[1];
sx q[1];
rz(-1.9932995) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1295706) q[0];
sx q[0];
rz(-0.14229933) q[0];
sx q[0];
rz(-2.9905031) q[0];
rz(1.7933153) q[2];
sx q[2];
rz(-2.1902124) q[2];
sx q[2];
rz(0.36048181) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1710098) q[1];
sx q[1];
rz(-2.2323771) q[1];
sx q[1];
rz(-0.2172825) q[1];
rz(-pi) q[2];
x q[2];
rz(0.34928068) q[3];
sx q[3];
rz(-1.619907) q[3];
sx q[3];
rz(1.9625361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.91223532) q[2];
sx q[2];
rz(-1.2939913) q[2];
sx q[2];
rz(-2.7764376) q[2];
rz(-3.0129516) q[3];
sx q[3];
rz(-1.235685) q[3];
sx q[3];
rz(0.45583367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(2.4776496) q[2];
sx q[2];
rz(-1.2524111) q[2];
sx q[2];
rz(-1.8128957) q[2];
rz(0.036988463) q[3];
sx q[3];
rz(-1.0109517) q[3];
sx q[3];
rz(3.0299822) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
