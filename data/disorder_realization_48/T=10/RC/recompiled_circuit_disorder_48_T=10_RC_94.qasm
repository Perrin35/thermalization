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
rz(-2.8074582) q[0];
rz(-0.45733991) q[1];
sx q[1];
rz(5.338905) q[1];
sx q[1];
rz(10.64325) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1760362) q[0];
sx q[0];
rz(-2.3210225) q[0];
sx q[0];
rz(1.5216989) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46484868) q[2];
sx q[2];
rz(-1.6072304) q[2];
sx q[2];
rz(0.51793098) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4301181) q[1];
sx q[1];
rz(-0.700044) q[1];
sx q[1];
rz(2.3858566) q[1];
rz(-0.18913194) q[3];
sx q[3];
rz(-1.7889708) q[3];
sx q[3];
rz(1.3355586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.618764) q[2];
sx q[2];
rz(-2.6601807) q[2];
sx q[2];
rz(-2.5640326) q[2];
rz(-1.1497568) q[3];
sx q[3];
rz(-1.3883608) q[3];
sx q[3];
rz(2.4770588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.98786551) q[0];
sx q[0];
rz(-2.5550714) q[0];
sx q[0];
rz(2.7541449) q[0];
rz(-2.2024343) q[1];
sx q[1];
rz(-2.1444131) q[1];
sx q[1];
rz(-1.4025677) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5706022) q[0];
sx q[0];
rz(-1.5413069) q[0];
sx q[0];
rz(1.5667314) q[0];
rz(-pi) q[1];
rz(-1.7895133) q[2];
sx q[2];
rz(-2.1195076) q[2];
sx q[2];
rz(2.1825841) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2991997) q[1];
sx q[1];
rz(-1.0789011) q[1];
sx q[1];
rz(-1.1115587) q[1];
rz(-3.0121147) q[3];
sx q[3];
rz(-0.63614142) q[3];
sx q[3];
rz(-2.5667218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42276057) q[2];
sx q[2];
rz(-1.3964802) q[2];
sx q[2];
rz(2.823901) q[2];
rz(2.9348532) q[3];
sx q[3];
rz(-2.5419149) q[3];
sx q[3];
rz(-2.3247705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3690255) q[0];
sx q[0];
rz(-1.439753) q[0];
sx q[0];
rz(-1.7279708) q[0];
rz(-0.47779045) q[1];
sx q[1];
rz(-1.3505892) q[1];
sx q[1];
rz(-0.40107045) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51940489) q[0];
sx q[0];
rz(-1.6739968) q[0];
sx q[0];
rz(-1.897057) q[0];
x q[1];
rz(2.4972649) q[2];
sx q[2];
rz(-1.4187078) q[2];
sx q[2];
rz(-2.6459141) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.229278) q[1];
sx q[1];
rz(-2.4114128) q[1];
sx q[1];
rz(-0.0095403949) q[1];
rz(2.0173666) q[3];
sx q[3];
rz(-2.7105769) q[3];
sx q[3];
rz(2.9197846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.59427375) q[2];
sx q[2];
rz(-1.6197562) q[2];
sx q[2];
rz(2.5857914) q[2];
rz(0.9764955) q[3];
sx q[3];
rz(-2.591811) q[3];
sx q[3];
rz(0.78021375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34898409) q[0];
sx q[0];
rz(-1.6225092) q[0];
sx q[0];
rz(-1.697631) q[0];
rz(-1.5199039) q[1];
sx q[1];
rz(-0.65602055) q[1];
sx q[1];
rz(2.8881883) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69617535) q[0];
sx q[0];
rz(-1.4679969) q[0];
sx q[0];
rz(2.1131383) q[0];
rz(1.6035945) q[2];
sx q[2];
rz(-2.6555853) q[2];
sx q[2];
rz(-0.82644586) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1052103) q[1];
sx q[1];
rz(-1.6670334) q[1];
sx q[1];
rz(0.63016816) q[1];
x q[2];
rz(-0.13417379) q[3];
sx q[3];
rz(-2.4878256) q[3];
sx q[3];
rz(-1.8133481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0306586) q[2];
sx q[2];
rz(-1.3867644) q[2];
sx q[2];
rz(-2.3542662) q[2];
rz(-0.91286719) q[3];
sx q[3];
rz(-0.74936167) q[3];
sx q[3];
rz(2.1319938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71516365) q[0];
sx q[0];
rz(-0.63868317) q[0];
sx q[0];
rz(3.0786247) q[0];
rz(3.0175623) q[1];
sx q[1];
rz(-0.80563671) q[1];
sx q[1];
rz(-2.6834992) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3375181) q[0];
sx q[0];
rz(-0.44263698) q[0];
sx q[0];
rz(-0.0054365693) q[0];
rz(-0.45256726) q[2];
sx q[2];
rz(-2.1308225) q[2];
sx q[2];
rz(2.3134311) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8372247) q[1];
sx q[1];
rz(-0.75290426) q[1];
sx q[1];
rz(-1.3992845) q[1];
x q[2];
rz(-1.0701418) q[3];
sx q[3];
rz(-1.6541012) q[3];
sx q[3];
rz(-1.7617776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1725585) q[2];
sx q[2];
rz(-0.92323747) q[2];
sx q[2];
rz(2.9120973) q[2];
rz(-3.138792) q[3];
sx q[3];
rz(-0.87001785) q[3];
sx q[3];
rz(1.8026479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.5795508) q[0];
sx q[0];
rz(-2.8631449) q[0];
sx q[0];
rz(-2.2221185) q[0];
rz(-3.0793076) q[1];
sx q[1];
rz(-1.0039763) q[1];
sx q[1];
rz(1.8744291) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8454682) q[0];
sx q[0];
rz(-2.2757029) q[0];
sx q[0];
rz(1.0494997) q[0];
rz(-1.7153347) q[2];
sx q[2];
rz(-0.83206165) q[2];
sx q[2];
rz(-1.9369672) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3082723) q[1];
sx q[1];
rz(-0.38494021) q[1];
sx q[1];
rz(-1.0233378) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8980428) q[3];
sx q[3];
rz(-2.0381513) q[3];
sx q[3];
rz(1.5598284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1798114) q[2];
sx q[2];
rz(-2.2634025) q[2];
sx q[2];
rz(2.5578257) q[2];
rz(-2.4328655) q[3];
sx q[3];
rz(-1.830359) q[3];
sx q[3];
rz(-0.023199737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0631183) q[0];
sx q[0];
rz(-0.45409504) q[0];
sx q[0];
rz(1.0725347) q[0];
rz(0.5468927) q[1];
sx q[1];
rz(-1.2439589) q[1];
sx q[1];
rz(-2.0297208) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5522963) q[0];
sx q[0];
rz(-2.6483905) q[0];
sx q[0];
rz(-1.769355) q[0];
x q[1];
rz(-1.215559) q[2];
sx q[2];
rz(-1.0676427) q[2];
sx q[2];
rz(-0.1728729) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6560981) q[1];
sx q[1];
rz(-1.6472677) q[1];
sx q[1];
rz(-0.2904201) q[1];
rz(-pi) q[2];
rz(-2.0734378) q[3];
sx q[3];
rz(-1.2253237) q[3];
sx q[3];
rz(0.84514602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.83773461) q[2];
sx q[2];
rz(-1.6656817) q[2];
sx q[2];
rz(-1.973935) q[2];
rz(1.6052823) q[3];
sx q[3];
rz(-1.4669908) q[3];
sx q[3];
rz(1.4060098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7325571) q[0];
sx q[0];
rz(-1.5719825) q[0];
sx q[0];
rz(-0.73079601) q[0];
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
rz(-2.9478361) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63379143) q[2];
sx q[2];
rz(-0.83665028) q[2];
sx q[2];
rz(1.6939236) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.054143993) q[1];
sx q[1];
rz(-2.1463697) q[1];
sx q[1];
rz(-1.7051484) q[1];
rz(-0.49236492) q[3];
sx q[3];
rz(-2.1015321) q[3];
sx q[3];
rz(2.9363971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.76453152) q[2];
sx q[2];
rz(-1.3724644) q[2];
sx q[2];
rz(-0.63684741) q[2];
rz(0.26646715) q[3];
sx q[3];
rz(-2.0839432) q[3];
sx q[3];
rz(1.5554957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4144142) q[0];
sx q[0];
rz(-2.0120912) q[0];
sx q[0];
rz(1.138858) q[0];
rz(0.75421929) q[1];
sx q[1];
rz(-2.8051839) q[1];
sx q[1];
rz(3.1220904) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0560023) q[0];
sx q[0];
rz(-1.3645571) q[0];
sx q[0];
rz(0.081736728) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15140622) q[2];
sx q[2];
rz(-1.3726241) q[2];
sx q[2];
rz(2.1422447) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6448977) q[1];
sx q[1];
rz(-0.95458191) q[1];
sx q[1];
rz(0.22209054) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9762514) q[3];
sx q[3];
rz(-1.8792361) q[3];
sx q[3];
rz(-0.43499085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1372244) q[2];
sx q[2];
rz(-1.4164111) q[2];
sx q[2];
rz(-0.84890378) q[2];
rz(-0.38765872) q[3];
sx q[3];
rz(-2.013423) q[3];
sx q[3];
rz(1.5415812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7983109) q[0];
sx q[0];
rz(-2.9738975) q[0];
sx q[0];
rz(-2.6570901) q[0];
rz(1.3867406) q[1];
sx q[1];
rz(-1.4258899) q[1];
sx q[1];
rz(-1.1482931) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1295706) q[0];
sx q[0];
rz(-2.9992933) q[0];
sx q[0];
rz(2.9905031) q[0];
x q[1];
rz(2.8414367) q[2];
sx q[2];
rz(-2.488392) q[2];
sx q[2];
rz(0.73210994) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3155568) q[1];
sx q[1];
rz(-0.6912187) q[1];
sx q[1];
rz(1.300699) q[1];
x q[2];
rz(-0.14264543) q[3];
sx q[3];
rz(-2.7890165) q[3];
sx q[3];
rz(2.8838317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.91223532) q[2];
sx q[2];
rz(-1.2939913) q[2];
sx q[2];
rz(-0.36515507) q[2];
rz(3.0129516) q[3];
sx q[3];
rz(-1.235685) q[3];
sx q[3];
rz(2.685759) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1098332) q[0];
sx q[0];
rz(-0.84072996) q[0];
sx q[0];
rz(1.6054556) q[0];
rz(2.1784492) q[1];
sx q[1];
rz(-1.2711202) q[1];
sx q[1];
rz(-1.0585379) q[1];
rz(-0.66394304) q[2];
sx q[2];
rz(-1.2524111) q[2];
sx q[2];
rz(-1.8128957) q[2];
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