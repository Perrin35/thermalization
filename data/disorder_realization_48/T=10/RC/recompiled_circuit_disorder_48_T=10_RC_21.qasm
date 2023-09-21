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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89361184) q[0];
sx q[0];
rz(-2.3900744) q[0];
sx q[0];
rz(0.052608842) q[0];
rz(-1.6115509) q[2];
sx q[2];
rz(-2.0353122) q[2];
sx q[2];
rz(-2.0704616) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7114746) q[1];
sx q[1];
rz(-2.4415486) q[1];
sx q[1];
rz(-0.75573604) q[1];
rz(-0.18913194) q[3];
sx q[3];
rz(-1.7889708) q[3];
sx q[3];
rz(1.3355586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.618764) q[2];
sx q[2];
rz(-2.6601807) q[2];
sx q[2];
rz(2.5640326) q[2];
rz(1.9918359) q[3];
sx q[3];
rz(-1.3883608) q[3];
sx q[3];
rz(2.4770588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98786551) q[0];
sx q[0];
rz(-2.5550714) q[0];
sx q[0];
rz(2.7541449) q[0];
rz(2.2024343) q[1];
sx q[1];
rz(-0.99717957) q[1];
sx q[1];
rz(-1.4025677) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99968602) q[0];
sx q[0];
rz(-1.5667331) q[0];
sx q[0];
rz(-3.1121029) q[0];
rz(0.55949975) q[2];
sx q[2];
rz(-1.7569949) q[2];
sx q[2];
rz(0.49636832) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.84239292) q[1];
sx q[1];
rz(-2.0626915) q[1];
sx q[1];
rz(1.1115587) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4757231) q[3];
sx q[3];
rz(-2.2007696) q[3];
sx q[3];
rz(-2.7271987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42276057) q[2];
sx q[2];
rz(-1.3964802) q[2];
sx q[2];
rz(-2.823901) q[2];
rz(-0.20673949) q[3];
sx q[3];
rz(-2.5419149) q[3];
sx q[3];
rz(-2.3247705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
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
rz(-0.40107045) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51940489) q[0];
sx q[0];
rz(-1.6739968) q[0];
sx q[0];
rz(-1.897057) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3813854) q[2];
sx q[2];
rz(-0.93511287) q[2];
sx q[2];
rz(-0.96178255) q[2];
sx q[3];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
x q[2];
rz(1.1242261) q[3];
sx q[3];
rz(-0.43101573) q[3];
sx q[3];
rz(-0.22180804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.59427375) q[2];
sx q[2];
rz(-1.6197562) q[2];
sx q[2];
rz(-2.5857914) q[2];
rz(0.9764955) q[3];
sx q[3];
rz(-0.54978168) q[3];
sx q[3];
rz(2.3613789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34898409) q[0];
sx q[0];
rz(-1.6225092) q[0];
sx q[0];
rz(-1.697631) q[0];
rz(-1.6216888) q[1];
sx q[1];
rz(-0.65602055) q[1];
sx q[1];
rz(0.25340432) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4454173) q[0];
sx q[0];
rz(-1.4679969) q[0];
sx q[0];
rz(1.0284543) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0850111) q[2];
sx q[2];
rz(-1.5861142) q[2];
sx q[2];
rz(2.3682396) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6771486) q[1];
sx q[1];
rz(-0.943999) q[1];
sx q[1];
rz(1.4518751) q[1];
rz(-0.6494135) q[3];
sx q[3];
rz(-1.6522437) q[3];
sx q[3];
rz(2.7922975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.110934) q[2];
sx q[2];
rz(-1.3867644) q[2];
sx q[2];
rz(-2.3542662) q[2];
rz(2.2287255) q[3];
sx q[3];
rz(-0.74936167) q[3];
sx q[3];
rz(2.1319938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.426429) q[0];
sx q[0];
rz(-0.63868317) q[0];
sx q[0];
rz(-0.062967904) q[0];
rz(3.0175623) q[1];
sx q[1];
rz(-0.80563671) q[1];
sx q[1];
rz(-2.6834992) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8040745) q[0];
sx q[0];
rz(-0.44263698) q[0];
sx q[0];
rz(-0.0054365693) q[0];
rz(-pi) q[1];
rz(-2.1797921) q[2];
sx q[2];
rz(-0.70448175) q[2];
sx q[2];
rz(-1.5693762) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.304368) q[1];
sx q[1];
rz(-2.3886884) q[1];
sx q[1];
rz(1.3992845) q[1];
rz(-1.0701418) q[3];
sx q[3];
rz(-1.4874914) q[3];
sx q[3];
rz(-1.3798151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1725585) q[2];
sx q[2];
rz(-2.2183552) q[2];
sx q[2];
rz(-2.9120973) q[2];
rz(3.138792) q[3];
sx q[3];
rz(-2.2715748) q[3];
sx q[3];
rz(1.8026479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(-0.5795508) q[0];
sx q[0];
rz(-0.27844772) q[0];
sx q[0];
rz(-2.2221185) q[0];
rz(-3.0793076) q[1];
sx q[1];
rz(-2.1376164) q[1];
sx q[1];
rz(1.2671635) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7125268) q[0];
sx q[0];
rz(-0.8493087) q[0];
sx q[0];
rz(0.52961403) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15684776) q[2];
sx q[2];
rz(-0.75011293) q[2];
sx q[2];
rz(1.4175121) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.25176469) q[1];
sx q[1];
rz(-1.3740731) q[1];
sx q[1];
rz(-1.2377435) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1250161) q[3];
sx q[3];
rz(-0.52281724) q[3];
sx q[3];
rz(1.0558053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1798114) q[2];
sx q[2];
rz(-0.87819019) q[2];
sx q[2];
rz(-2.5578257) q[2];
rz(-2.4328655) q[3];
sx q[3];
rz(-1.3112336) q[3];
sx q[3];
rz(0.023199737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0631183) q[0];
sx q[0];
rz(-2.6874976) q[0];
sx q[0];
rz(-2.069058) q[0];
rz(-0.5468927) q[1];
sx q[1];
rz(-1.8976338) q[1];
sx q[1];
rz(1.1118719) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3355013) q[0];
sx q[0];
rz(-1.4772692) q[0];
sx q[0];
rz(2.0557687) q[0];
rz(-pi) q[1];
rz(2.6107437) q[2];
sx q[2];
rz(-1.2611654) q[2];
sx q[2];
rz(-1.9206778) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6560981) q[1];
sx q[1];
rz(-1.6472677) q[1];
sx q[1];
rz(0.2904201) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38970077) q[3];
sx q[3];
rz(-2.0412363) q[3];
sx q[3];
rz(-0.54159347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.303858) q[2];
sx q[2];
rz(-1.4759109) q[2];
sx q[2];
rz(1.1676577) q[2];
rz(-1.5363103) q[3];
sx q[3];
rz(-1.4669908) q[3];
sx q[3];
rz(-1.7355828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7325571) q[0];
sx q[0];
rz(-1.5696101) q[0];
sx q[0];
rz(2.4107966) q[0];
rz(-2.2413975) q[1];
sx q[1];
rz(-2.3370445) q[1];
sx q[1];
rz(2.3866167) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1599931) q[0];
sx q[0];
rz(-2.2845075) q[0];
sx q[0];
rz(-0.1937565) q[0];
rz(-pi) q[1];
x q[1];
rz(2.412699) q[2];
sx q[2];
rz(-1.1155827) q[2];
sx q[2];
rz(-2.560937) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.29754408) q[1];
sx q[1];
rz(-0.58931671) q[1];
sx q[1];
rz(-2.938016) q[1];
x q[2];
rz(-2.2488392) q[3];
sx q[3];
rz(-0.70749456) q[3];
sx q[3];
rz(-1.0196109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3770611) q[2];
sx q[2];
rz(-1.7691282) q[2];
sx q[2];
rz(2.5047452) q[2];
rz(0.26646715) q[3];
sx q[3];
rz(-2.0839432) q[3];
sx q[3];
rz(1.5554957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.4144142) q[0];
sx q[0];
rz(-1.1295015) q[0];
sx q[0];
rz(1.138858) q[0];
rz(2.3873734) q[1];
sx q[1];
rz(-2.8051839) q[1];
sx q[1];
rz(-3.1220904) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50197983) q[0];
sx q[0];
rz(-1.6507971) q[0];
sx q[0];
rz(1.3638858) q[0];
rz(0.15140622) q[2];
sx q[2];
rz(-1.7689686) q[2];
sx q[2];
rz(2.1422447) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1972678) q[1];
sx q[1];
rz(-1.3900583) q[1];
sx q[1];
rz(0.9428057) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2583624) q[3];
sx q[3];
rz(-1.4133246) q[3];
sx q[3];
rz(-2.0563994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1372244) q[2];
sx q[2];
rz(-1.4164111) q[2];
sx q[2];
rz(2.2926889) q[2];
rz(-0.38765872) q[3];
sx q[3];
rz(-1.1281697) q[3];
sx q[3];
rz(-1.5415812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3432817) q[0];
sx q[0];
rz(-2.9738975) q[0];
sx q[0];
rz(-0.48450255) q[0];
rz(1.7548521) q[1];
sx q[1];
rz(-1.7157028) q[1];
sx q[1];
rz(-1.1482931) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1185547) q[0];
sx q[0];
rz(-1.430129) q[0];
sx q[0];
rz(1.592357) q[0];
rz(-pi) q[1];
rz(-1.3482773) q[2];
sx q[2];
rz(-0.95138022) q[2];
sx q[2];
rz(-0.36048181) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3155568) q[1];
sx q[1];
rz(-0.6912187) q[1];
sx q[1];
rz(-1.8408937) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5185353) q[3];
sx q[3];
rz(-1.2219547) q[3];
sx q[3];
rz(-2.7319752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.91223532) q[2];
sx q[2];
rz(-1.2939913) q[2];
sx q[2];
rz(2.7764376) q[2];
rz(-3.0129516) q[3];
sx q[3];
rz(-1.235685) q[3];
sx q[3];
rz(0.45583367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1098332) q[0];
sx q[0];
rz(-2.3008627) q[0];
sx q[0];
rz(-1.536137) q[0];
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
rz(-0.036988463) q[3];
sx q[3];
rz(-2.130641) q[3];
sx q[3];
rz(-0.11161042) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];