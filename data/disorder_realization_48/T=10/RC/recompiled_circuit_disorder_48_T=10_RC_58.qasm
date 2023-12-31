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
rz(-0.45733991) q[1];
sx q[1];
rz(5.338905) q[1];
sx q[1];
rz(10.64325) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96555644) q[0];
sx q[0];
rz(-2.3210225) q[0];
sx q[0];
rz(1.6198938) q[0];
rz(-pi) q[1];
rz(0.081131447) q[2];
sx q[2];
rz(-0.46617026) q[2];
sx q[2];
rz(-0.98035882) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.76518607) q[1];
sx q[1];
rz(-2.0284488) q[1];
sx q[1];
rz(-0.54995579) q[1];
rz(-pi) q[2];
rz(0.18913194) q[3];
sx q[3];
rz(-1.3526219) q[3];
sx q[3];
rz(1.3355586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.618764) q[2];
sx q[2];
rz(-0.4814119) q[2];
sx q[2];
rz(-0.5775601) q[2];
rz(1.9918359) q[3];
sx q[3];
rz(-1.7532319) q[3];
sx q[3];
rz(-2.4770588) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1537271) q[0];
sx q[0];
rz(-2.5550714) q[0];
sx q[0];
rz(2.7541449) q[0];
rz(-2.2024343) q[1];
sx q[1];
rz(-0.99717957) q[1];
sx q[1];
rz(1.4025677) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1419066) q[0];
sx q[0];
rz(-1.5667331) q[0];
sx q[0];
rz(3.1121029) q[0];
x q[1];
rz(-2.5820929) q[2];
sx q[2];
rz(-1.3845978) q[2];
sx q[2];
rz(-0.49636832) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1837511) q[1];
sx q[1];
rz(-1.1693923) q[1];
sx q[1];
rz(-2.6028231) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4757231) q[3];
sx q[3];
rz(-0.94082309) q[3];
sx q[3];
rz(-2.7271987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.42276057) q[2];
sx q[2];
rz(-1.3964802) q[2];
sx q[2];
rz(2.823901) q[2];
rz(-0.20673949) q[3];
sx q[3];
rz(-2.5419149) q[3];
sx q[3];
rz(0.81682214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7725672) q[0];
sx q[0];
rz(-1.439753) q[0];
sx q[0];
rz(1.4136219) q[0];
rz(0.47779045) q[1];
sx q[1];
rz(-1.7910035) q[1];
sx q[1];
rz(-0.40107045) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51940489) q[0];
sx q[0];
rz(-1.4675958) q[0];
sx q[0];
rz(-1.2445356) q[0];
rz(-pi) q[1];
rz(-1.7602073) q[2];
sx q[2];
rz(-0.93511287) q[2];
sx q[2];
rz(-2.1798101) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2420826) q[1];
sx q[1];
rz(-2.3009355) q[1];
sx q[1];
rz(-1.5622557) q[1];
x q[2];
rz(-1.1776351) q[3];
sx q[3];
rz(-1.3893681) q[3];
sx q[3];
rz(-1.7593311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.59427375) q[2];
sx q[2];
rz(-1.6197562) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7926086) q[0];
sx q[0];
rz(-1.5190834) q[0];
sx q[0];
rz(1.4439616) q[0];
rz(-1.5199039) q[1];
sx q[1];
rz(-0.65602055) q[1];
sx q[1];
rz(2.8881883) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4454173) q[0];
sx q[0];
rz(-1.6735958) q[0];
sx q[0];
rz(1.0284543) q[0];
x q[1];
rz(-0.01732145) q[2];
sx q[2];
rz(-2.0565196) q[2];
sx q[2];
rz(-0.7893562) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.464444) q[1];
sx q[1];
rz(-0.943999) q[1];
sx q[1];
rz(1.4518751) q[1];
rz(-0.13417379) q[3];
sx q[3];
rz(-0.65376702) q[3];
sx q[3];
rz(-1.3282446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.110934) q[2];
sx q[2];
rz(-1.3867644) q[2];
sx q[2];
rz(0.78732642) q[2];
rz(-2.2287255) q[3];
sx q[3];
rz(-0.74936167) q[3];
sx q[3];
rz(-2.1319938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71516365) q[0];
sx q[0];
rz(-0.63868317) q[0];
sx q[0];
rz(3.0786247) q[0];
rz(-0.12403034) q[1];
sx q[1];
rz(-0.80563671) q[1];
sx q[1];
rz(-2.6834992) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3315017) q[0];
sx q[0];
rz(-2.0134263) q[0];
sx q[0];
rz(-1.5733733) q[0];
rz(-2.6890254) q[2];
sx q[2];
rz(-2.1308225) q[2];
sx q[2];
rz(0.82816154) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.071306989) q[1];
sx q[1];
rz(-2.3100393) q[1];
sx q[1];
rz(2.9830095) q[1];
rz(-pi) q[2];
rz(0.094893806) q[3];
sx q[3];
rz(-2.0695544) q[3];
sx q[3];
rz(-2.9961078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9690341) q[2];
sx q[2];
rz(-2.2183552) q[2];
sx q[2];
rz(2.9120973) q[2];
rz(-0.0028006639) q[3];
sx q[3];
rz(-2.2715748) q[3];
sx q[3];
rz(1.8026479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5795508) q[0];
sx q[0];
rz(-2.8631449) q[0];
sx q[0];
rz(2.2221185) q[0];
rz(-3.0793076) q[1];
sx q[1];
rz(-1.0039763) q[1];
sx q[1];
rz(1.8744291) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6309109) q[0];
sx q[0];
rz(-1.1817389) q[0];
sx q[0];
rz(-2.3657777) q[0];
rz(-pi) q[1];
rz(-0.15684776) q[2];
sx q[2];
rz(-0.75011293) q[2];
sx q[2];
rz(1.7240806) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.889828) q[1];
sx q[1];
rz(-1.7675195) q[1];
sx q[1];
rz(-1.2377435) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8980428) q[3];
sx q[3];
rz(-1.1034414) q[3];
sx q[3];
rz(1.5817643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.9617812) q[2];
sx q[2];
rz(-0.87819019) q[2];
sx q[2];
rz(0.58376694) q[2];
rz(-0.70872712) q[3];
sx q[3];
rz(-1.3112336) q[3];
sx q[3];
rz(3.1183929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07847438) q[0];
sx q[0];
rz(-0.45409504) q[0];
sx q[0];
rz(-2.069058) q[0];
rz(2.5947) q[1];
sx q[1];
rz(-1.2439589) q[1];
sx q[1];
rz(-1.1118719) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81387732) q[0];
sx q[0];
rz(-2.0534671) q[0];
sx q[0];
rz(3.0359603) q[0];
rz(-pi) q[1];
rz(-2.578031) q[2];
sx q[2];
rz(-0.60699082) q[2];
sx q[2];
rz(-2.3129472) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0791196) q[1];
sx q[1];
rz(-1.8603431) q[1];
sx q[1];
rz(-1.6505961) q[1];
rz(-0.92915793) q[3];
sx q[3];
rz(-2.540179) q[3];
sx q[3];
rz(-1.2777559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.303858) q[2];
sx q[2];
rz(-1.6656817) q[2];
sx q[2];
rz(-1.1676577) q[2];
rz(1.6052823) q[3];
sx q[3];
rz(-1.6746018) q[3];
sx q[3];
rz(1.7355828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4090356) q[0];
sx q[0];
rz(-1.5719825) q[0];
sx q[0];
rz(-0.73079601) q[0];
rz(2.2413975) q[1];
sx q[1];
rz(-0.80454818) q[1];
sx q[1];
rz(-0.75497595) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1599931) q[0];
sx q[0];
rz(-2.2845075) q[0];
sx q[0];
rz(2.9478361) q[0];
rz(-pi) q[1];
rz(0.63379143) q[2];
sx q[2];
rz(-0.83665028) q[2];
sx q[2];
rz(1.447669) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0874487) q[1];
sx q[1];
rz(-0.99522299) q[1];
sx q[1];
rz(-1.7051484) q[1];
x q[2];
rz(2.6492277) q[3];
sx q[3];
rz(-1.0400606) q[3];
sx q[3];
rz(-2.9363971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.76453152) q[2];
sx q[2];
rz(-1.3724644) q[2];
sx q[2];
rz(-2.5047452) q[2];
rz(2.8751255) q[3];
sx q[3];
rz(-1.0576495) q[3];
sx q[3];
rz(-1.586097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4144142) q[0];
sx q[0];
rz(-1.1295015) q[0];
sx q[0];
rz(2.0027347) q[0];
rz(-2.3873734) q[1];
sx q[1];
rz(-0.33640877) q[1];
sx q[1];
rz(0.019502217) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50197983) q[0];
sx q[0];
rz(-1.4907955) q[0];
sx q[0];
rz(1.3638858) q[0];
rz(-pi) q[1];
rz(-0.15140622) q[2];
sx q[2];
rz(-1.3726241) q[2];
sx q[2];
rz(-0.99934794) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2724185) q[1];
sx q[1];
rz(-2.4915016) q[1];
sx q[1];
rz(1.8723349) q[1];
rz(-pi) q[2];
rz(-2.0476258) q[3];
sx q[3];
rz(-0.34871021) q[3];
sx q[3];
rz(-0.93758067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.0043682178) q[2];
sx q[2];
rz(-1.7251816) q[2];
sx q[2];
rz(2.2926889) q[2];
rz(-0.38765872) q[3];
sx q[3];
rz(-1.1281697) q[3];
sx q[3];
rz(1.6000115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(1.7983109) q[0];
sx q[0];
rz(-0.16769519) q[0];
sx q[0];
rz(0.48450255) q[0];
rz(-1.3867406) q[1];
sx q[1];
rz(-1.7157028) q[1];
sx q[1];
rz(1.9932995) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023037993) q[0];
sx q[0];
rz(-1.430129) q[0];
sx q[0];
rz(1.5492357) q[0];
rz(-pi) q[1];
rz(-1.3482773) q[2];
sx q[2];
rz(-2.1902124) q[2];
sx q[2];
rz(0.36048181) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.53459586) q[1];
sx q[1];
rz(-1.3998704) q[1];
sx q[1];
rz(-0.89764885) q[1];
rz(-pi) q[2];
rz(-1.6230574) q[3];
sx q[3];
rz(-1.9196379) q[3];
sx q[3];
rz(0.40961743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2293573) q[2];
sx q[2];
rz(-1.2939913) q[2];
sx q[2];
rz(2.7764376) q[2];
rz(-3.0129516) q[3];
sx q[3];
rz(-1.9059076) q[3];
sx q[3];
rz(2.685759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(3.0317595) q[0];
sx q[0];
rz(-2.3008627) q[0];
sx q[0];
rz(-1.536137) q[0];
rz(-2.1784492) q[1];
sx q[1];
rz(-1.8704725) q[1];
sx q[1];
rz(2.0830547) q[1];
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
