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
rz(3.9324023) q[0];
sx q[0];
rz(12.232236) q[0];
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
rz(0.63872913) q[0];
sx q[0];
rz(-1.5348866) q[0];
sx q[0];
rz(0.7508276) q[0];
x q[1];
rz(-0.46484868) q[2];
sx q[2];
rz(-1.6072304) q[2];
sx q[2];
rz(-2.6236617) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.54108809) q[1];
sx q[1];
rz(-2.0588015) q[1];
sx q[1];
rz(1.0469251) q[1];
rz(-pi) q[2];
rz(0.18913194) q[3];
sx q[3];
rz(-1.7889708) q[3];
sx q[3];
rz(1.8060341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.618764) q[2];
sx q[2];
rz(-0.4814119) q[2];
sx q[2];
rz(0.5775601) q[2];
rz(1.1497568) q[3];
sx q[3];
rz(-1.7532319) q[3];
sx q[3];
rz(2.4770588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1537271) q[0];
sx q[0];
rz(-0.58652121) q[0];
sx q[0];
rz(-0.38744774) q[0];
rz(-2.2024343) q[1];
sx q[1];
rz(-0.99717957) q[1];
sx q[1];
rz(-1.739025) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1419066) q[0];
sx q[0];
rz(-1.5667331) q[0];
sx q[0];
rz(3.1121029) q[0];
rz(-pi) q[1];
rz(-2.5820929) q[2];
sx q[2];
rz(-1.3845978) q[2];
sx q[2];
rz(-0.49636832) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2991997) q[1];
sx q[1];
rz(-2.0626915) q[1];
sx q[1];
rz(-2.030034) q[1];
rz(0.6321294) q[3];
sx q[3];
rz(-1.4940133) q[3];
sx q[3];
rz(2.0413105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7188321) q[2];
sx q[2];
rz(-1.7451124) q[2];
sx q[2];
rz(-0.31769162) q[2];
rz(2.9348532) q[3];
sx q[3];
rz(-0.59967774) q[3];
sx q[3];
rz(-0.81682214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7725672) q[0];
sx q[0];
rz(-1.439753) q[0];
sx q[0];
rz(1.7279708) q[0];
rz(-0.47779045) q[1];
sx q[1];
rz(-1.7910035) q[1];
sx q[1];
rz(0.40107045) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7558407) q[0];
sx q[0];
rz(-2.7999561) q[0];
sx q[0];
rz(1.2582448) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4972649) q[2];
sx q[2];
rz(-1.4187078) q[2];
sx q[2];
rz(-2.6459141) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8071825) q[1];
sx q[1];
rz(-1.5644329) q[1];
sx q[1];
rz(-0.73015726) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19604711) q[3];
sx q[3];
rz(-1.184433) q[3];
sx q[3];
rz(-0.26323174) q[3];
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
rz(-2.3613789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34898409) q[0];
sx q[0];
rz(-1.5190834) q[0];
sx q[0];
rz(-1.697631) q[0];
rz(1.6216888) q[1];
sx q[1];
rz(-0.65602055) q[1];
sx q[1];
rz(2.8881883) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4454173) q[0];
sx q[0];
rz(-1.4679969) q[0];
sx q[0];
rz(2.1131383) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5379982) q[2];
sx q[2];
rz(-2.6555853) q[2];
sx q[2];
rz(-2.3151468) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.464444) q[1];
sx q[1];
rz(-2.1975937) q[1];
sx q[1];
rz(-1.6897175) q[1];
x q[2];
rz(1.4686618) q[3];
sx q[3];
rz(-2.2176952) q[3];
sx q[3];
rz(1.981786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.110934) q[2];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.426429) q[0];
sx q[0];
rz(-2.5029095) q[0];
sx q[0];
rz(-0.062967904) q[0];
rz(-0.12403034) q[1];
sx q[1];
rz(-2.3359559) q[1];
sx q[1];
rz(2.6834992) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3375181) q[0];
sx q[0];
rz(-2.6989557) q[0];
sx q[0];
rz(3.1361561) q[0];
rz(-2.1797921) q[2];
sx q[2];
rz(-2.4371109) q[2];
sx q[2];
rz(1.5693762) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.749436) q[1];
sx q[1];
rz(-1.4538308) q[1];
sx q[1];
rz(2.3163296) q[1];
rz(-1.3985653) q[3];
sx q[3];
rz(-0.50695626) q[3];
sx q[3];
rz(-2.7996922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1725585) q[2];
sx q[2];
rz(-0.92323747) q[2];
sx q[2];
rz(0.22949533) q[2];
rz(3.138792) q[3];
sx q[3];
rz(-0.87001785) q[3];
sx q[3];
rz(1.3389448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5795508) q[0];
sx q[0];
rz(-0.27844772) q[0];
sx q[0];
rz(2.2221185) q[0];
rz(-3.0793076) q[1];
sx q[1];
rz(-2.1376164) q[1];
sx q[1];
rz(1.2671635) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5106817) q[0];
sx q[0];
rz(-1.9598538) q[0];
sx q[0];
rz(-2.3657777) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7153347) q[2];
sx q[2];
rz(-2.309531) q[2];
sx q[2];
rz(1.2046255) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2515182) q[1];
sx q[1];
rz(-1.2444082) q[1];
sx q[1];
rz(-0.20784394) q[1];
rz(-pi) q[2];
rz(2.8980428) q[3];
sx q[3];
rz(-2.0381513) q[3];
sx q[3];
rz(-1.5817643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.9617812) q[2];
sx q[2];
rz(-2.2634025) q[2];
sx q[2];
rz(0.58376694) q[2];
rz(0.70872712) q[3];
sx q[3];
rz(-1.830359) q[3];
sx q[3];
rz(-0.023199737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07847438) q[0];
sx q[0];
rz(-2.6874976) q[0];
sx q[0];
rz(-2.069058) q[0];
rz(0.5468927) q[1];
sx q[1];
rz(-1.2439589) q[1];
sx q[1];
rz(1.1118719) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81387732) q[0];
sx q[0];
rz(-2.0534671) q[0];
sx q[0];
rz(-3.0359603) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6107437) q[2];
sx q[2];
rz(-1.2611654) q[2];
sx q[2];
rz(-1.2209148) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8060311) q[1];
sx q[1];
rz(-0.30004382) q[1];
sx q[1];
rz(2.8801444) q[1];
x q[2];
rz(2.7518919) q[3];
sx q[3];
rz(-2.0412363) q[3];
sx q[3];
rz(-0.54159347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.303858) q[2];
sx q[2];
rz(-1.6656817) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4090356) q[0];
sx q[0];
rz(-1.5696101) q[0];
sx q[0];
rz(-2.4107966) q[0];
rz(2.2413975) q[1];
sx q[1];
rz(-2.3370445) q[1];
sx q[1];
rz(0.75497595) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1599931) q[0];
sx q[0];
rz(-0.85708517) q[0];
sx q[0];
rz(-2.9478361) q[0];
x q[1];
rz(0.99004284) q[2];
sx q[2];
rz(-0.92951894) q[2];
sx q[2];
rz(-2.5255447) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6983812) q[1];
sx q[1];
rz(-1.6834007) q[1];
sx q[1];
rz(-0.57971445) q[1];
x q[2];
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
rz(2.3770611) q[2];
sx q[2];
rz(-1.3724644) q[2];
sx q[2];
rz(2.5047452) q[2];
rz(-2.8751255) q[3];
sx q[3];
rz(-2.0839432) q[3];
sx q[3];
rz(1.5554957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4144142) q[0];
sx q[0];
rz(-2.0120912) q[0];
sx q[0];
rz(2.0027347) q[0];
rz(-0.75421929) q[1];
sx q[1];
rz(-0.33640877) q[1];
sx q[1];
rz(3.1220904) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0560023) q[0];
sx q[0];
rz(-1.7770355) q[0];
sx q[0];
rz(-3.0598559) q[0];
rz(2.2150061) q[2];
sx q[2];
rz(-2.8928061) q[2];
sx q[2];
rz(-2.8015346) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2724185) q[1];
sx q[1];
rz(-2.4915016) q[1];
sx q[1];
rz(-1.2692578) q[1];
rz(-pi) q[2];
rz(-2.0476258) q[3];
sx q[3];
rz(-0.34871021) q[3];
sx q[3];
rz(2.204012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.0043682178) q[2];
sx q[2];
rz(-1.4164111) q[2];
sx q[2];
rz(-0.84890378) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7983109) q[0];
sx q[0];
rz(-2.9738975) q[0];
sx q[0];
rz(0.48450255) q[0];
rz(-1.3867406) q[1];
sx q[1];
rz(-1.7157028) q[1];
sx q[1];
rz(-1.1482931) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.590811) q[0];
sx q[0];
rz(-1.5494487) q[0];
sx q[0];
rz(-3.000893) q[0];
rz(-pi) q[1];
rz(-1.3482773) q[2];
sx q[2];
rz(-2.1902124) q[2];
sx q[2];
rz(0.36048181) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.53459586) q[1];
sx q[1];
rz(-1.3998704) q[1];
sx q[1];
rz(2.2439438) q[1];
rz(-0.34928068) q[3];
sx q[3];
rz(-1.619907) q[3];
sx q[3];
rz(1.1790566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91223532) q[2];
sx q[2];
rz(-1.8476013) q[2];
sx q[2];
rz(-0.36515507) q[2];
rz(-0.12864104) q[3];
sx q[3];
rz(-1.9059076) q[3];
sx q[3];
rz(-2.685759) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0317595) q[0];
sx q[0];
rz(-2.3008627) q[0];
sx q[0];
rz(-1.536137) q[0];
rz(0.96314349) q[1];
sx q[1];
rz(-1.8704725) q[1];
sx q[1];
rz(2.0830547) q[1];
rz(1.1744432) q[2];
sx q[2];
rz(-2.195993) q[2];
sx q[2];
rz(3.1396951) q[2];
rz(2.1309489) q[3];
sx q[3];
rz(-1.602136) q[3];
sx q[3];
rz(-1.7020561) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
