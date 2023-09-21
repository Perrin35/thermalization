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
rz(-0.9442803) q[1];
sx q[1];
rz(1.2184719) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5028635) q[0];
sx q[0];
rz(-1.5348866) q[0];
sx q[0];
rz(0.7508276) q[0];
rz(-pi) q[1];
rz(-2.676744) q[2];
sx q[2];
rz(-1.6072304) q[2];
sx q[2];
rz(-0.51793098) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7114746) q[1];
sx q[1];
rz(-0.700044) q[1];
sx q[1];
rz(2.3858566) q[1];
rz(-pi) q[2];
rz(-2.9524607) q[3];
sx q[3];
rz(-1.7889708) q[3];
sx q[3];
rz(1.8060341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5228287) q[2];
sx q[2];
rz(-2.6601807) q[2];
sx q[2];
rz(0.5775601) q[2];
rz(1.1497568) q[3];
sx q[3];
rz(-1.3883608) q[3];
sx q[3];
rz(-2.4770588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
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
rz(-2.1444131) q[1];
sx q[1];
rz(-1.4025677) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5706022) q[0];
sx q[0];
rz(-1.5413069) q[0];
sx q[0];
rz(1.5667314) q[0];
rz(-pi) q[1];
rz(-0.34105532) q[2];
sx q[2];
rz(-2.555072) q[2];
sx q[2];
rz(-1.7797433) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1837511) q[1];
sx q[1];
rz(-1.9722003) q[1];
sx q[1];
rz(-0.5387696) q[1];
x q[2];
rz(-0.12947793) q[3];
sx q[3];
rz(-2.5054512) q[3];
sx q[3];
rz(0.57487088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
sx q[1];
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
rz(-1.7725672) q[0];
sx q[0];
rz(-1.7018397) q[0];
sx q[0];
rz(1.4136219) q[0];
rz(2.6638022) q[1];
sx q[1];
rz(-1.3505892) q[1];
sx q[1];
rz(-0.40107045) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.385752) q[0];
sx q[0];
rz(-2.7999561) q[0];
sx q[0];
rz(1.2582448) q[0];
x q[1];
rz(-1.7602073) q[2];
sx q[2];
rz(-0.93511287) q[2];
sx q[2];
rz(-2.1798101) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.229278) q[1];
sx q[1];
rz(-0.73017987) q[1];
sx q[1];
rz(3.1320523) q[1];
rz(-pi) q[2];
rz(1.1776351) q[3];
sx q[3];
rz(-1.7522246) q[3];
sx q[3];
rz(-1.7593311) q[3];
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
rz(2.5857914) q[2];
rz(-0.9764955) q[3];
sx q[3];
rz(-2.591811) q[3];
sx q[3];
rz(-0.78021375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7926086) q[0];
sx q[0];
rz(-1.6225092) q[0];
sx q[0];
rz(-1.4439616) q[0];
rz(1.5199039) q[1];
sx q[1];
rz(-0.65602055) q[1];
sx q[1];
rz(-2.8881883) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70595104) q[0];
sx q[0];
rz(-2.5905529) q[0];
sx q[0];
rz(-1.7680697) q[0];
rz(3.1242712) q[2];
sx q[2];
rz(-2.0565196) q[2];
sx q[2];
rz(2.3522365) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6654012) q[1];
sx q[1];
rz(-0.63648495) q[1];
sx q[1];
rz(0.16237662) q[1];
rz(-pi) q[2];
rz(0.13417379) q[3];
sx q[3];
rz(-2.4878256) q[3];
sx q[3];
rz(1.8133481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0306586) q[2];
sx q[2];
rz(-1.3867644) q[2];
sx q[2];
rz(0.78732642) q[2];
rz(-2.2287255) q[3];
sx q[3];
rz(-2.392231) q[3];
sx q[3];
rz(-1.0095989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71516365) q[0];
sx q[0];
rz(-2.5029095) q[0];
sx q[0];
rz(-0.062967904) q[0];
rz(-3.0175623) q[1];
sx q[1];
rz(-2.3359559) q[1];
sx q[1];
rz(-2.6834992) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3375181) q[0];
sx q[0];
rz(-0.44263698) q[0];
sx q[0];
rz(-0.0054365693) q[0];
x q[1];
rz(-0.96180054) q[2];
sx q[2];
rz(-2.4371109) q[2];
sx q[2];
rz(1.5722164) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.071306989) q[1];
sx q[1];
rz(-0.83155337) q[1];
sx q[1];
rz(0.15858312) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3985653) q[3];
sx q[3];
rz(-2.6346364) q[3];
sx q[3];
rz(2.7996922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9690341) q[2];
sx q[2];
rz(-0.92323747) q[2];
sx q[2];
rz(-0.22949533) q[2];
rz(-3.138792) q[3];
sx q[3];
rz(-2.2715748) q[3];
sx q[3];
rz(-1.8026479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7125268) q[0];
sx q[0];
rz(-0.8493087) q[0];
sx q[0];
rz(2.6119786) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.426258) q[2];
sx q[2];
rz(-0.83206165) q[2];
sx q[2];
rz(-1.2046255) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.889828) q[1];
sx q[1];
rz(-1.7675195) q[1];
sx q[1];
rz(1.9038492) q[1];
rz(-2.0165765) q[3];
sx q[3];
rz(-0.52281724) q[3];
sx q[3];
rz(1.0558053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.9617812) q[2];
sx q[2];
rz(-2.2634025) q[2];
sx q[2];
rz(-0.58376694) q[2];
rz(0.70872712) q[3];
sx q[3];
rz(-1.830359) q[3];
sx q[3];
rz(3.1183929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0631183) q[0];
sx q[0];
rz(-2.6874976) q[0];
sx q[0];
rz(1.0725347) q[0];
rz(0.5468927) q[1];
sx q[1];
rz(-1.2439589) q[1];
sx q[1];
rz(1.1118719) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5522963) q[0];
sx q[0];
rz(-0.49320212) q[0];
sx q[0];
rz(1.769355) q[0];
rz(-pi) q[1];
rz(-0.53084897) q[2];
sx q[2];
rz(-1.2611654) q[2];
sx q[2];
rz(1.2209148) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.48549451) q[1];
sx q[1];
rz(-1.494325) q[1];
sx q[1];
rz(-0.2904201) q[1];
x q[2];
rz(2.7518919) q[3];
sx q[3];
rz(-2.0412363) q[3];
sx q[3];
rz(2.5999992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.83773461) q[2];
sx q[2];
rz(-1.6656817) q[2];
sx q[2];
rz(-1.1676577) q[2];
rz(-1.5363103) q[3];
sx q[3];
rz(-1.6746018) q[3];
sx q[3];
rz(-1.4060098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4090356) q[0];
sx q[0];
rz(-1.5696101) q[0];
sx q[0];
rz(2.4107966) q[0];
rz(-0.90019512) q[1];
sx q[1];
rz(-0.80454818) q[1];
sx q[1];
rz(-0.75497595) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6801493) q[0];
sx q[0];
rz(-1.4247243) q[0];
sx q[0];
rz(-2.2938674) q[0];
rz(-0.72889363) q[2];
sx q[2];
rz(-2.0260099) q[2];
sx q[2];
rz(2.560937) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8440486) q[1];
sx q[1];
rz(-0.58931671) q[1];
sx q[1];
rz(-2.938016) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98324361) q[3];
sx q[3];
rz(-1.1508815) q[3];
sx q[3];
rz(-2.0411232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3770611) q[2];
sx q[2];
rz(-1.3724644) q[2];
sx q[2];
rz(-0.63684741) q[2];
rz(-2.8751255) q[3];
sx q[3];
rz(-2.0839432) q[3];
sx q[3];
rz(-1.586097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72717845) q[0];
sx q[0];
rz(-2.0120912) q[0];
sx q[0];
rz(-1.138858) q[0];
rz(2.3873734) q[1];
sx q[1];
rz(-0.33640877) q[1];
sx q[1];
rz(-0.019502217) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4365387) q[0];
sx q[0];
rz(-2.9199613) q[0];
sx q[0];
rz(-1.9428695) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9901864) q[2];
sx q[2];
rz(-1.7689686) q[2];
sx q[2];
rz(2.1422447) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6448977) q[1];
sx q[1];
rz(-2.1870107) q[1];
sx q[1];
rz(2.9195021) q[1];
rz(-pi) q[2];
rz(-2.9762514) q[3];
sx q[3];
rz(-1.2623566) q[3];
sx q[3];
rz(2.7066018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0043682178) q[2];
sx q[2];
rz(-1.4164111) q[2];
sx q[2];
rz(0.84890378) q[2];
rz(-0.38765872) q[3];
sx q[3];
rz(-2.013423) q[3];
sx q[3];
rz(1.5415812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.4258899) q[1];
sx q[1];
rz(-1.9932995) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023037993) q[0];
sx q[0];
rz(-1.430129) q[0];
sx q[0];
rz(1.592357) q[0];
rz(0.63126385) q[2];
sx q[2];
rz(-1.751465) q[2];
sx q[2];
rz(1.0797015) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9705829) q[1];
sx q[1];
rz(-2.2323771) q[1];
sx q[1];
rz(0.2172825) q[1];
x q[2];
rz(-2.9989472) q[3];
sx q[3];
rz(-2.7890165) q[3];
sx q[3];
rz(0.25776097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2293573) q[2];
sx q[2];
rz(-1.8476013) q[2];
sx q[2];
rz(0.36515507) q[2];
rz(0.12864104) q[3];
sx q[3];
rz(-1.235685) q[3];
sx q[3];
rz(0.45583367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0317595) q[0];
sx q[0];
rz(-2.3008627) q[0];
sx q[0];
rz(-1.536137) q[0];
rz(2.1784492) q[1];
sx q[1];
rz(-1.2711202) q[1];
sx q[1];
rz(-1.0585379) q[1];
rz(0.49114901) q[2];
sx q[2];
rz(-2.4158203) q[2];
sx q[2];
rz(-0.62266785) q[2];
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
