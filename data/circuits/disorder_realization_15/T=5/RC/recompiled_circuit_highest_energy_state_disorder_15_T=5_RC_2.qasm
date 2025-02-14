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
rz(1.6871356) q[0];
sx q[0];
rz(-1.1633101) q[0];
sx q[0];
rz(1.8736725) q[0];
rz(-1.4866225) q[1];
sx q[1];
rz(-1.2868737) q[1];
sx q[1];
rz(-0.16965228) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96972966) q[0];
sx q[0];
rz(-2.5262882) q[0];
sx q[0];
rz(-2.2236373) q[0];
rz(-pi) q[1];
rz(-2.0889766) q[2];
sx q[2];
rz(-1.3498326) q[2];
sx q[2];
rz(-0.25928869) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.84089609) q[1];
sx q[1];
rz(-1.5920326) q[1];
sx q[1];
rz(-1.9663215) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21244399) q[3];
sx q[3];
rz(-1.2646958) q[3];
sx q[3];
rz(0.82992947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9156645) q[2];
sx q[2];
rz(-1.516284) q[2];
sx q[2];
rz(-0.020326745) q[2];
rz(-0.23975553) q[3];
sx q[3];
rz(-2.8862947) q[3];
sx q[3];
rz(-2.3026626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13424419) q[0];
sx q[0];
rz(-0.61779314) q[0];
sx q[0];
rz(1.0963305) q[0];
rz(1.6657375) q[1];
sx q[1];
rz(-1.219039) q[1];
sx q[1];
rz(0.79442564) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5223094) q[0];
sx q[0];
rz(-1.2566152) q[0];
sx q[0];
rz(-0.35616014) q[0];
rz(2.1486304) q[2];
sx q[2];
rz(-2.2374399) q[2];
sx q[2];
rz(1.7844328) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5880044) q[1];
sx q[1];
rz(-0.4376227) q[1];
sx q[1];
rz(2.8904084) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3268324) q[3];
sx q[3];
rz(-1.7069121) q[3];
sx q[3];
rz(0.55908113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0558001) q[2];
sx q[2];
rz(-1.5726568) q[2];
sx q[2];
rz(0.22573486) q[2];
rz(-0.24383946) q[3];
sx q[3];
rz(-0.95832458) q[3];
sx q[3];
rz(0.17914151) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8826411) q[0];
sx q[0];
rz(-1.3132341) q[0];
sx q[0];
rz(2.7146345) q[0];
rz(-0.34307617) q[1];
sx q[1];
rz(-1.1767358) q[1];
sx q[1];
rz(-0.25161904) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7386603) q[0];
sx q[0];
rz(-2.1485747) q[0];
sx q[0];
rz(-1.8389804) q[0];
rz(-pi) q[1];
rz(-0.80406739) q[2];
sx q[2];
rz(-1.40925) q[2];
sx q[2];
rz(1.3093914) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.15089825) q[1];
sx q[1];
rz(-2.7643235) q[1];
sx q[1];
rz(-2.6713085) q[1];
x q[2];
rz(-3.0915458) q[3];
sx q[3];
rz(-0.86090961) q[3];
sx q[3];
rz(1.0753701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.130927) q[2];
sx q[2];
rz(-1.6341354) q[2];
sx q[2];
rz(-2.6962386) q[2];
rz(-1.2238067) q[3];
sx q[3];
rz(-1.8755951) q[3];
sx q[3];
rz(2.7538917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39684108) q[0];
sx q[0];
rz(-0.80533177) q[0];
sx q[0];
rz(-1.1591563) q[0];
rz(1.2027488) q[1];
sx q[1];
rz(-1.1801327) q[1];
sx q[1];
rz(-0.73701) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0452249) q[0];
sx q[0];
rz(-2.8857433) q[0];
sx q[0];
rz(-2.8907828) q[0];
rz(-2.0194172) q[2];
sx q[2];
rz(-0.72690847) q[2];
sx q[2];
rz(2.2688933) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.23589489) q[1];
sx q[1];
rz(-1.7445843) q[1];
sx q[1];
rz(0.64069616) q[1];
rz(2.9748671) q[3];
sx q[3];
rz(-1.5136216) q[3];
sx q[3];
rz(-1.0852607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1272588) q[2];
sx q[2];
rz(-2.2517683) q[2];
sx q[2];
rz(-2.3298402) q[2];
rz(0.74770606) q[3];
sx q[3];
rz(-2.2816608) q[3];
sx q[3];
rz(3.0809793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49374813) q[0];
sx q[0];
rz(-2.0525377) q[0];
sx q[0];
rz(1.7895948) q[0];
rz(0.79212517) q[1];
sx q[1];
rz(-2.242576) q[1];
sx q[1];
rz(1.0101213) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2628233) q[0];
sx q[0];
rz(-2.4973625) q[0];
sx q[0];
rz(0.99272002) q[0];
x q[1];
rz(1.5158122) q[2];
sx q[2];
rz(-1.9628982) q[2];
sx q[2];
rz(-0.34811172) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6163075) q[1];
sx q[1];
rz(-1.0152363) q[1];
sx q[1];
rz(-1.6643334) q[1];
rz(-pi) q[2];
rz(2.185818) q[3];
sx q[3];
rz(-0.20032756) q[3];
sx q[3];
rz(2.1708787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.036756214) q[2];
sx q[2];
rz(-0.75883055) q[2];
sx q[2];
rz(1.7581615) q[2];
rz(-0.029959921) q[3];
sx q[3];
rz(-0.99830097) q[3];
sx q[3];
rz(-1.8647319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5112011) q[0];
sx q[0];
rz(-2.7250405) q[0];
sx q[0];
rz(-0.1314441) q[0];
rz(-0.99880544) q[1];
sx q[1];
rz(-1.9824332) q[1];
sx q[1];
rz(-1.7703895) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4192062) q[0];
sx q[0];
rz(-0.50438577) q[0];
sx q[0];
rz(0.76865159) q[0];
x q[1];
rz(-2.4553301) q[2];
sx q[2];
rz(-1.4126083) q[2];
sx q[2];
rz(-0.37294086) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.35624337) q[1];
sx q[1];
rz(-1.2977072) q[1];
sx q[1];
rz(-2.9778984) q[1];
rz(3.1375299) q[3];
sx q[3];
rz(-0.71089478) q[3];
sx q[3];
rz(1.2234185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.78352952) q[2];
sx q[2];
rz(-0.72856599) q[2];
sx q[2];
rz(-1.0895458) q[2];
rz(-0.36695668) q[3];
sx q[3];
rz(-2.7373382) q[3];
sx q[3];
rz(0.75016108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.628196) q[0];
sx q[0];
rz(-0.80444002) q[0];
sx q[0];
rz(-0.87271571) q[0];
rz(1.577042) q[1];
sx q[1];
rz(-1.4955474) q[1];
sx q[1];
rz(2.4094792) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4756389) q[0];
sx q[0];
rz(-1.7717965) q[0];
sx q[0];
rz(2.4253286) q[0];
rz(-pi) q[1];
rz(1.6953515) q[2];
sx q[2];
rz(-0.8890748) q[2];
sx q[2];
rz(2.9602697) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.32953507) q[1];
sx q[1];
rz(-1.1033711) q[1];
sx q[1];
rz(2.9653564) q[1];
rz(-pi) q[2];
rz(-1.1752582) q[3];
sx q[3];
rz(-2.8962781) q[3];
sx q[3];
rz(0.79110629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.5319891) q[2];
sx q[2];
rz(-1.0010109) q[2];
sx q[2];
rz(-0.58094376) q[2];
rz(1.1254492) q[3];
sx q[3];
rz(-1.1939129) q[3];
sx q[3];
rz(-1.1552936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4693212) q[0];
sx q[0];
rz(-0.10675616) q[0];
sx q[0];
rz(0.17053764) q[0];
rz(-1.8960309) q[1];
sx q[1];
rz(-1.5252557) q[1];
sx q[1];
rz(2.4571498) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3869707) q[0];
sx q[0];
rz(-0.59387654) q[0];
sx q[0];
rz(2.519849) q[0];
rz(-0.23875321) q[2];
sx q[2];
rz(-2.646822) q[2];
sx q[2];
rz(0.35428167) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9065793) q[1];
sx q[1];
rz(-2.7827127) q[1];
sx q[1];
rz(2.3770664) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6489469) q[3];
sx q[3];
rz(-1.6187549) q[3];
sx q[3];
rz(-0.063604442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.30274063) q[2];
sx q[2];
rz(-1.799823) q[2];
sx q[2];
rz(-1.1203241) q[2];
rz(0.66425792) q[3];
sx q[3];
rz(-0.40226007) q[3];
sx q[3];
rz(-2.6334488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6054194) q[0];
sx q[0];
rz(-0.37261951) q[0];
sx q[0];
rz(-2.5033409) q[0];
rz(3.1328746) q[1];
sx q[1];
rz(-1.1161085) q[1];
sx q[1];
rz(-0.4062103) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2531949) q[0];
sx q[0];
rz(-0.29650819) q[0];
sx q[0];
rz(-1.967708) q[0];
x q[1];
rz(-2.0044486) q[2];
sx q[2];
rz(-1.6182476) q[2];
sx q[2];
rz(2.673867) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6510424) q[1];
sx q[1];
rz(-0.79572751) q[1];
sx q[1];
rz(0.1601689) q[1];
rz(2.8088517) q[3];
sx q[3];
rz(-3.0638998) q[3];
sx q[3];
rz(-1.2487703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.093988769) q[2];
sx q[2];
rz(-1.1270019) q[2];
sx q[2];
rz(2.7809533) q[2];
rz(-2.4141198) q[3];
sx q[3];
rz(-2.45939) q[3];
sx q[3];
rz(2.8235161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8135081) q[0];
sx q[0];
rz(-2.1961975) q[0];
sx q[0];
rz(-0.16950053) q[0];
rz(-0.70973474) q[1];
sx q[1];
rz(-1.7252445) q[1];
sx q[1];
rz(-2.0555326) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5616578) q[0];
sx q[0];
rz(-0.57889639) q[0];
sx q[0];
rz(0.41335841) q[0];
rz(-1.1342088) q[2];
sx q[2];
rz(-1.2870169) q[2];
sx q[2];
rz(-1.5320317) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5063613) q[1];
sx q[1];
rz(-1.6690134) q[1];
sx q[1];
rz(2.1327095) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38137718) q[3];
sx q[3];
rz(-1.3998195) q[3];
sx q[3];
rz(-1.387763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.78709948) q[2];
sx q[2];
rz(-3.0566065) q[2];
sx q[2];
rz(-0.3698012) q[2];
rz(1.1981111) q[3];
sx q[3];
rz(-1.5912278) q[3];
sx q[3];
rz(0.91355356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13854606) q[0];
sx q[0];
rz(-1.9771165) q[0];
sx q[0];
rz(1.8856915) q[0];
rz(-2.1239602) q[1];
sx q[1];
rz(-1.3980649) q[1];
sx q[1];
rz(-1.3921888) q[1];
rz(1.5003149) q[2];
sx q[2];
rz(-1.4983836) q[2];
sx q[2];
rz(-0.56624779) q[2];
rz(-1.7012588) q[3];
sx q[3];
rz(-1.4080019) q[3];
sx q[3];
rz(-0.27346591) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
