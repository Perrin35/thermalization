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
rz(1.6549702) q[1];
sx q[1];
rz(4.4284664) q[1];
sx q[1];
rz(9.5944302) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9823224) q[0];
sx q[0];
rz(-1.9290315) q[0];
sx q[0];
rz(2.0824177) q[0];
rz(-pi) q[1];
rz(1.9965788) q[2];
sx q[2];
rz(-2.5822431) q[2];
sx q[2];
rz(0.9445136) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4028266) q[1];
sx q[1];
rz(-1.1753653) q[1];
sx q[1];
rz(-3.1185802) q[1];
rz(0.98244169) q[3];
sx q[3];
rz(-0.37068493) q[3];
sx q[3];
rz(1.6904168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.22592813) q[2];
sx q[2];
rz(-1.6253086) q[2];
sx q[2];
rz(0.020326745) q[2];
rz(-2.9018371) q[3];
sx q[3];
rz(-0.25529796) q[3];
sx q[3];
rz(0.83893004) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0073485) q[0];
sx q[0];
rz(-2.5237995) q[0];
sx q[0];
rz(1.0963305) q[0];
rz(1.4758551) q[1];
sx q[1];
rz(-1.219039) q[1];
sx q[1];
rz(2.347167) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3045469) q[0];
sx q[0];
rz(-1.908806) q[0];
sx q[0];
rz(-1.2370626) q[0];
rz(-pi) q[1];
rz(0.6068318) q[2];
sx q[2];
rz(-2.2893659) q[2];
sx q[2];
rz(-2.5963714) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2455755) q[1];
sx q[1];
rz(-1.676325) q[1];
sx q[1];
rz(0.42550931) q[1];
x q[2];
rz(1.7678472) q[3];
sx q[3];
rz(-0.76579491) q[3];
sx q[3];
rz(1.1546419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0558001) q[2];
sx q[2];
rz(-1.5689359) q[2];
sx q[2];
rz(2.9158578) q[2];
rz(0.24383946) q[3];
sx q[3];
rz(-0.95832458) q[3];
sx q[3];
rz(-0.17914151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
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
rz(-1.8826411) q[0];
sx q[0];
rz(-1.8283586) q[0];
sx q[0];
rz(-0.42695811) q[0];
rz(0.34307617) q[1];
sx q[1];
rz(-1.9648569) q[1];
sx q[1];
rz(-0.25161904) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7386603) q[0];
sx q[0];
rz(-2.1485747) q[0];
sx q[0];
rz(1.3026122) q[0];
rz(0.80406739) q[2];
sx q[2];
rz(-1.40925) q[2];
sx q[2];
rz(-1.3093914) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15089825) q[1];
sx q[1];
rz(-0.37726918) q[1];
sx q[1];
rz(-2.6713085) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0915458) q[3];
sx q[3];
rz(-0.86090961) q[3];
sx q[3];
rz(2.0662226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0106657) q[2];
sx q[2];
rz(-1.5074573) q[2];
sx q[2];
rz(-0.44535401) q[2];
rz(-1.917786) q[3];
sx q[3];
rz(-1.8755951) q[3];
sx q[3];
rz(0.38770097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39684108) q[0];
sx q[0];
rz(-0.80533177) q[0];
sx q[0];
rz(-1.1591563) q[0];
rz(1.9388439) q[1];
sx q[1];
rz(-1.9614599) q[1];
sx q[1];
rz(-0.73701) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16251462) q[0];
sx q[0];
rz(-1.8184699) q[0];
sx q[0];
rz(-1.635627) q[0];
rz(-pi) q[1];
rz(1.1221755) q[2];
sx q[2];
rz(-0.72690847) q[2];
sx q[2];
rz(2.2688933) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1070246) q[1];
sx q[1];
rz(-0.66063297) q[1];
sx q[1];
rz(-0.28566499) q[1];
x q[2];
rz(-1.5128194) q[3];
sx q[3];
rz(-1.4043458) q[3];
sx q[3];
rz(-2.6464407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1272588) q[2];
sx q[2];
rz(-0.88982439) q[2];
sx q[2];
rz(0.8117525) q[2];
rz(-2.3938866) q[3];
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
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6478445) q[0];
sx q[0];
rz(-1.0890549) q[0];
sx q[0];
rz(1.3519979) q[0];
rz(-2.3494675) q[1];
sx q[1];
rz(-2.242576) q[1];
sx q[1];
rz(1.0101213) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3528004) q[0];
sx q[0];
rz(-1.236434) q[0];
sx q[0];
rz(1.0092606) q[0];
rz(-pi) q[1];
rz(0.3926362) q[2];
sx q[2];
rz(-1.5199888) q[2];
sx q[2];
rz(-1.2016553) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7923342) q[1];
sx q[1];
rz(-2.5790303) q[1];
sx q[1];
rz(2.9922561) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11662401) q[3];
sx q[3];
rz(-1.734048) q[3];
sx q[3];
rz(1.5463055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.036756214) q[2];
sx q[2];
rz(-2.3827621) q[2];
sx q[2];
rz(-1.7581615) q[2];
rz(3.1116327) q[3];
sx q[3];
rz(-0.99830097) q[3];
sx q[3];
rz(1.2768607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6303915) q[0];
sx q[0];
rz(-2.7250405) q[0];
sx q[0];
rz(-0.1314441) q[0];
rz(0.99880544) q[1];
sx q[1];
rz(-1.1591594) q[1];
sx q[1];
rz(-1.7703895) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11267647) q[0];
sx q[0];
rz(-1.2160016) q[0];
sx q[0];
rz(1.2043857) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68626257) q[2];
sx q[2];
rz(-1.7289844) q[2];
sx q[2];
rz(-0.37294086) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2590717) q[1];
sx q[1];
rz(-1.4132199) q[1];
sx q[1];
rz(-1.8474008) q[1];
rz(-pi) q[2];
rz(-1.5672979) q[3];
sx q[3];
rz(-0.85990866) q[3];
sx q[3];
rz(1.2287799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.78352952) q[2];
sx q[2];
rz(-0.72856599) q[2];
sx q[2];
rz(-1.0895458) q[2];
rz(2.774636) q[3];
sx q[3];
rz(-2.7373382) q[3];
sx q[3];
rz(-2.3914316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.628196) q[0];
sx q[0];
rz(-2.3371526) q[0];
sx q[0];
rz(-0.87271571) q[0];
rz(-1.5645507) q[1];
sx q[1];
rz(-1.4955474) q[1];
sx q[1];
rz(-0.73211342) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66595378) q[0];
sx q[0];
rz(-1.7717965) q[0];
sx q[0];
rz(-0.71626407) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68552981) q[2];
sx q[2];
rz(-1.6674124) q[2];
sx q[2];
rz(-1.8308507) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8120576) q[1];
sx q[1];
rz(-1.1033711) q[1];
sx q[1];
rz(-2.9653564) q[1];
rz(-pi) q[2];
rz(-1.3437531) q[3];
sx q[3];
rz(-1.6645091) q[3];
sx q[3];
rz(-0.39484398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6096036) q[2];
sx q[2];
rz(-2.1405818) q[2];
sx q[2];
rz(-2.5606489) q[2];
rz(-2.0161435) q[3];
sx q[3];
rz(-1.1939129) q[3];
sx q[3];
rz(-1.1552936) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4693212) q[0];
sx q[0];
rz(-3.0348365) q[0];
sx q[0];
rz(2.971055) q[0];
rz(-1.8960309) q[1];
sx q[1];
rz(-1.6163369) q[1];
sx q[1];
rz(-2.4571498) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3869707) q[0];
sx q[0];
rz(-2.5477161) q[0];
sx q[0];
rz(0.62174364) q[0];
x q[1];
rz(-2.9028394) q[2];
sx q[2];
rz(-0.49477067) q[2];
sx q[2];
rz(0.35428167) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.43757968) q[1];
sx q[1];
rz(-1.3145169) q[1];
sx q[1];
rz(1.3167472) q[1];
rz(-pi) q[2];
rz(1.6309647) q[3];
sx q[3];
rz(-2.2188713) q[3];
sx q[3];
rz(-1.5435406) q[3];
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
rz(2.0212685) q[2];
rz(-0.66425792) q[3];
sx q[3];
rz(-0.40226007) q[3];
sx q[3];
rz(-0.50814381) q[3];
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
rz(-2.6054194) q[0];
sx q[0];
rz(-0.37261951) q[0];
sx q[0];
rz(0.63825178) q[0];
rz(-3.1328746) q[1];
sx q[1];
rz(-1.1161085) q[1];
sx q[1];
rz(0.4062103) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063651872) q[0];
sx q[0];
rz(-1.6839875) q[0];
sx q[0];
rz(-1.8454396) q[0];
rz(-2.0044486) q[2];
sx q[2];
rz(-1.523345) q[2];
sx q[2];
rz(0.46772568) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.26365852) q[1];
sx q[1];
rz(-2.3535427) q[1];
sx q[1];
rz(1.4093983) q[1];
x q[2];
rz(-1.5962192) q[3];
sx q[3];
rz(-1.4973728) q[3];
sx q[3];
rz(1.5591476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.093988769) q[2];
sx q[2];
rz(-2.0145907) q[2];
sx q[2];
rz(-2.7809533) q[2];
rz(-2.4141198) q[3];
sx q[3];
rz(-2.45939) q[3];
sx q[3];
rz(-0.31807652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3280846) q[0];
sx q[0];
rz(-2.1961975) q[0];
sx q[0];
rz(2.9720921) q[0];
rz(-0.70973474) q[1];
sx q[1];
rz(-1.4163481) q[1];
sx q[1];
rz(-1.08606) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0789886) q[0];
sx q[0];
rz(-2.0955756) q[0];
sx q[0];
rz(1.3140509) q[0];
rz(-pi) q[1];
rz(-2.8302221) q[2];
sx q[2];
rz(-1.9888005) q[2];
sx q[2];
rz(3.0504464) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5063613) q[1];
sx q[1];
rz(-1.6690134) q[1];
sx q[1];
rz(1.0088831) q[1];
rz(1.754722) q[3];
sx q[3];
rz(-1.1952595) q[3];
sx q[3];
rz(0.25115764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3544932) q[2];
sx q[2];
rz(-0.084986173) q[2];
sx q[2];
rz(-0.3698012) q[2];
rz(-1.9434816) q[3];
sx q[3];
rz(-1.5912278) q[3];
sx q[3];
rz(0.91355356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13854606) q[0];
sx q[0];
rz(-1.9771165) q[0];
sx q[0];
rz(1.8856915) q[0];
rz(2.1239602) q[1];
sx q[1];
rz(-1.7435278) q[1];
sx q[1];
rz(1.7494038) q[1];
rz(1.5003149) q[2];
sx q[2];
rz(-1.4983836) q[2];
sx q[2];
rz(-0.56624779) q[2];
rz(-2.9774278) q[3];
sx q[3];
rz(-1.4420684) q[3];
sx q[3];
rz(1.3185929) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
