OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6686749) q[0];
sx q[0];
rz(-0.023107419) q[0];
sx q[0];
rz(0.90149108) q[0];
rz(1.3540406) q[1];
sx q[1];
rz(4.7572264) q[1];
sx q[1];
rz(4.9487257) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5164099) q[0];
sx q[0];
rz(-1.7889272) q[0];
sx q[0];
rz(-1.6263917) q[0];
x q[1];
rz(-0.15057474) q[2];
sx q[2];
rz(-1.625744) q[2];
sx q[2];
rz(-0.51412941) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1062231) q[1];
sx q[1];
rz(-1.6236042) q[1];
sx q[1];
rz(-1.04671) q[1];
rz(0.24231916) q[3];
sx q[3];
rz(-1.524914) q[3];
sx q[3];
rz(2.322779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2287075) q[2];
sx q[2];
rz(-0.097147377) q[2];
sx q[2];
rz(2.482282) q[2];
rz(-0.77624503) q[3];
sx q[3];
rz(-3.1228437) q[3];
sx q[3];
rz(-0.68500486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.2289675) q[0];
sx q[0];
rz(-1.9321059) q[0];
sx q[0];
rz(1.9483161) q[0];
rz(-0.044366447) q[1];
sx q[1];
rz(-0.012244789) q[1];
sx q[1];
rz(2.9150229) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8346342) q[0];
sx q[0];
rz(-0.0037841664) q[0];
sx q[0];
rz(-1.1192805) q[0];
rz(2.7625257) q[2];
sx q[2];
rz(-1.5866536) q[2];
sx q[2];
rz(0.011034688) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5331796) q[1];
sx q[1];
rz(-2.7044109) q[1];
sx q[1];
rz(-2.9730148) q[1];
rz(-pi) q[2];
rz(-1.8227406) q[3];
sx q[3];
rz(-2.3980707) q[3];
sx q[3];
rz(-1.8094667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.442753) q[2];
sx q[2];
rz(-1.5696462) q[2];
sx q[2];
rz(-1.6119831) q[2];
rz(-2.2085341) q[3];
sx q[3];
rz(-1.482684) q[3];
sx q[3];
rz(-0.23811594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9642692) q[0];
sx q[0];
rz(-0.015559109) q[0];
sx q[0];
rz(2.9413057) q[0];
rz(3.1411723) q[1];
sx q[1];
rz(-0.93685189) q[1];
sx q[1];
rz(-0.013484152) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5258085) q[0];
sx q[0];
rz(-1.1466525) q[0];
sx q[0];
rz(3.0706617) q[0];
x q[1];
rz(2.5691367) q[2];
sx q[2];
rz(-3.0623263) q[2];
sx q[2];
rz(-2.5545189) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.84220195) q[1];
sx q[1];
rz(-2.9889571) q[1];
sx q[1];
rz(2.045689) q[1];
rz(-pi) q[2];
rz(1.5192658) q[3];
sx q[3];
rz(-1.4084219) q[3];
sx q[3];
rz(-2.3414617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1991594) q[2];
sx q[2];
rz(-1.556267) q[2];
sx q[2];
rz(-1.5996492) q[2];
rz(2.13983) q[3];
sx q[3];
rz(-0.21654138) q[3];
sx q[3];
rz(-2.969363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37429419) q[0];
sx q[0];
rz(-2.9758487) q[0];
sx q[0];
rz(0.34749183) q[0];
rz(-0.63684288) q[1];
sx q[1];
rz(-3.1359735) q[1];
sx q[1];
rz(1.9510795) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4012642) q[0];
sx q[0];
rz(-1.51121) q[0];
sx q[0];
rz(-1.7028832) q[0];
rz(-pi) q[1];
rz(-1.6830446) q[2];
sx q[2];
rz(-3.0722716) q[2];
sx q[2];
rz(1.6905418) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.61100436) q[1];
sx q[1];
rz(-1.3827033) q[1];
sx q[1];
rz(-0.76716073) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7905964) q[3];
sx q[3];
rz(-1.2145208) q[3];
sx q[3];
rz(1.364546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5845329) q[2];
sx q[2];
rz(-3.1019042) q[2];
sx q[2];
rz(-1.3603127) q[2];
rz(-1.6654061) q[3];
sx q[3];
rz(-1.5675661) q[3];
sx q[3];
rz(0.40569693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0535468) q[0];
sx q[0];
rz(-2.4021554) q[0];
sx q[0];
rz(-0.0060225688) q[0];
rz(-1.7209523) q[1];
sx q[1];
rz(-0.050844897) q[1];
sx q[1];
rz(-3.0555225) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6688441) q[0];
sx q[0];
rz(-1.5875582) q[0];
sx q[0];
rz(-3.0344323) q[0];
rz(1.5447247) q[2];
sx q[2];
rz(-1.6258473) q[2];
sx q[2];
rz(-1.2616273) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4115178) q[1];
sx q[1];
rz(-3.0429384) q[1];
sx q[1];
rz(-2.7559571) q[1];
x q[2];
rz(-0.045343355) q[3];
sx q[3];
rz(-1.6312977) q[3];
sx q[3];
rz(-2.7013472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.033279557) q[2];
sx q[2];
rz(-0.72282183) q[2];
sx q[2];
rz(-1.7576199) q[2];
rz(2.7464187) q[3];
sx q[3];
rz(-3.0925909) q[3];
sx q[3];
rz(1.9743732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5351717) q[0];
sx q[0];
rz(-0.21927729) q[0];
sx q[0];
rz(-1.0004591) q[0];
rz(-0.74042997) q[1];
sx q[1];
rz(-2.705997) q[1];
sx q[1];
rz(0.37954095) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.767547) q[0];
sx q[0];
rz(-1.5963608) q[0];
sx q[0];
rz(0.095439571) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8431115) q[2];
sx q[2];
rz(-2.0305995) q[2];
sx q[2];
rz(1.0155755) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4011456) q[1];
sx q[1];
rz(-1.5441431) q[1];
sx q[1];
rz(2.074764) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0197755) q[3];
sx q[3];
rz(-1.4662577) q[3];
sx q[3];
rz(3.0825101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.015811054) q[2];
sx q[2];
rz(-2.8534079) q[2];
sx q[2];
rz(0.077032653) q[2];
rz(0.020126255) q[3];
sx q[3];
rz(-3.088981) q[3];
sx q[3];
rz(-2.3441815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.035148419) q[0];
sx q[0];
rz(-0.012520944) q[0];
sx q[0];
rz(1.7092108) q[0];
rz(-2.8445981) q[1];
sx q[1];
rz(-0.15428267) q[1];
sx q[1];
rz(-3.0800379) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2765574) q[0];
sx q[0];
rz(-1.1662959) q[0];
sx q[0];
rz(0.96488953) q[0];
rz(-pi) q[1];
rz(-0.02587895) q[2];
sx q[2];
rz(-0.21458438) q[2];
sx q[2];
rz(-1.1550922) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0377215) q[1];
sx q[1];
rz(-1.6391409) q[1];
sx q[1];
rz(-2.8496024) q[1];
rz(-pi) q[2];
rz(2.706884) q[3];
sx q[3];
rz(-1.604933) q[3];
sx q[3];
rz(-3.1354836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.84294549) q[2];
sx q[2];
rz(-0.25235287) q[2];
sx q[2];
rz(1.7227777) q[2];
rz(1.8055387) q[3];
sx q[3];
rz(-0.027848363) q[3];
sx q[3];
rz(1.4068039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1143188) q[0];
sx q[0];
rz(-2.424746) q[0];
sx q[0];
rz(1.640821) q[0];
rz(-2.0188792) q[1];
sx q[1];
rz(-0.32674679) q[1];
sx q[1];
rz(1.2935125) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0121078) q[0];
sx q[0];
rz(-2.2177758) q[0];
sx q[0];
rz(2.7478474) q[0];
rz(-pi) q[1];
rz(0.40385623) q[2];
sx q[2];
rz(-2.4060898) q[2];
sx q[2];
rz(-2.1063358) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.499524) q[1];
sx q[1];
rz(-1.4022938) q[1];
sx q[1];
rz(-2.8649855) q[1];
x q[2];
rz(0.82610749) q[3];
sx q[3];
rz(-0.87658054) q[3];
sx q[3];
rz(0.40761596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5844172) q[2];
sx q[2];
rz(-0.36089218) q[2];
sx q[2];
rz(-1.7822251) q[2];
rz(2.6946097) q[3];
sx q[3];
rz(-0.042526571) q[3];
sx q[3];
rz(2.9744448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3249224) q[0];
sx q[0];
rz(-1.8055547) q[0];
sx q[0];
rz(1.100612) q[0];
rz(-1.0758411) q[1];
sx q[1];
rz(-0.64972076) q[1];
sx q[1];
rz(-2.4646387) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9788584) q[0];
sx q[0];
rz(-1.7679201) q[0];
sx q[0];
rz(-2.7480024) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4482037) q[2];
sx q[2];
rz(-1.685678) q[2];
sx q[2];
rz(-0.5960532) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9857603) q[1];
sx q[1];
rz(-1.5712156) q[1];
sx q[1];
rz(-3.1414933) q[1];
rz(-1.7221801) q[3];
sx q[3];
rz(-1.1202235) q[3];
sx q[3];
rz(-2.7344658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1415652) q[2];
sx q[2];
rz(-0.0024777369) q[2];
sx q[2];
rz(2.4139717) q[2];
rz(-1.242312) q[3];
sx q[3];
rz(-0.036402313) q[3];
sx q[3];
rz(-1.1718933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90945554) q[0];
sx q[0];
rz(-1.0219034) q[0];
sx q[0];
rz(2.452028) q[0];
rz(-1.6652971) q[1];
sx q[1];
rz(-0.25054014) q[1];
sx q[1];
rz(-2.9857059) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7321701) q[0];
sx q[0];
rz(-0.66930938) q[0];
sx q[0];
rz(-0.21440345) q[0];
rz(-pi) q[1];
x q[1];
rz(0.25145289) q[2];
sx q[2];
rz(-1.6864459) q[2];
sx q[2];
rz(-0.72249362) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0956679) q[1];
sx q[1];
rz(-1.5692488) q[1];
sx q[1];
rz(-3.1374088) q[1];
rz(1.4149551) q[3];
sx q[3];
rz(-1.4771802) q[3];
sx q[3];
rz(-2.7839212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.186824) q[2];
sx q[2];
rz(-0.12207741) q[2];
sx q[2];
rz(-0.99511498) q[2];
rz(-2.9156445) q[3];
sx q[3];
rz(-0.048361691) q[3];
sx q[3];
rz(2.3154955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7786998) q[0];
sx q[0];
rz(-0.97465546) q[0];
sx q[0];
rz(1.4018651) q[0];
rz(1.4317935) q[1];
sx q[1];
rz(-1.2964389) q[1];
sx q[1];
rz(-2.5242205) q[1];
rz(-2.4564248) q[2];
sx q[2];
rz(-1.6772267) q[2];
sx q[2];
rz(-0.75135372) q[2];
rz(0.13629436) q[3];
sx q[3];
rz(-1.1850428) q[3];
sx q[3];
rz(-0.94120126) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
