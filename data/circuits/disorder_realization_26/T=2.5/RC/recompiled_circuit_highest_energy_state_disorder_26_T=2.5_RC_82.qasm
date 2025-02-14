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
rz(0.8910203) q[0];
sx q[0];
rz(1.855259) q[0];
sx q[0];
rz(8.5387908) q[0];
rz(0.52357829) q[1];
sx q[1];
rz(3.7012586) q[1];
sx q[1];
rz(8.1044365) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1170066) q[0];
sx q[0];
rz(-1.3947233) q[0];
sx q[0];
rz(2.1137289) q[0];
rz(-pi) q[1];
rz(-2.7971326) q[2];
sx q[2];
rz(-2.3340324) q[2];
sx q[2];
rz(-1.5429614) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7976171) q[1];
sx q[1];
rz(-1.5098176) q[1];
sx q[1];
rz(1.0073376) q[1];
x q[2];
rz(-2.376972) q[3];
sx q[3];
rz(-0.9249827) q[3];
sx q[3];
rz(0.82987204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3628799) q[2];
sx q[2];
rz(-1.3211297) q[2];
sx q[2];
rz(2.2134181) q[2];
rz(1.5859531) q[3];
sx q[3];
rz(-1.0471683) q[3];
sx q[3];
rz(1.9716523) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8856119) q[0];
sx q[0];
rz(-0.33807999) q[0];
sx q[0];
rz(2.0804491) q[0];
rz(-1.2127009) q[1];
sx q[1];
rz(-0.78293982) q[1];
sx q[1];
rz(1.3139668) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7018258) q[0];
sx q[0];
rz(-1.9226388) q[0];
sx q[0];
rz(0.74811305) q[0];
rz(-0.92483123) q[2];
sx q[2];
rz(-1.926117) q[2];
sx q[2];
rz(0.79910797) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3281341) q[1];
sx q[1];
rz(-2.5812006) q[1];
sx q[1];
rz(1.8916393) q[1];
rz(2.2960194) q[3];
sx q[3];
rz(-1.9872287) q[3];
sx q[3];
rz(1.4786947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5709915) q[2];
sx q[2];
rz(-0.97390318) q[2];
sx q[2];
rz(0.904733) q[2];
rz(1.5500801) q[3];
sx q[3];
rz(-1.855987) q[3];
sx q[3];
rz(1.2929644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0631113) q[0];
sx q[0];
rz(-3.0516629) q[0];
sx q[0];
rz(-0.91823804) q[0];
rz(2.4507484) q[1];
sx q[1];
rz(-1.7450688) q[1];
sx q[1];
rz(2.3450559) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18609234) q[0];
sx q[0];
rz(-1.4791795) q[0];
sx q[0];
rz(-1.7797526) q[0];
rz(-pi) q[1];
rz(-2.1941625) q[2];
sx q[2];
rz(-3.0683225) q[2];
sx q[2];
rz(0.87033349) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5480876) q[1];
sx q[1];
rz(-1.6789762) q[1];
sx q[1];
rz(-0.87203474) q[1];
rz(-pi) q[2];
rz(-0.34557202) q[3];
sx q[3];
rz(-2.1949147) q[3];
sx q[3];
rz(2.6081599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.010926509) q[2];
sx q[2];
rz(-2.1911759) q[2];
sx q[2];
rz(-1.2476904) q[2];
rz(0.55825663) q[3];
sx q[3];
rz(-1.6853251) q[3];
sx q[3];
rz(-0.021350967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70304865) q[0];
sx q[0];
rz(-3.092364) q[0];
sx q[0];
rz(2.4141648) q[0];
rz(-0.1618596) q[1];
sx q[1];
rz(-1.9700123) q[1];
sx q[1];
rz(1.9836099) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25140554) q[0];
sx q[0];
rz(-1.5047538) q[0];
sx q[0];
rz(-1.4736045) q[0];
x q[1];
rz(-1.5776872) q[2];
sx q[2];
rz(-3.0022923) q[2];
sx q[2];
rz(1.3819288) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.53044415) q[1];
sx q[1];
rz(-1.6940677) q[1];
sx q[1];
rz(1.7632496) q[1];
rz(-pi) q[2];
rz(-1.4876266) q[3];
sx q[3];
rz(-1.0462049) q[3];
sx q[3];
rz(2.7980141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8689279) q[2];
sx q[2];
rz(-1.2020035) q[2];
sx q[2];
rz(1.391601) q[2];
rz(0.23022716) q[3];
sx q[3];
rz(-1.9672491) q[3];
sx q[3];
rz(-2.9496884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72628438) q[0];
sx q[0];
rz(-3.0400161) q[0];
sx q[0];
rz(2.0368077) q[0];
rz(0.11100189) q[1];
sx q[1];
rz(-1.1155201) q[1];
sx q[1];
rz(1.8288137) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3462145) q[0];
sx q[0];
rz(-1.7722478) q[0];
sx q[0];
rz(0.18152118) q[0];
rz(-pi) q[1];
rz(3.0022125) q[2];
sx q[2];
rz(-0.5682033) q[2];
sx q[2];
rz(0.85124337) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1089563) q[1];
sx q[1];
rz(-1.3707038) q[1];
sx q[1];
rz(1.1632827) q[1];
rz(-0.41683414) q[3];
sx q[3];
rz(-0.89589828) q[3];
sx q[3];
rz(-1.0639497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.50292) q[2];
sx q[2];
rz(-2.0827677) q[2];
sx q[2];
rz(2.1144833) q[2];
rz(2.1549759) q[3];
sx q[3];
rz(-0.28237453) q[3];
sx q[3];
rz(1.4237039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8091549) q[0];
sx q[0];
rz(-1.2323392) q[0];
sx q[0];
rz(-2.3409081) q[0];
rz(-2.370131) q[1];
sx q[1];
rz(-1.0779251) q[1];
sx q[1];
rz(1.3763743) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9405917) q[0];
sx q[0];
rz(-1.1744813) q[0];
sx q[0];
rz(-0.70756377) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1065953) q[2];
sx q[2];
rz(-2.0366324) q[2];
sx q[2];
rz(-0.83110561) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.034198) q[1];
sx q[1];
rz(-1.5840816) q[1];
sx q[1];
rz(2.4360869) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8326859) q[3];
sx q[3];
rz(-2.247274) q[3];
sx q[3];
rz(0.075084075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.8419522) q[2];
sx q[2];
rz(-1.2066634) q[2];
sx q[2];
rz(-0.017814962) q[2];
rz(2.8282015) q[3];
sx q[3];
rz(-2.119795) q[3];
sx q[3];
rz(0.71948403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3571091) q[0];
sx q[0];
rz(-1.5743558) q[0];
sx q[0];
rz(-2.9743279) q[0];
rz(2.3978865) q[1];
sx q[1];
rz(-0.88631648) q[1];
sx q[1];
rz(3.0053265) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3534451) q[0];
sx q[0];
rz(-2.2290725) q[0];
sx q[0];
rz(2.7892755) q[0];
rz(-pi) q[1];
rz(2.0799147) q[2];
sx q[2];
rz(-0.3442685) q[2];
sx q[2];
rz(-1.4835675) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4748342) q[1];
sx q[1];
rz(-1.6377565) q[1];
sx q[1];
rz(0.35891846) q[1];
x q[2];
rz(1.5149917) q[3];
sx q[3];
rz(-1.9161092) q[3];
sx q[3];
rz(-1.2962411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.83884376) q[2];
sx q[2];
rz(-0.15041298) q[2];
sx q[2];
rz(2.0602843) q[2];
rz(0.14351621) q[3];
sx q[3];
rz(-1.84294) q[3];
sx q[3];
rz(1.4474086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3926587) q[0];
sx q[0];
rz(-1.8086139) q[0];
sx q[0];
rz(2.4984711) q[0];
rz(2.2195623) q[1];
sx q[1];
rz(-2.8755201) q[1];
sx q[1];
rz(1.6304852) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63207549) q[0];
sx q[0];
rz(-2.7537573) q[0];
sx q[0];
rz(1.0799506) q[0];
rz(-pi) q[1];
x q[1];
rz(1.950932) q[2];
sx q[2];
rz(-0.30001773) q[2];
sx q[2];
rz(-2.0643108) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0865759) q[1];
sx q[1];
rz(-2.9447703) q[1];
sx q[1];
rz(1.6114525) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.891252) q[3];
sx q[3];
rz(-2.0600187) q[3];
sx q[3];
rz(-1.7699522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1954605) q[2];
sx q[2];
rz(-1.5617424) q[2];
sx q[2];
rz(1.3004318) q[2];
rz(-2.2293633) q[3];
sx q[3];
rz(-1.2765086) q[3];
sx q[3];
rz(0.31360489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6533971) q[0];
sx q[0];
rz(-0.47976872) q[0];
sx q[0];
rz(-2.524014) q[0];
rz(0.081534475) q[1];
sx q[1];
rz(-2.640994) q[1];
sx q[1];
rz(-0.68731442) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9195557) q[0];
sx q[0];
rz(-1.3341781) q[0];
sx q[0];
rz(3.0927004) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26726802) q[2];
sx q[2];
rz(-1.4630001) q[2];
sx q[2];
rz(-1.9330618) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7973104) q[1];
sx q[1];
rz(-2.629529) q[1];
sx q[1];
rz(-0.6060098) q[1];
x q[2];
rz(-3.0932759) q[3];
sx q[3];
rz(-1.028065) q[3];
sx q[3];
rz(2.2996772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0240137) q[2];
sx q[2];
rz(-1.1774096) q[2];
sx q[2];
rz(3.0702316) q[2];
rz(-1.2355545) q[3];
sx q[3];
rz(-0.33696285) q[3];
sx q[3];
rz(2.5753042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5679034) q[0];
sx q[0];
rz(-2.8901143) q[0];
sx q[0];
rz(-3.0036744) q[0];
rz(0.57016405) q[1];
sx q[1];
rz(-2.0591996) q[1];
sx q[1];
rz(-0.015448419) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.041067657) q[0];
sx q[0];
rz(-2.3599632) q[0];
sx q[0];
rz(2.0360721) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1006946) q[2];
sx q[2];
rz(-2.2743871) q[2];
sx q[2];
rz(-0.64794618) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5970832) q[1];
sx q[1];
rz(-0.5209777) q[1];
sx q[1];
rz(-2.8057666) q[1];
x q[2];
rz(1.3872765) q[3];
sx q[3];
rz(-1.5158815) q[3];
sx q[3];
rz(-2.4650989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40921679) q[2];
sx q[2];
rz(-2.3541383) q[2];
sx q[2];
rz(2.709205) q[2];
rz(-2.2435097) q[3];
sx q[3];
rz(-0.87528527) q[3];
sx q[3];
rz(-1.5842277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.4846004) q[0];
sx q[0];
rz(-2.0130172) q[0];
sx q[0];
rz(0.76242557) q[0];
rz(0.55105974) q[1];
sx q[1];
rz(-1.8012128) q[1];
sx q[1];
rz(2.6691379) q[1];
rz(-1.6161935) q[2];
sx q[2];
rz(-1.207796) q[2];
sx q[2];
rz(2.7784204) q[2];
rz(-0.9394218) q[3];
sx q[3];
rz(-0.71790725) q[3];
sx q[3];
rz(0.090286615) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
