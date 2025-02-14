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
rz(0.34613553) q[0];
sx q[0];
rz(-1.5032285) q[0];
sx q[0];
rz(0.7315973) q[0];
rz(-2.6919964) q[1];
sx q[1];
rz(-0.1875339) q[1];
sx q[1];
rz(-0.43391689) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2303378) q[0];
sx q[0];
rz(-0.83090913) q[0];
sx q[0];
rz(0.047538443) q[0];
rz(-pi) q[1];
rz(-2.163274) q[2];
sx q[2];
rz(-1.0772395) q[2];
sx q[2];
rz(1.6101642) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0632532) q[1];
sx q[1];
rz(-1.1800449) q[1];
sx q[1];
rz(-0.88576742) q[1];
x q[2];
rz(-1.3796666) q[3];
sx q[3];
rz(-1.5401296) q[3];
sx q[3];
rz(-0.14740869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4291541) q[2];
sx q[2];
rz(-1.6753847) q[2];
sx q[2];
rz(2.1377371) q[2];
rz(-2.6058274) q[3];
sx q[3];
rz(-2.2771213) q[3];
sx q[3];
rz(-0.0011477688) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67343229) q[0];
sx q[0];
rz(-2.4591481) q[0];
sx q[0];
rz(-2.4144507) q[0];
rz(-0.06761059) q[1];
sx q[1];
rz(-1.3550242) q[1];
sx q[1];
rz(2.0179857) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.362181) q[0];
sx q[0];
rz(-1.0955278) q[0];
sx q[0];
rz(-1.2527466) q[0];
rz(-0.065326377) q[2];
sx q[2];
rz(-1.4484754) q[2];
sx q[2];
rz(2.9768012) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1084845) q[1];
sx q[1];
rz(-0.29769167) q[1];
sx q[1];
rz(-1.2254524) q[1];
rz(-pi) q[2];
rz(-0.83060925) q[3];
sx q[3];
rz(-0.6081444) q[3];
sx q[3];
rz(2.6839395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2850538) q[2];
sx q[2];
rz(-1.5736138) q[2];
sx q[2];
rz(0.48173586) q[2];
rz(0.5528062) q[3];
sx q[3];
rz(-2.063664) q[3];
sx q[3];
rz(0.089574561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8168617) q[0];
sx q[0];
rz(-0.4751927) q[0];
sx q[0];
rz(0.09356308) q[0];
rz(1.3803253) q[1];
sx q[1];
rz(-2.1880136) q[1];
sx q[1];
rz(0.63724744) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4555511) q[0];
sx q[0];
rz(-1.0555869) q[0];
sx q[0];
rz(2.0705219) q[0];
rz(-pi) q[1];
rz(0.60814823) q[2];
sx q[2];
rz(-1.7575193) q[2];
sx q[2];
rz(3.1352459) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2356469) q[1];
sx q[1];
rz(-1.596758) q[1];
sx q[1];
rz(1.9683377) q[1];
x q[2];
rz(-1.6138541) q[3];
sx q[3];
rz(-1.8612487) q[3];
sx q[3];
rz(1.2424948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9351585) q[2];
sx q[2];
rz(-0.55525246) q[2];
sx q[2];
rz(-0.68378249) q[2];
rz(-0.23009662) q[3];
sx q[3];
rz(-1.7347387) q[3];
sx q[3];
rz(-0.42207119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049659599) q[0];
sx q[0];
rz(-2.6347418) q[0];
sx q[0];
rz(-1.7012713) q[0];
rz(-2.6662042) q[1];
sx q[1];
rz(-2.5071867) q[1];
sx q[1];
rz(-2.5604274) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0770532) q[0];
sx q[0];
rz(-2.4145899) q[0];
sx q[0];
rz(-1.6995656) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23431205) q[2];
sx q[2];
rz(-2.0221124) q[2];
sx q[2];
rz(-2.6036604) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.95422259) q[1];
sx q[1];
rz(-2.0285602) q[1];
sx q[1];
rz(-1.2332031) q[1];
x q[2];
rz(-1.4122333) q[3];
sx q[3];
rz(-2.057029) q[3];
sx q[3];
rz(-2.9610046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0541957) q[2];
sx q[2];
rz(-0.18614686) q[2];
sx q[2];
rz(-2.6206214) q[2];
rz(-1.6446796) q[3];
sx q[3];
rz(-1.4119586) q[3];
sx q[3];
rz(1.6269256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(1.5882551) q[0];
sx q[0];
rz(-0.27565685) q[0];
sx q[0];
rz(2.0113373) q[0];
rz(-2.0626119) q[1];
sx q[1];
rz(-1.9422453) q[1];
sx q[1];
rz(-2.030453) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95100194) q[0];
sx q[0];
rz(-1.7365121) q[0];
sx q[0];
rz(2.034305) q[0];
x q[1];
rz(-3.0658998) q[2];
sx q[2];
rz(-2.6106129) q[2];
sx q[2];
rz(-0.39226433) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9359253) q[1];
sx q[1];
rz(-1.4135409) q[1];
sx q[1];
rz(2.2155511) q[1];
x q[2];
rz(-2.2258198) q[3];
sx q[3];
rz(-2.0213599) q[3];
sx q[3];
rz(-2.5289867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.40372103) q[2];
sx q[2];
rz(-2.626494) q[2];
sx q[2];
rz(-2.1007288) q[2];
rz(-2.3216085) q[3];
sx q[3];
rz(-1.2330202) q[3];
sx q[3];
rz(-0.6238873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2417004) q[0];
sx q[0];
rz(-2.5427759) q[0];
sx q[0];
rz(2.4851121) q[0];
rz(1.7680602) q[1];
sx q[1];
rz(-0.7936002) q[1];
sx q[1];
rz(2.0097282) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0431932) q[0];
sx q[0];
rz(-0.88623673) q[0];
sx q[0];
rz(1.8266023) q[0];
x q[1];
rz(2.491729) q[2];
sx q[2];
rz(-2.3520326) q[2];
sx q[2];
rz(-2.2672578) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.060185) q[1];
sx q[1];
rz(-0.8601195) q[1];
sx q[1];
rz(-0.036545444) q[1];
rz(-pi) q[2];
rz(-0.15040654) q[3];
sx q[3];
rz(-1.5751221) q[3];
sx q[3];
rz(-1.9459917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.446283) q[2];
sx q[2];
rz(-0.60574564) q[2];
sx q[2];
rz(0.31965762) q[2];
rz(-2.9122635) q[3];
sx q[3];
rz(-1.0234443) q[3];
sx q[3];
rz(2.2314609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0730154) q[0];
sx q[0];
rz(-1.1443161) q[0];
sx q[0];
rz(1.2314433) q[0];
rz(-0.07864174) q[1];
sx q[1];
rz(-1.6784724) q[1];
sx q[1];
rz(2.9873649) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.820136) q[0];
sx q[0];
rz(-0.38561441) q[0];
sx q[0];
rz(0.48954757) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1604105) q[2];
sx q[2];
rz(-1.2912116) q[2];
sx q[2];
rz(-1.2321763) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.63111607) q[1];
sx q[1];
rz(-1.0923903) q[1];
sx q[1];
rz(1.9892742) q[1];
rz(-pi) q[2];
rz(2.6811872) q[3];
sx q[3];
rz(-0.45259991) q[3];
sx q[3];
rz(-0.79595882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.33588931) q[2];
sx q[2];
rz(-1.5259537) q[2];
sx q[2];
rz(1.4914782) q[2];
rz(1.7283745) q[3];
sx q[3];
rz(-1.7468942) q[3];
sx q[3];
rz(1.570805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042628057) q[0];
sx q[0];
rz(-1.093981) q[0];
sx q[0];
rz(1.4549103) q[0];
rz(-2.1658354) q[1];
sx q[1];
rz(-1.059633) q[1];
sx q[1];
rz(-1.3386493) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8737301) q[0];
sx q[0];
rz(-2.3288245) q[0];
sx q[0];
rz(-1.0229179) q[0];
rz(-pi) q[1];
rz(1.067131) q[2];
sx q[2];
rz(-1.2758453) q[2];
sx q[2];
rz(0.3857715) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.83709675) q[1];
sx q[1];
rz(-0.66667914) q[1];
sx q[1];
rz(-1.2231989) q[1];
rz(1.9429132) q[3];
sx q[3];
rz(-0.94148472) q[3];
sx q[3];
rz(0.019364186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6300388) q[2];
sx q[2];
rz(-2.126882) q[2];
sx q[2];
rz(-0.83941984) q[2];
rz(-1.9073585) q[3];
sx q[3];
rz(-1.0890361) q[3];
sx q[3];
rz(0.031410005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79427528) q[0];
sx q[0];
rz(-1.7592156) q[0];
sx q[0];
rz(-0.41912249) q[0];
rz(1.4979111) q[1];
sx q[1];
rz(-1.3587911) q[1];
sx q[1];
rz(-2.7083414) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33746359) q[0];
sx q[0];
rz(-0.47352284) q[0];
sx q[0];
rz(1.8329404) q[0];
rz(1.5574607) q[2];
sx q[2];
rz(-0.91382256) q[2];
sx q[2];
rz(-1.2648392) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.097523263) q[1];
sx q[1];
rz(-1.2882003) q[1];
sx q[1];
rz(1.0364012) q[1];
x q[2];
rz(1.786539) q[3];
sx q[3];
rz(-2.1387055) q[3];
sx q[3];
rz(1.9892094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7822632) q[2];
sx q[2];
rz(-1.1507582) q[2];
sx q[2];
rz(2.9456054) q[2];
rz(1.1228784) q[3];
sx q[3];
rz(-0.71176088) q[3];
sx q[3];
rz(-1.6897197) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6561683) q[0];
sx q[0];
rz(-2.4805785) q[0];
sx q[0];
rz(1.0957023) q[0];
rz(2.6307259) q[1];
sx q[1];
rz(-0.26900649) q[1];
sx q[1];
rz(-0.21024545) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71111403) q[0];
sx q[0];
rz(-2.1415882) q[0];
sx q[0];
rz(-0.23352233) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0265507) q[2];
sx q[2];
rz(-0.91918531) q[2];
sx q[2];
rz(0.63331214) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2481459) q[1];
sx q[1];
rz(-2.4706744) q[1];
sx q[1];
rz(-2.4144961) q[1];
rz(-pi) q[2];
rz(-0.56783592) q[3];
sx q[3];
rz(-2.46008) q[3];
sx q[3];
rz(-2.5208254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.11253396) q[2];
sx q[2];
rz(-1.5052648) q[2];
sx q[2];
rz(1.0851592) q[2];
rz(2.8339913) q[3];
sx q[3];
rz(-0.31809536) q[3];
sx q[3];
rz(-0.65348452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6899684) q[0];
sx q[0];
rz(-1.987048) q[0];
sx q[0];
rz(-1.2183627) q[0];
rz(-2.4901509) q[1];
sx q[1];
rz(-2.1919498) q[1];
sx q[1];
rz(0.776074) q[1];
rz(1.5522926) q[2];
sx q[2];
rz(-2.450569) q[2];
sx q[2];
rz(2.8995502) q[2];
rz(-2.2023946) q[3];
sx q[3];
rz(-1.9892577) q[3];
sx q[3];
rz(1.6792959) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
