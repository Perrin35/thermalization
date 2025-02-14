OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71896267) q[0];
sx q[0];
rz(-0.29932061) q[0];
sx q[0];
rz(0.49459767) q[0];
rz(1.142113) q[1];
sx q[1];
rz(-1.0057058) q[1];
sx q[1];
rz(-2.0118654) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62837961) q[0];
sx q[0];
rz(-2.2404379) q[0];
sx q[0];
rz(-1.7731401) q[0];
x q[1];
rz(-2.0689993) q[2];
sx q[2];
rz(-1.1668201) q[2];
sx q[2];
rz(-2.7369268) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3898824) q[1];
sx q[1];
rz(-2.2851508) q[1];
sx q[1];
rz(-2.8044279) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7143728) q[3];
sx q[3];
rz(-2.1764206) q[3];
sx q[3];
rz(2.1744414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.31892458) q[2];
sx q[2];
rz(-1.319004) q[2];
sx q[2];
rz(0.85282105) q[2];
rz(1.8114113) q[3];
sx q[3];
rz(-0.69555247) q[3];
sx q[3];
rz(0.046796355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80832076) q[0];
sx q[0];
rz(-0.051017314) q[0];
sx q[0];
rz(-1.6774696) q[0];
rz(1.6630215) q[1];
sx q[1];
rz(-1.1590978) q[1];
sx q[1];
rz(1.2082072) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20044662) q[0];
sx q[0];
rz(-2.1028554) q[0];
sx q[0];
rz(-1.0131628) q[0];
x q[1];
rz(-2.0457532) q[2];
sx q[2];
rz(-1.7647226) q[2];
sx q[2];
rz(2.3467807) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3390926) q[1];
sx q[1];
rz(-0.29401699) q[1];
sx q[1];
rz(2.0014928) q[1];
x q[2];
rz(-0.34388108) q[3];
sx q[3];
rz(-2.6626427) q[3];
sx q[3];
rz(0.44287455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6055484) q[2];
sx q[2];
rz(-2.7216585) q[2];
sx q[2];
rz(1.4844683) q[2];
rz(0.015497192) q[3];
sx q[3];
rz(-1.213538) q[3];
sx q[3];
rz(-1.1789471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.991796) q[0];
sx q[0];
rz(-1.8114256) q[0];
sx q[0];
rz(-2.8699744) q[0];
rz(-0.89667165) q[1];
sx q[1];
rz(-0.50183693) q[1];
sx q[1];
rz(1.8439878) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5974242) q[0];
sx q[0];
rz(-0.46314683) q[0];
sx q[0];
rz(1.5482272) q[0];
rz(-pi) q[1];
rz(-2.6372725) q[2];
sx q[2];
rz(-1.2337304) q[2];
sx q[2];
rz(0.45318174) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.43186346) q[1];
sx q[1];
rz(-1.6846906) q[1];
sx q[1];
rz(1.3224056) q[1];
rz(-1.5657229) q[3];
sx q[3];
rz(-0.80883615) q[3];
sx q[3];
rz(3.0938597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.74035949) q[2];
sx q[2];
rz(-0.77784246) q[2];
sx q[2];
rz(3.0878301) q[2];
rz(1.2939804) q[3];
sx q[3];
rz(-0.79837489) q[3];
sx q[3];
rz(2.9062041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9409222) q[0];
sx q[0];
rz(-2.5735452) q[0];
sx q[0];
rz(-1.4642375) q[0];
rz(-2.0109406) q[1];
sx q[1];
rz(-0.73431763) q[1];
sx q[1];
rz(3.0272223) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76246156) q[0];
sx q[0];
rz(-1.987559) q[0];
sx q[0];
rz(2.5083191) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0691772) q[2];
sx q[2];
rz(-0.97626462) q[2];
sx q[2];
rz(-1.6785113) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.76941608) q[1];
sx q[1];
rz(-1.5443364) q[1];
sx q[1];
rz(-0.89511223) q[1];
x q[2];
rz(-3.1088282) q[3];
sx q[3];
rz(-0.19812852) q[3];
sx q[3];
rz(1.1968544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6681246) q[2];
sx q[2];
rz(-1.6056085) q[2];
sx q[2];
rz(-0.075695666) q[2];
rz(0.43478742) q[3];
sx q[3];
rz(-1.9861168) q[3];
sx q[3];
rz(-3.0619612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39744034) q[0];
sx q[0];
rz(-1.3120774) q[0];
sx q[0];
rz(-2.721526) q[0];
rz(2.9282667) q[1];
sx q[1];
rz(-1.408564) q[1];
sx q[1];
rz(-1.2423645) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39109502) q[0];
sx q[0];
rz(-0.89078442) q[0];
sx q[0];
rz(2.1260421) q[0];
rz(2.5630997) q[2];
sx q[2];
rz(-0.50225406) q[2];
sx q[2];
rz(-1.945221) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8458357) q[1];
sx q[1];
rz(-1.6713665) q[1];
sx q[1];
rz(-0.015458903) q[1];
x q[2];
rz(-0.22132921) q[3];
sx q[3];
rz(-1.6175272) q[3];
sx q[3];
rz(2.5525301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3199978) q[2];
sx q[2];
rz(-2.2264806) q[2];
sx q[2];
rz(0.37459174) q[2];
rz(-2.1290667) q[3];
sx q[3];
rz(-2.7495224) q[3];
sx q[3];
rz(-0.64468002) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0300765) q[0];
sx q[0];
rz(-2.6426297) q[0];
sx q[0];
rz(-0.43055713) q[0];
rz(2.997609) q[1];
sx q[1];
rz(-1.0824243) q[1];
sx q[1];
rz(-0.63124257) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4743606) q[0];
sx q[0];
rz(-2.1788414) q[0];
sx q[0];
rz(1.6365746) q[0];
x q[1];
rz(2.9955533) q[2];
sx q[2];
rz(-1.1542873) q[2];
sx q[2];
rz(-3.0768968) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5908554) q[1];
sx q[1];
rz(-0.8181347) q[1];
sx q[1];
rz(2.4356151) q[1];
rz(2.1926375) q[3];
sx q[3];
rz(-2.6027711) q[3];
sx q[3];
rz(-0.99297374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1708019) q[2];
sx q[2];
rz(-1.0589212) q[2];
sx q[2];
rz(1.7010752) q[2];
rz(-1.3614281) q[3];
sx q[3];
rz(-1.1037339) q[3];
sx q[3];
rz(2.6543999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.454527) q[0];
sx q[0];
rz(-2.8519958) q[0];
sx q[0];
rz(0.92426306) q[0];
rz(-2.848792) q[1];
sx q[1];
rz(-1.2288789) q[1];
sx q[1];
rz(2.3407095) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070904254) q[0];
sx q[0];
rz(-1.3137806) q[0];
sx q[0];
rz(1.1852253) q[0];
rz(-pi) q[1];
rz(2.7960294) q[2];
sx q[2];
rz(-0.93831944) q[2];
sx q[2];
rz(-1.8008302) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.427721) q[1];
sx q[1];
rz(-1.6000976) q[1];
sx q[1];
rz(-3.0364365) q[1];
x q[2];
rz(0.8893251) q[3];
sx q[3];
rz(-3.1066537) q[3];
sx q[3];
rz(0.068471758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.822829) q[2];
sx q[2];
rz(-1.2206565) q[2];
sx q[2];
rz(1.2804383) q[2];
rz(-0.66796962) q[3];
sx q[3];
rz(-0.48798713) q[3];
sx q[3];
rz(-2.7282696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8845344) q[0];
sx q[0];
rz(-1.9030544) q[0];
sx q[0];
rz(1.9588233) q[0];
rz(-0.23712748) q[1];
sx q[1];
rz(-3.0393937) q[1];
sx q[1];
rz(3.0822486) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8069495) q[0];
sx q[0];
rz(-2.84253) q[0];
sx q[0];
rz(1.1537667) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58299033) q[2];
sx q[2];
rz(-0.81524476) q[2];
sx q[2];
rz(-1.5992523) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4044734) q[1];
sx q[1];
rz(-1.8343907) q[1];
sx q[1];
rz(0.10837676) q[1];
rz(-1.343325) q[3];
sx q[3];
rz(-2.8099647) q[3];
sx q[3];
rz(2.22081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5817029) q[2];
sx q[2];
rz(-2.6010332) q[2];
sx q[2];
rz(-1.7891368) q[2];
rz(1.7743568) q[3];
sx q[3];
rz(-1.7188027) q[3];
sx q[3];
rz(2.2860315) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2818114) q[0];
sx q[0];
rz(-0.78524041) q[0];
sx q[0];
rz(2.6819041) q[0];
rz(0.028060878) q[1];
sx q[1];
rz(-1.1496081) q[1];
sx q[1];
rz(1.185816) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5494649) q[0];
sx q[0];
rz(-0.26927265) q[0];
sx q[0];
rz(-0.16945355) q[0];
x q[1];
rz(-2.2345951) q[2];
sx q[2];
rz(-0.12842783) q[2];
sx q[2];
rz(-0.80824404) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1255155) q[1];
sx q[1];
rz(-1.0927769) q[1];
sx q[1];
rz(2.1090871) q[1];
x q[2];
rz(1.8529296) q[3];
sx q[3];
rz(-2.1698639) q[3];
sx q[3];
rz(-1.4144858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.49843732) q[2];
sx q[2];
rz(-2.2648621) q[2];
sx q[2];
rz(0.15677491) q[2];
rz(-0.59569851) q[3];
sx q[3];
rz(-0.34606338) q[3];
sx q[3];
rz(-0.025207635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6237685) q[0];
sx q[0];
rz(-2.1110004) q[0];
sx q[0];
rz(0.18381707) q[0];
rz(0.078016438) q[1];
sx q[1];
rz(-2.0842431) q[1];
sx q[1];
rz(0.26430166) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3027356) q[0];
sx q[0];
rz(-2.502901) q[0];
sx q[0];
rz(-0.46052082) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61923142) q[2];
sx q[2];
rz(-0.78699099) q[2];
sx q[2];
rz(-0.84112043) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7122927) q[1];
sx q[1];
rz(-1.7279062) q[1];
sx q[1];
rz(0.18204851) q[1];
x q[2];
rz(-1.4749583) q[3];
sx q[3];
rz(-2.0639827) q[3];
sx q[3];
rz(-3.1214023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.8514303) q[2];
sx q[2];
rz(-0.94238472) q[2];
sx q[2];
rz(-2.9128722) q[2];
rz(1.4043572) q[3];
sx q[3];
rz(-0.93969932) q[3];
sx q[3];
rz(0.082988113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9901154) q[0];
sx q[0];
rz(-0.84828068) q[0];
sx q[0];
rz(-1.620851) q[0];
rz(0.71113853) q[1];
sx q[1];
rz(-1.8322721) q[1];
sx q[1];
rz(-2.4087404) q[1];
rz(-0.044471519) q[2];
sx q[2];
rz(-1.0723249) q[2];
sx q[2];
rz(1.9981801) q[2];
rz(0.76821297) q[3];
sx q[3];
rz(-2.429001) q[3];
sx q[3];
rz(2.902239) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
