OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3747342) q[0];
sx q[0];
rz(-2.2280966) q[0];
sx q[0];
rz(-1.7537533) q[0];
rz(0.18985441) q[1];
sx q[1];
rz(-1.8228276) q[1];
sx q[1];
rz(0.28611046) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6494516) q[0];
sx q[0];
rz(-1.1104662) q[0];
sx q[0];
rz(1.5440617) q[0];
rz(0.95979664) q[2];
sx q[2];
rz(-1.6382484) q[2];
sx q[2];
rz(-2.5416201) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0148791) q[1];
sx q[1];
rz(-0.61903799) q[1];
sx q[1];
rz(-2.6837885) q[1];
rz(-pi) q[2];
x q[2];
rz(0.055130868) q[3];
sx q[3];
rz(-1.1305439) q[3];
sx q[3];
rz(2.3326479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.046595786) q[2];
sx q[2];
rz(-1.6511714) q[2];
sx q[2];
rz(-0.68520927) q[2];
rz(1.4677216) q[3];
sx q[3];
rz(-1.0823715) q[3];
sx q[3];
rz(2.8732324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.14107038) q[0];
sx q[0];
rz(-1.5474316) q[0];
sx q[0];
rz(2.1970314) q[0];
rz(-3.0682849) q[1];
sx q[1];
rz(-2.3088375) q[1];
sx q[1];
rz(-1.8836969) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8507754) q[0];
sx q[0];
rz(-0.75277872) q[0];
sx q[0];
rz(-0.87429177) q[0];
x q[1];
rz(2.6075105) q[2];
sx q[2];
rz(-1.3943496) q[2];
sx q[2];
rz(-0.29613972) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.91852616) q[1];
sx q[1];
rz(-0.41689532) q[1];
sx q[1];
rz(-2.7781163) q[1];
rz(2.2444751) q[3];
sx q[3];
rz(-0.059487933) q[3];
sx q[3];
rz(0.18732303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.32626095) q[2];
sx q[2];
rz(-1.7504642) q[2];
sx q[2];
rz(-2.1168671) q[2];
rz(-2.4552086) q[3];
sx q[3];
rz(-2.4095583) q[3];
sx q[3];
rz(2.2288403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82655418) q[0];
sx q[0];
rz(-1.8678366) q[0];
sx q[0];
rz(-2.3231373) q[0];
rz(1.2089027) q[1];
sx q[1];
rz(-1.6345638) q[1];
sx q[1];
rz(-2.2817629) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6085109) q[0];
sx q[0];
rz(-1.2405335) q[0];
sx q[0];
rz(0.94655605) q[0];
x q[1];
rz(2.0468065) q[2];
sx q[2];
rz(-0.90887672) q[2];
sx q[2];
rz(0.87023416) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2261994) q[1];
sx q[1];
rz(-0.92297727) q[1];
sx q[1];
rz(-2.1658705) q[1];
x q[2];
rz(1.1016261) q[3];
sx q[3];
rz(-2.8896324) q[3];
sx q[3];
rz(-0.89706883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7439421) q[2];
sx q[2];
rz(-0.33471477) q[2];
sx q[2];
rz(-0.115455) q[2];
rz(-0.3913106) q[3];
sx q[3];
rz(-2.1152928) q[3];
sx q[3];
rz(1.1439884) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7855969) q[0];
sx q[0];
rz(-1.1859256) q[0];
sx q[0];
rz(-0.24056973) q[0];
rz(-1.7754414) q[1];
sx q[1];
rz(-1.6484478) q[1];
sx q[1];
rz(-1.8919401) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6281462) q[0];
sx q[0];
rz(-1.8246987) q[0];
sx q[0];
rz(-2.2280919) q[0];
x q[1];
rz(2.0555105) q[2];
sx q[2];
rz(-2.4668985) q[2];
sx q[2];
rz(2.7759068) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.74339691) q[1];
sx q[1];
rz(-2.2155688) q[1];
sx q[1];
rz(1.866524) q[1];
rz(-1.0650738) q[3];
sx q[3];
rz(-0.95129644) q[3];
sx q[3];
rz(2.2175585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.25632855) q[2];
sx q[2];
rz(-2.5958131) q[2];
sx q[2];
rz(-1.3679999) q[2];
rz(-2.8374425) q[3];
sx q[3];
rz(-1.5757898) q[3];
sx q[3];
rz(0.69653851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2279385) q[0];
sx q[0];
rz(-2.9253687) q[0];
sx q[0];
rz(0.53713334) q[0];
rz(-0.74983239) q[1];
sx q[1];
rz(-1.0822108) q[1];
sx q[1];
rz(-2.4044663) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1340807) q[0];
sx q[0];
rz(-1.5069557) q[0];
sx q[0];
rz(-3.0456187) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3959048) q[2];
sx q[2];
rz(-1.6564257) q[2];
sx q[2];
rz(1.081664) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9530216) q[1];
sx q[1];
rz(-0.28975315) q[1];
sx q[1];
rz(2.6798331) q[1];
rz(0.26925663) q[3];
sx q[3];
rz(-2.807121) q[3];
sx q[3];
rz(2.0967029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4445112) q[2];
sx q[2];
rz(-2.7883174) q[2];
sx q[2];
rz(-2.9949761) q[2];
rz(-1.4491436) q[3];
sx q[3];
rz(-1.861898) q[3];
sx q[3];
rz(0.52391887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.523282) q[0];
sx q[0];
rz(-2.1250516) q[0];
sx q[0];
rz(1.7200394) q[0];
rz(1.8824185) q[1];
sx q[1];
rz(-0.6764532) q[1];
sx q[1];
rz(-2.7377985) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97966563) q[0];
sx q[0];
rz(-2.296519) q[0];
sx q[0];
rz(2.7620074) q[0];
rz(-pi) q[1];
rz(-2.5228259) q[2];
sx q[2];
rz(-0.31224373) q[2];
sx q[2];
rz(-1.1340326) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.67298573) q[1];
sx q[1];
rz(-0.98277175) q[1];
sx q[1];
rz(2.468416) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3293224) q[3];
sx q[3];
rz(-0.89697402) q[3];
sx q[3];
rz(2.633147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8518565) q[2];
sx q[2];
rz(-1.7278262) q[2];
sx q[2];
rz(0.16656052) q[2];
rz(1.5038495) q[3];
sx q[3];
rz(-2.6681191) q[3];
sx q[3];
rz(3.04305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2378167) q[0];
sx q[0];
rz(-0.38919583) q[0];
sx q[0];
rz(2.26407) q[0];
rz(-2.1794043) q[1];
sx q[1];
rz(-2.018237) q[1];
sx q[1];
rz(0.35433623) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49119774) q[0];
sx q[0];
rz(-1.808486) q[0];
sx q[0];
rz(1.5525805) q[0];
rz(-pi) q[1];
rz(-1.7298752) q[2];
sx q[2];
rz(-2.9949778) q[2];
sx q[2];
rz(2.7377812) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-8*pi/15) q[1];
sx q[1];
rz(-2.1051268) q[1];
sx q[1];
rz(-2.2855845) q[1];
x q[2];
rz(-1.6463424) q[3];
sx q[3];
rz(-1.0620688) q[3];
sx q[3];
rz(1.4010324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.80498901) q[2];
sx q[2];
rz(-0.73770928) q[2];
sx q[2];
rz(-1.0450411) q[2];
rz(-2.7392144) q[3];
sx q[3];
rz(-1.1674403) q[3];
sx q[3];
rz(1.6654525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7104915) q[0];
sx q[0];
rz(-3.0248108) q[0];
sx q[0];
rz(-2.2905599) q[0];
rz(2.7878413) q[1];
sx q[1];
rz(-1.4101135) q[1];
sx q[1];
rz(-0.85465777) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8806649) q[0];
sx q[0];
rz(-1.3890084) q[0];
sx q[0];
rz(-2.8559235) q[0];
rz(-2.607367) q[2];
sx q[2];
rz(-2.8534128) q[2];
sx q[2];
rz(2.2792918) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.51616299) q[1];
sx q[1];
rz(-0.99402797) q[1];
sx q[1];
rz(-1.7712084) q[1];
rz(1.6519671) q[3];
sx q[3];
rz(-1.9725091) q[3];
sx q[3];
rz(0.69854004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.19055584) q[2];
sx q[2];
rz(-0.64728105) q[2];
sx q[2];
rz(-0.62038842) q[2];
rz(-2.9151211) q[3];
sx q[3];
rz(-1.7445931) q[3];
sx q[3];
rz(-0.58102077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1802375) q[0];
sx q[0];
rz(-1.2018452) q[0];
sx q[0];
rz(-0.015137976) q[0];
rz(2.119428) q[1];
sx q[1];
rz(-2.731555) q[1];
sx q[1];
rz(-1.2000363) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2101636) q[0];
sx q[0];
rz(-2.4151122) q[0];
sx q[0];
rz(-1.8379712) q[0];
x q[1];
rz(-2.5483918) q[2];
sx q[2];
rz(-1.6199379) q[2];
sx q[2];
rz(-1.0047439) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2851499) q[1];
sx q[1];
rz(-1.0593482) q[1];
sx q[1];
rz(-3.1072381) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0807812) q[3];
sx q[3];
rz(-1.1378764) q[3];
sx q[3];
rz(-2.1149389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6156561) q[2];
sx q[2];
rz(-0.91292149) q[2];
sx q[2];
rz(-2.9021662) q[2];
rz(-3.0827403) q[3];
sx q[3];
rz(-0.81164304) q[3];
sx q[3];
rz(-2.8656901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9489768) q[0];
sx q[0];
rz(-1.2139576) q[0];
sx q[0];
rz(-2.1666727) q[0];
rz(0.67543593) q[1];
sx q[1];
rz(-0.91921872) q[1];
sx q[1];
rz(-0.98178896) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4083207) q[0];
sx q[0];
rz(-1.8309203) q[0];
sx q[0];
rz(-2.1704546) q[0];
rz(-2.9862829) q[2];
sx q[2];
rz(-0.41602817) q[2];
sx q[2];
rz(1.1890026) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2416545) q[1];
sx q[1];
rz(-1.3119427) q[1];
sx q[1];
rz(2.3090825) q[1];
x q[2];
rz(-2.0174867) q[3];
sx q[3];
rz(-1.5367931) q[3];
sx q[3];
rz(-1.407935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5741817) q[2];
sx q[2];
rz(-0.97102204) q[2];
sx q[2];
rz(1.9166454) q[2];
rz(-2.775906) q[3];
sx q[3];
rz(-2.1920125) q[3];
sx q[3];
rz(2.001781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0443307) q[0];
sx q[0];
rz(-1.8987569) q[0];
sx q[0];
rz(1.4415997) q[0];
rz(3.0615831) q[1];
sx q[1];
rz(-1.7623822) q[1];
sx q[1];
rz(-0.0026127876) q[1];
rz(0.2856725) q[2];
sx q[2];
rz(-2.6798669) q[2];
sx q[2];
rz(0.2017022) q[2];
rz(-0.046924165) q[3];
sx q[3];
rz(-2.6886446) q[3];
sx q[3];
rz(-0.038084134) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
