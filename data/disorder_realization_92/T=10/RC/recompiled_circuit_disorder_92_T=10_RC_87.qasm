OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6264412) q[0];
sx q[0];
rz(-3.1080973) q[0];
sx q[0];
rz(-1.7749696) q[0];
rz(-1.051149) q[1];
sx q[1];
rz(-1.4895952) q[1];
sx q[1];
rz(1.1319914) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7819408) q[0];
sx q[0];
rz(-1.9160761) q[0];
sx q[0];
rz(-2.3495673) q[0];
x q[1];
rz(-2.3425383) q[2];
sx q[2];
rz(-0.22944268) q[2];
sx q[2];
rz(0.45861751) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.72343091) q[1];
sx q[1];
rz(-1.9730113) q[1];
sx q[1];
rz(0.9546141) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2855929) q[3];
sx q[3];
rz(-2.1072901) q[3];
sx q[3];
rz(-2.5453018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.91036096) q[2];
sx q[2];
rz(-1.3277162) q[2];
sx q[2];
rz(1.988391) q[2];
rz(-2.6575346) q[3];
sx q[3];
rz(-0.81367937) q[3];
sx q[3];
rz(2.6822065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83950481) q[0];
sx q[0];
rz(-2.8110101) q[0];
sx q[0];
rz(-0.50305811) q[0];
rz(-1.5548276) q[1];
sx q[1];
rz(-2.4350872) q[1];
sx q[1];
rz(2.9876626) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.670823) q[0];
sx q[0];
rz(-1.5564859) q[0];
sx q[0];
rz(1.0679246) q[0];
rz(-pi) q[1];
rz(2.2203127) q[2];
sx q[2];
rz(-1.3037455) q[2];
sx q[2];
rz(2.9608179) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0347552) q[1];
sx q[1];
rz(-1.991786) q[1];
sx q[1];
rz(1.0616598) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1400938) q[3];
sx q[3];
rz(-2.631819) q[3];
sx q[3];
rz(-1.203376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8872035) q[2];
sx q[2];
rz(-0.35787359) q[2];
sx q[2];
rz(1.8187693) q[2];
rz(-1.4860738) q[3];
sx q[3];
rz(-1.5406939) q[3];
sx q[3];
rz(0.71050182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94674295) q[0];
sx q[0];
rz(-1.1179593) q[0];
sx q[0];
rz(-0.90993607) q[0];
rz(2.3643156) q[1];
sx q[1];
rz(-0.83559075) q[1];
sx q[1];
rz(0.98532239) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1811718) q[0];
sx q[0];
rz(-2.4665678) q[0];
sx q[0];
rz(1.8280562) q[0];
rz(-2.8172242) q[2];
sx q[2];
rz(-0.93967059) q[2];
sx q[2];
rz(2.1567675) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9039771) q[1];
sx q[1];
rz(-1.8198038) q[1];
sx q[1];
rz(-1.6275089) q[1];
rz(0.91089532) q[3];
sx q[3];
rz(-1.0064555) q[3];
sx q[3];
rz(0.49163715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2788006) q[2];
sx q[2];
rz(-1.9868439) q[2];
sx q[2];
rz(-1.8939691) q[2];
rz(-2.6990081) q[3];
sx q[3];
rz(-1.3675888) q[3];
sx q[3];
rz(2.8201593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9001532) q[0];
sx q[0];
rz(-2.1223919) q[0];
sx q[0];
rz(-2.0181657) q[0];
rz(-0.57888794) q[1];
sx q[1];
rz(-1.6826948) q[1];
sx q[1];
rz(1.7480063) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9959895) q[0];
sx q[0];
rz(-1.9847426) q[0];
sx q[0];
rz(-2.4664509) q[0];
rz(-pi) q[1];
x q[1];
rz(0.43487866) q[2];
sx q[2];
rz(-1.7231427) q[2];
sx q[2];
rz(-0.16532126) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3634998) q[1];
sx q[1];
rz(-2.2908195) q[1];
sx q[1];
rz(-2.5219445) q[1];
x q[2];
rz(1.8279151) q[3];
sx q[3];
rz(-2.3895279) q[3];
sx q[3];
rz(2.5079692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1018155) q[2];
sx q[2];
rz(-1.5856051) q[2];
sx q[2];
rz(-0.0021136443) q[2];
rz(-0.56143108) q[3];
sx q[3];
rz(-2.1241472) q[3];
sx q[3];
rz(0.58825618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.026022) q[0];
sx q[0];
rz(-2.0176812) q[0];
sx q[0];
rz(-1.2258688) q[0];
rz(-1.6745802) q[1];
sx q[1];
rz(-1.2743228) q[1];
sx q[1];
rz(1.7747169) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0910049) q[0];
sx q[0];
rz(-2.707621) q[0];
sx q[0];
rz(2.6758053) q[0];
rz(-pi) q[1];
rz(-0.032614313) q[2];
sx q[2];
rz(-0.74748749) q[2];
sx q[2];
rz(-0.05687296) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.67147672) q[1];
sx q[1];
rz(-1.4911545) q[1];
sx q[1];
rz(-2.2786042) q[1];
rz(-pi) q[2];
rz(0.50232356) q[3];
sx q[3];
rz(-1.2548496) q[3];
sx q[3];
rz(-1.5918921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.86429578) q[2];
sx q[2];
rz(-1.0884476) q[2];
sx q[2];
rz(2.7362291) q[2];
rz(2.6799485) q[3];
sx q[3];
rz(-0.82563892) q[3];
sx q[3];
rz(1.595114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49155238) q[0];
sx q[0];
rz(-1.7164282) q[0];
sx q[0];
rz(-2.6382085) q[0];
rz(2.9227496) q[1];
sx q[1];
rz(-1.8914521) q[1];
sx q[1];
rz(0.65178451) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1056021) q[0];
sx q[0];
rz(-2.833617) q[0];
sx q[0];
rz(0.67291798) q[0];
x q[1];
rz(-1.4621554) q[2];
sx q[2];
rz(-2.3896653) q[2];
sx q[2];
rz(1.7679364) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8004605) q[1];
sx q[1];
rz(-0.144185) q[1];
sx q[1];
rz(-1.9393117) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7995293) q[3];
sx q[3];
rz(-0.29933375) q[3];
sx q[3];
rz(0.39810668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.51263222) q[2];
sx q[2];
rz(-1.7444538) q[2];
sx q[2];
rz(0.39166489) q[2];
rz(0.20646778) q[3];
sx q[3];
rz(-0.73409096) q[3];
sx q[3];
rz(-2.4519043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3702635) q[0];
sx q[0];
rz(-1.4498793) q[0];
sx q[0];
rz(-2.1719334) q[0];
rz(2.5580653) q[1];
sx q[1];
rz(-2.0136166) q[1];
sx q[1];
rz(-0.13024174) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1286344) q[0];
sx q[0];
rz(-2.7854714) q[0];
sx q[0];
rz(0.14333368) q[0];
x q[1];
rz(2.4728751) q[2];
sx q[2];
rz(-1.0895184) q[2];
sx q[2];
rz(0.95578335) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3245557) q[1];
sx q[1];
rz(-0.55700028) q[1];
sx q[1];
rz(2.7213098) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98826615) q[3];
sx q[3];
rz(-1.8997972) q[3];
sx q[3];
rz(2.695801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.51817259) q[2];
sx q[2];
rz(-1.5815846) q[2];
sx q[2];
rz(-2.0557892) q[2];
rz(-0.052224934) q[3];
sx q[3];
rz(-1.4445392) q[3];
sx q[3];
rz(2.0558555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4645585) q[0];
sx q[0];
rz(-0.98419404) q[0];
sx q[0];
rz(-1.1664671) q[0];
rz(0.72558609) q[1];
sx q[1];
rz(-1.2936932) q[1];
sx q[1];
rz(0.74277791) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60177875) q[0];
sx q[0];
rz(-0.83866461) q[0];
sx q[0];
rz(-1.8307122) q[0];
rz(2.4648415) q[2];
sx q[2];
rz(-1.1427715) q[2];
sx q[2];
rz(1.8523491) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.038866) q[1];
sx q[1];
rz(-2.3042149) q[1];
sx q[1];
rz(2.1329761) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5524213) q[3];
sx q[3];
rz(-0.99654752) q[3];
sx q[3];
rz(1.1816927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.62721884) q[2];
sx q[2];
rz(-1.5292995) q[2];
sx q[2];
rz(0.26088866) q[2];
rz(-2.0339113) q[3];
sx q[3];
rz(-1.8543782) q[3];
sx q[3];
rz(-2.0598944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7850007) q[0];
sx q[0];
rz(-2.3605425) q[0];
sx q[0];
rz(2.2055431) q[0];
rz(3.127457) q[1];
sx q[1];
rz(-1.3052669) q[1];
sx q[1];
rz(0.65151185) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0885568) q[0];
sx q[0];
rz(-3.1175107) q[0];
sx q[0];
rz(1.3911029) q[0];
x q[1];
rz(0.65666171) q[2];
sx q[2];
rz(-1.0495532) q[2];
sx q[2];
rz(1.9929287) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7266255) q[1];
sx q[1];
rz(-1.8888998) q[1];
sx q[1];
rz(-0.12673641) q[1];
rz(-pi) q[2];
rz(-0.43990073) q[3];
sx q[3];
rz(-1.2370046) q[3];
sx q[3];
rz(-1.4378289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2953879) q[2];
sx q[2];
rz(-2.442895) q[2];
sx q[2];
rz(-2.6182168) q[2];
rz(-2.8335617) q[3];
sx q[3];
rz(-2.568646) q[3];
sx q[3];
rz(-1.8035005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(1.734252) q[0];
sx q[0];
rz(-2.9946406) q[0];
sx q[0];
rz(1.7379606) q[0];
rz(-2.667528) q[1];
sx q[1];
rz(-1.6876551) q[1];
sx q[1];
rz(1.7636991) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4353367) q[0];
sx q[0];
rz(-2.113111) q[0];
sx q[0];
rz(-0.52973912) q[0];
rz(-2.2950298) q[2];
sx q[2];
rz(-0.91869527) q[2];
sx q[2];
rz(0.50154274) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.024152012) q[1];
sx q[1];
rz(-1.1785058) q[1];
sx q[1];
rz(-1.3082318) q[1];
rz(-pi) q[2];
x q[2];
rz(0.5573427) q[3];
sx q[3];
rz(-1.7094269) q[3];
sx q[3];
rz(-1.0728474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3906117) q[2];
sx q[2];
rz(-1.5079974) q[2];
sx q[2];
rz(-1.5314468) q[2];
rz(-1.4577929) q[3];
sx q[3];
rz(-1.0554353) q[3];
sx q[3];
rz(-2.2916601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4703341) q[0];
sx q[0];
rz(-0.27161921) q[0];
sx q[0];
rz(-2.0441396) q[0];
rz(-2.509027) q[1];
sx q[1];
rz(-1.0610776) q[1];
sx q[1];
rz(-2.9656596) q[1];
rz(-2.73545) q[2];
sx q[2];
rz(-0.71229013) q[2];
sx q[2];
rz(-2.9980414) q[2];
rz(-0.73344161) q[3];
sx q[3];
rz(-1.0167132) q[3];
sx q[3];
rz(0.50501962) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
