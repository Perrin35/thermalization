OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7005641) q[0];
sx q[0];
rz(-1.9987885) q[0];
sx q[0];
rz(-1.9300652) q[0];
rz(2.9149574) q[1];
sx q[1];
rz(-1.5645138) q[1];
sx q[1];
rz(-0.29830631) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5814712) q[0];
sx q[0];
rz(-1.5157962) q[0];
sx q[0];
rz(2.4762857) q[0];
rz(-pi) q[1];
rz(-0.093703336) q[2];
sx q[2];
rz(-1.9022577) q[2];
sx q[2];
rz(1.9625488) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.71516192) q[1];
sx q[1];
rz(-0.75725812) q[1];
sx q[1];
rz(-2.2992579) q[1];
rz(-pi) q[2];
rz(-2.2536623) q[3];
sx q[3];
rz(-3.008932) q[3];
sx q[3];
rz(2.5761029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6478708) q[2];
sx q[2];
rz(-1.2080668) q[2];
sx q[2];
rz(-2.1477264) q[2];
rz(-2.1422051) q[3];
sx q[3];
rz(-1.2402273) q[3];
sx q[3];
rz(0.58888155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64269972) q[0];
sx q[0];
rz(-1.9357341) q[0];
sx q[0];
rz(3.1233741) q[0];
rz(-0.81623626) q[1];
sx q[1];
rz(-2.1111592) q[1];
sx q[1];
rz(-0.47168628) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1175849) q[0];
sx q[0];
rz(-1.4837259) q[0];
sx q[0];
rz(-2.9136806) q[0];
rz(-pi) q[1];
rz(2.6366028) q[2];
sx q[2];
rz(-1.8171176) q[2];
sx q[2];
rz(-1.9976975) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.7447409) q[1];
sx q[1];
rz(-2.5137915) q[1];
sx q[1];
rz(2.6746034) q[1];
rz(-pi) q[2];
rz(0.15700335) q[3];
sx q[3];
rz(-0.96884851) q[3];
sx q[3];
rz(-0.70037819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.54962426) q[2];
sx q[2];
rz(-1.1529808) q[2];
sx q[2];
rz(-0.96898752) q[2];
rz(2.5668868) q[3];
sx q[3];
rz(-0.55137268) q[3];
sx q[3];
rz(2.1000752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7085003) q[0];
sx q[0];
rz(-2.0860465) q[0];
sx q[0];
rz(0.18181268) q[0];
rz(-2.0388942) q[1];
sx q[1];
rz(-1.5010553) q[1];
sx q[1];
rz(-1.4556494) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2244814) q[0];
sx q[0];
rz(-2.3492976) q[0];
sx q[0];
rz(0.86865058) q[0];
rz(-2.4007912) q[2];
sx q[2];
rz(-1.2503137) q[2];
sx q[2];
rz(1.6522811) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0297444) q[1];
sx q[1];
rz(-1.6964456) q[1];
sx q[1];
rz(1.5602342) q[1];
rz(-2.9946795) q[3];
sx q[3];
rz(-2.5111755) q[3];
sx q[3];
rz(-0.86236766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.87749798) q[2];
sx q[2];
rz(-1.654518) q[2];
sx q[2];
rz(-2.1739615) q[2];
rz(2.4140221) q[3];
sx q[3];
rz(-1.2604159) q[3];
sx q[3];
rz(-0.23770604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85686344) q[0];
sx q[0];
rz(-2.6155222) q[0];
sx q[0];
rz(-0.63823429) q[0];
rz(1.1278641) q[1];
sx q[1];
rz(-2.3141839) q[1];
sx q[1];
rz(1.2329873) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5947943) q[0];
sx q[0];
rz(-1.6398755) q[0];
sx q[0];
rz(-2.6462206) q[0];
x q[1];
rz(2.4013176) q[2];
sx q[2];
rz(-2.4057655) q[2];
sx q[2];
rz(-2.1528113) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.14703688) q[1];
sx q[1];
rz(-1.5347267) q[1];
sx q[1];
rz(1.7493164) q[1];
x q[2];
rz(-1.8423555) q[3];
sx q[3];
rz(-2.4304667) q[3];
sx q[3];
rz(-1.7600876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2234852) q[2];
sx q[2];
rz(-1.2780259) q[2];
sx q[2];
rz(2.3275862) q[2];
rz(2.0984086) q[3];
sx q[3];
rz(-0.63101763) q[3];
sx q[3];
rz(-1.957318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84213132) q[0];
sx q[0];
rz(-1.3548387) q[0];
sx q[0];
rz(2.2498851) q[0];
rz(1.2437598) q[1];
sx q[1];
rz(-1.3777106) q[1];
sx q[1];
rz(2.9290501) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5355797) q[0];
sx q[0];
rz(-1.5545462) q[0];
sx q[0];
rz(-1.5505126) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96617713) q[2];
sx q[2];
rz(-1.6024616) q[2];
sx q[2];
rz(2.4633173) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6032431) q[1];
sx q[1];
rz(-1.0298567) q[1];
sx q[1];
rz(-0.69640883) q[1];
rz(2.0810633) q[3];
sx q[3];
rz(-0.79499309) q[3];
sx q[3];
rz(1.0802964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.033096878) q[2];
sx q[2];
rz(-0.97110811) q[2];
sx q[2];
rz(0.5212211) q[2];
rz(1.7565578) q[3];
sx q[3];
rz(-1.9135467) q[3];
sx q[3];
rz(1.2683755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2824771) q[0];
sx q[0];
rz(-0.92882597) q[0];
sx q[0];
rz(-2.916472) q[0];
rz(1.3549995) q[1];
sx q[1];
rz(-1.0083818) q[1];
sx q[1];
rz(-2.7640142) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9281884) q[0];
sx q[0];
rz(-1.8367935) q[0];
sx q[0];
rz(3.1187952) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5224783) q[2];
sx q[2];
rz(-1.2836576) q[2];
sx q[2];
rz(1.5887807) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1212335) q[1];
sx q[1];
rz(-1.6755783) q[1];
sx q[1];
rz(-1.1377513) q[1];
rz(-pi) q[2];
x q[2];
rz(0.79490957) q[3];
sx q[3];
rz(-1.5921633) q[3];
sx q[3];
rz(-1.4710609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.98465115) q[2];
sx q[2];
rz(-2.3539383) q[2];
sx q[2];
rz(2.6605576) q[2];
rz(0.40361079) q[3];
sx q[3];
rz(-1.0389046) q[3];
sx q[3];
rz(0.31479442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0890546) q[0];
sx q[0];
rz(-2.5367694) q[0];
sx q[0];
rz(-0.19454923) q[0];
rz(2.9220707) q[1];
sx q[1];
rz(-1.4621282) q[1];
sx q[1];
rz(2.887168) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1628111) q[0];
sx q[0];
rz(-2.8023976) q[0];
sx q[0];
rz(-1.9041054) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9872586) q[2];
sx q[2];
rz(-2.1774204) q[2];
sx q[2];
rz(-0.91147214) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.19941482) q[1];
sx q[1];
rz(-0.53655469) q[1];
sx q[1];
rz(0.55120991) q[1];
rz(0.80612225) q[3];
sx q[3];
rz(-0.61180173) q[3];
sx q[3];
rz(-0.34477371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1640132) q[2];
sx q[2];
rz(-1.3914725) q[2];
sx q[2];
rz(-1.3158201) q[2];
rz(-2.2655462) q[3];
sx q[3];
rz(-3.0026569) q[3];
sx q[3];
rz(-2.1379437) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9309689) q[0];
sx q[0];
rz(-0.3759149) q[0];
sx q[0];
rz(1.6865431) q[0];
rz(-0.82398206) q[1];
sx q[1];
rz(-2.1897557) q[1];
sx q[1];
rz(1.5664068) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2306838) q[0];
sx q[0];
rz(-2.3678603) q[0];
sx q[0];
rz(0.45549972) q[0];
x q[1];
rz(2.1112061) q[2];
sx q[2];
rz(-1.724913) q[2];
sx q[2];
rz(-0.87481462) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.52289256) q[1];
sx q[1];
rz(-0.52378264) q[1];
sx q[1];
rz(-0.79843847) q[1];
rz(-2.4616562) q[3];
sx q[3];
rz(-2.5222062) q[3];
sx q[3];
rz(-3.1114651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.87551293) q[2];
sx q[2];
rz(-1.9078887) q[2];
sx q[2];
rz(0.27754647) q[2];
rz(-1.6905486) q[3];
sx q[3];
rz(-0.45193672) q[3];
sx q[3];
rz(1.014876) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9797416) q[0];
sx q[0];
rz(-2.6888872) q[0];
sx q[0];
rz(1.45654) q[0];
rz(-0.62943554) q[1];
sx q[1];
rz(-1.1673735) q[1];
sx q[1];
rz(-1.1368407) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0658543) q[0];
sx q[0];
rz(-0.59959164) q[0];
sx q[0];
rz(2.0980741) q[0];
rz(-2.2898229) q[2];
sx q[2];
rz(-0.53817828) q[2];
sx q[2];
rz(0.9966419) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.34502621) q[1];
sx q[1];
rz(-1.476164) q[1];
sx q[1];
rz(1.8887397) q[1];
rz(-pi) q[2];
rz(1.9433446) q[3];
sx q[3];
rz(-2.2714943) q[3];
sx q[3];
rz(-0.31853279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.64951605) q[2];
sx q[2];
rz(-1.4294383) q[2];
sx q[2];
rz(-1.0409522) q[2];
rz(0.0020290931) q[3];
sx q[3];
rz(-0.40015951) q[3];
sx q[3];
rz(-0.31203976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3025538) q[0];
sx q[0];
rz(-2.9170687) q[0];
sx q[0];
rz(0.94605207) q[0];
rz(0.91167766) q[1];
sx q[1];
rz(-1.9263575) q[1];
sx q[1];
rz(-2.5295703) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8004868) q[0];
sx q[0];
rz(-0.90011156) q[0];
sx q[0];
rz(2.1675046) q[0];
rz(2.7962748) q[2];
sx q[2];
rz(-1.7974263) q[2];
sx q[2];
rz(-0.40030865) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.260173) q[1];
sx q[1];
rz(-0.91671413) q[1];
sx q[1];
rz(1.7263078) q[1];
rz(2.6066783) q[3];
sx q[3];
rz(-1.2423008) q[3];
sx q[3];
rz(1.9096149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0570021) q[2];
sx q[2];
rz(-0.64208639) q[2];
sx q[2];
rz(0.65336147) q[2];
rz(-2.7907794) q[3];
sx q[3];
rz(-1.6143129) q[3];
sx q[3];
rz(2.4408834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6012797) q[0];
sx q[0];
rz(-1.606034) q[0];
sx q[0];
rz(0.10869797) q[0];
rz(2.3868949) q[1];
sx q[1];
rz(-1.3194059) q[1];
sx q[1];
rz(-1.5059765) q[1];
rz(2.4621261) q[2];
sx q[2];
rz(-2.0308528) q[2];
sx q[2];
rz(0.46926342) q[2];
rz(-3.0572206) q[3];
sx q[3];
rz(-1.2481239) q[3];
sx q[3];
rz(0.39831755) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];