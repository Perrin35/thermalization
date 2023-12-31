OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.055846) q[0];
sx q[0];
rz(-3.0598109) q[0];
sx q[0];
rz(-0.50146377) q[0];
rz(-1.6429098) q[1];
sx q[1];
rz(-0.39615762) q[1];
sx q[1];
rz(0.3224386) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0816406) q[0];
sx q[0];
rz(-1.2498706) q[0];
sx q[0];
rz(0.089846213) q[0];
rz(-pi) q[1];
rz(2.2519977) q[2];
sx q[2];
rz(-2.0686364) q[2];
sx q[2];
rz(1.7878469) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1134125) q[1];
sx q[1];
rz(-2.1537188) q[1];
sx q[1];
rz(0.41863538) q[1];
x q[2];
rz(2.1665855) q[3];
sx q[3];
rz(-2.4132204) q[3];
sx q[3];
rz(-1.7155852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.50513187) q[2];
sx q[2];
rz(-0.59288609) q[2];
sx q[2];
rz(0.55603975) q[2];
rz(-2.3089144) q[3];
sx q[3];
rz(-1.4913538) q[3];
sx q[3];
rz(-0.94579831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6933724) q[0];
sx q[0];
rz(-1.4602666) q[0];
sx q[0];
rz(-0.15727501) q[0];
rz(2.8804624) q[1];
sx q[1];
rz(-1.3477247) q[1];
sx q[1];
rz(-0.10903407) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029023829) q[0];
sx q[0];
rz(-1.7984263) q[0];
sx q[0];
rz(0.17507041) q[0];
rz(-0.26308665) q[2];
sx q[2];
rz(-1.7619942) q[2];
sx q[2];
rz(-0.69743246) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1034531) q[1];
sx q[1];
rz(-0.98587576) q[1];
sx q[1];
rz(-1.000688) q[1];
x q[2];
rz(1.7440967) q[3];
sx q[3];
rz(-2.0850075) q[3];
sx q[3];
rz(-1.7931995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9033501) q[2];
sx q[2];
rz(-1.976333) q[2];
sx q[2];
rz(1.2634574) q[2];
rz(0.3271099) q[3];
sx q[3];
rz(-1.5771022) q[3];
sx q[3];
rz(1.2143149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0771714) q[0];
sx q[0];
rz(-0.049296878) q[0];
sx q[0];
rz(1.3431312) q[0];
rz(2.893977) q[1];
sx q[1];
rz(-0.74664465) q[1];
sx q[1];
rz(0.48167357) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7307229) q[0];
sx q[0];
rz(-2.0559089) q[0];
sx q[0];
rz(-2.7184125) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1842314) q[2];
sx q[2];
rz(-1.295225) q[2];
sx q[2];
rz(0.36188175) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2783918) q[1];
sx q[1];
rz(-0.87676261) q[1];
sx q[1];
rz(-2.4064526) q[1];
rz(-2.2504911) q[3];
sx q[3];
rz(-0.84935969) q[3];
sx q[3];
rz(-0.47618714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8032916) q[2];
sx q[2];
rz(-0.81739601) q[2];
sx q[2];
rz(2.6417007) q[2];
rz(2.5806184) q[3];
sx q[3];
rz(-1.8818972) q[3];
sx q[3];
rz(-1.6803754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19701476) q[0];
sx q[0];
rz(-0.16600969) q[0];
sx q[0];
rz(-2.5894077) q[0];
rz(-1.588297) q[1];
sx q[1];
rz(-2.242656) q[1];
sx q[1];
rz(1.8968556) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23918505) q[0];
sx q[0];
rz(-2.8021078) q[0];
sx q[0];
rz(0.0048588077) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7719901) q[2];
sx q[2];
rz(-1.6987213) q[2];
sx q[2];
rz(0.97845562) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2326395) q[1];
sx q[1];
rz(-2.5563572) q[1];
sx q[1];
rz(-1.9613683) q[1];
x q[2];
rz(-0.57085412) q[3];
sx q[3];
rz(-1.7613162) q[3];
sx q[3];
rz(-1.966147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2923979) q[2];
sx q[2];
rz(-1.8820102) q[2];
sx q[2];
rz(-1.1506895) q[2];
rz(-1.4771279) q[3];
sx q[3];
rz(-1.5090347) q[3];
sx q[3];
rz(-0.48294827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0634336) q[0];
sx q[0];
rz(-0.76197356) q[0];
sx q[0];
rz(-3.0601236) q[0];
rz(-0.062462656) q[1];
sx q[1];
rz(-2.000258) q[1];
sx q[1];
rz(1.6385471) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3969288) q[0];
sx q[0];
rz(-1.443112) q[0];
sx q[0];
rz(-0.98608195) q[0];
rz(-1.4625174) q[2];
sx q[2];
rz(-1.4881954) q[2];
sx q[2];
rz(-0.29512197) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4245783) q[1];
sx q[1];
rz(-1.6728405) q[1];
sx q[1];
rz(-1.8841519) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53346975) q[3];
sx q[3];
rz(-2.5662078) q[3];
sx q[3];
rz(1.433978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.5082671) q[2];
sx q[2];
rz(-2.1990364) q[2];
sx q[2];
rz(1.2379237) q[2];
rz(2.0189019) q[3];
sx q[3];
rz(-2.4653547) q[3];
sx q[3];
rz(0.52156633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6102819) q[0];
sx q[0];
rz(-2.1370482) q[0];
sx q[0];
rz(0.26671985) q[0];
rz(-0.56089127) q[1];
sx q[1];
rz(-1.8436878) q[1];
sx q[1];
rz(-0.7985324) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079433867) q[0];
sx q[0];
rz(-1.2957797) q[0];
sx q[0];
rz(0.14607231) q[0];
x q[1];
rz(-1.2217667) q[2];
sx q[2];
rz(-1.0997699) q[2];
sx q[2];
rz(0.40086056) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.861607) q[1];
sx q[1];
rz(-1.0000739) q[1];
sx q[1];
rz(-0.062203783) q[1];
rz(0.800662) q[3];
sx q[3];
rz(-2.0293529) q[3];
sx q[3];
rz(-0.82160219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9810527) q[2];
sx q[2];
rz(-1.8926228) q[2];
sx q[2];
rz(0.36671656) q[2];
rz(1.2612873) q[3];
sx q[3];
rz(-1.6882608) q[3];
sx q[3];
rz(-3.0453851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40903184) q[0];
sx q[0];
rz(-0.92027396) q[0];
sx q[0];
rz(2.5352056) q[0];
rz(2.9442893) q[1];
sx q[1];
rz(-1.1261255) q[1];
sx q[1];
rz(-2.6775449) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0021792) q[0];
sx q[0];
rz(-2.3521949) q[0];
sx q[0];
rz(-0.95285691) q[0];
x q[1];
rz(3.1214141) q[2];
sx q[2];
rz(-1.8115461) q[2];
sx q[2];
rz(1.9764331) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.11394994) q[1];
sx q[1];
rz(-2.0538035) q[1];
sx q[1];
rz(-0.34160683) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23630996) q[3];
sx q[3];
rz(-1.650562) q[3];
sx q[3];
rz(-2.7304756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2074034) q[2];
sx q[2];
rz(-2.1384017) q[2];
sx q[2];
rz(0.25804538) q[2];
rz(-1.9559654) q[3];
sx q[3];
rz(-1.5153171) q[3];
sx q[3];
rz(-0.0035088249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.946452) q[0];
sx q[0];
rz(-1.2807245) q[0];
sx q[0];
rz(2.7602957) q[0];
rz(0.095245846) q[1];
sx q[1];
rz(-2.1691599) q[1];
sx q[1];
rz(1.7000748) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38948108) q[0];
sx q[0];
rz(-1.507326) q[0];
sx q[0];
rz(1.5556637) q[0];
rz(-pi) q[1];
rz(-2.5289815) q[2];
sx q[2];
rz(-1.6022041) q[2];
sx q[2];
rz(-0.14234662) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.86130202) q[1];
sx q[1];
rz(-0.18096563) q[1];
sx q[1];
rz(-0.7678395) q[1];
rz(-pi) q[2];
rz(-0.12318792) q[3];
sx q[3];
rz(-2.8052969) q[3];
sx q[3];
rz(-1.5817225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9986481) q[2];
sx q[2];
rz(-2.7286077) q[2];
sx q[2];
rz(2.9150035) q[2];
rz(0.44858027) q[3];
sx q[3];
rz(-1.6058763) q[3];
sx q[3];
rz(-2.3118238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7609693) q[0];
sx q[0];
rz(-2.3801104) q[0];
sx q[0];
rz(-1.7425849) q[0];
rz(-0.31708583) q[1];
sx q[1];
rz(-1.6665019) q[1];
sx q[1];
rz(2.1549966) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2664468) q[0];
sx q[0];
rz(-0.52030021) q[0];
sx q[0];
rz(0.32169028) q[0];
rz(-pi) q[1];
rz(-2.1836906) q[2];
sx q[2];
rz(-2.3447678) q[2];
sx q[2];
rz(-2.5755142) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4793195) q[1];
sx q[1];
rz(-2.5179177) q[1];
sx q[1];
rz(0.24346607) q[1];
rz(-pi) q[2];
rz(-1.9916233) q[3];
sx q[3];
rz(-1.4126561) q[3];
sx q[3];
rz(-1.6588253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2150779) q[2];
sx q[2];
rz(-0.72700095) q[2];
sx q[2];
rz(-0.40965664) q[2];
rz(0.26327291) q[3];
sx q[3];
rz(-1.3116838) q[3];
sx q[3];
rz(-0.51945654) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39919329) q[0];
sx q[0];
rz(-3.0629459) q[0];
sx q[0];
rz(1.4051399) q[0];
rz(2.3204904) q[1];
sx q[1];
rz(-0.91870538) q[1];
sx q[1];
rz(1.4155037) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9916519) q[0];
sx q[0];
rz(-1.9016001) q[0];
sx q[0];
rz(-1.3209692) q[0];
rz(-pi) q[1];
rz(1.9344994) q[2];
sx q[2];
rz(-1.3377681) q[2];
sx q[2];
rz(1.9050913) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.532383) q[1];
sx q[1];
rz(-1.4192974) q[1];
sx q[1];
rz(-1.7864368) q[1];
rz(-pi) q[2];
rz(2.9726082) q[3];
sx q[3];
rz(-1.0645234) q[3];
sx q[3];
rz(-0.81912012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4225509) q[2];
sx q[2];
rz(-0.30095235) q[2];
sx q[2];
rz(0.12410513) q[2];
rz(0.96578807) q[3];
sx q[3];
rz(-1.6671168) q[3];
sx q[3];
rz(0.66108274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.538095) q[0];
sx q[0];
rz(-2.8932543) q[0];
sx q[0];
rz(2.2809991) q[0];
rz(-2.8339236) q[1];
sx q[1];
rz(-1.888702) q[1];
sx q[1];
rz(-1.9370334) q[1];
rz(-0.260367) q[2];
sx q[2];
rz(-1.7454864) q[2];
sx q[2];
rz(-2.7228552) q[2];
rz(-1.7210759) q[3];
sx q[3];
rz(-2.1223304) q[3];
sx q[3];
rz(0.35029678) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
