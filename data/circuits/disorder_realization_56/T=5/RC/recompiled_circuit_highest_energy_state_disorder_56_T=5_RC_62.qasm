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
rz(1.4752969) q[0];
sx q[0];
rz(1.8721606) q[0];
sx q[0];
rz(8.7526487) q[0];
rz(0.44828662) q[1];
sx q[1];
rz(-1.5223794) q[1];
sx q[1];
rz(0.7575922) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31296988) q[0];
sx q[0];
rz(-1.4951473) q[0];
sx q[0];
rz(1.6806169) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.640921) q[2];
sx q[2];
rz(-1.8953875) q[2];
sx q[2];
rz(2.2364834) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.014908286) q[1];
sx q[1];
rz(-1.9796625) q[1];
sx q[1];
rz(1.8168628) q[1];
x q[2];
rz(0.69656482) q[3];
sx q[3];
rz(-1.0961354) q[3];
sx q[3];
rz(-1.153095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0932833) q[2];
sx q[2];
rz(-1.3292686) q[2];
sx q[2];
rz(-1.7993571) q[2];
rz(-1.3679158) q[3];
sx q[3];
rz(-2.048309) q[3];
sx q[3];
rz(-1.5526519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20464483) q[0];
sx q[0];
rz(-1.0138252) q[0];
sx q[0];
rz(-2.192705) q[0];
rz(-0.10143796) q[1];
sx q[1];
rz(-2.0815492) q[1];
sx q[1];
rz(0.92996517) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0596493) q[0];
sx q[0];
rz(-2.905272) q[0];
sx q[0];
rz(0.44775072) q[0];
rz(0.013596046) q[2];
sx q[2];
rz(-2.6363723) q[2];
sx q[2];
rz(-1.1824974) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9252214) q[1];
sx q[1];
rz(-2.9489046) q[1];
sx q[1];
rz(-2.5641901) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92949246) q[3];
sx q[3];
rz(-2.2713823) q[3];
sx q[3];
rz(1.6456347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0019504) q[2];
sx q[2];
rz(-1.4652239) q[2];
sx q[2];
rz(0.742221) q[2];
rz(-2.6774075) q[3];
sx q[3];
rz(-1.7964541) q[3];
sx q[3];
rz(1.5884885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4699698) q[0];
sx q[0];
rz(-2.9556584) q[0];
sx q[0];
rz(-1.4416913) q[0];
rz(0.16547671) q[1];
sx q[1];
rz(-2.4022357) q[1];
sx q[1];
rz(-2.2742719) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4246068) q[0];
sx q[0];
rz(-3.1410257) q[0];
sx q[0];
rz(-1.6419771) q[0];
x q[1];
rz(1.0652415) q[2];
sx q[2];
rz(-1.8095513) q[2];
sx q[2];
rz(2.3468034) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8738385) q[1];
sx q[1];
rz(-1.8513362) q[1];
sx q[1];
rz(-0.33919097) q[1];
x q[2];
rz(-0.8115143) q[3];
sx q[3];
rz(-1.9841188) q[3];
sx q[3];
rz(2.1275525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9637588) q[2];
sx q[2];
rz(-1.9671665) q[2];
sx q[2];
rz(-2.3826694) q[2];
rz(2.3376076) q[3];
sx q[3];
rz(-0.94954973) q[3];
sx q[3];
rz(2.4132531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8300962) q[0];
sx q[0];
rz(-1.281597) q[0];
sx q[0];
rz(-0.48402825) q[0];
rz(-1.6054035) q[1];
sx q[1];
rz(-1.4151298) q[1];
sx q[1];
rz(1.4253634) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7361476) q[0];
sx q[0];
rz(-0.63129866) q[0];
sx q[0];
rz(2.5456356) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31230782) q[2];
sx q[2];
rz(-1.3538401) q[2];
sx q[2];
rz(1.4258476) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0565383) q[1];
sx q[1];
rz(-0.94496545) q[1];
sx q[1];
rz(2.8686348) q[1];
rz(-pi) q[2];
rz(1.1953765) q[3];
sx q[3];
rz(-1.1470801) q[3];
sx q[3];
rz(2.8590607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.385685) q[2];
sx q[2];
rz(-2.0734831) q[2];
sx q[2];
rz(0.29328406) q[2];
rz(3.0934603) q[3];
sx q[3];
rz(-1.8354514) q[3];
sx q[3];
rz(2.8369246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3047979) q[0];
sx q[0];
rz(-1.4076819) q[0];
sx q[0];
rz(-1.1037214) q[0];
rz(-0.91148218) q[1];
sx q[1];
rz(-1.6171004) q[1];
sx q[1];
rz(-3.0326861) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8965204) q[0];
sx q[0];
rz(-2.2811416) q[0];
sx q[0];
rz(0.65585391) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2997024) q[2];
sx q[2];
rz(-0.55520081) q[2];
sx q[2];
rz(-0.19367684) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8629294) q[1];
sx q[1];
rz(-2.4554774) q[1];
sx q[1];
rz(-1.3115694) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2690684) q[3];
sx q[3];
rz(-1.0142039) q[3];
sx q[3];
rz(0.8100141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.64985046) q[2];
sx q[2];
rz(-1.3567341) q[2];
sx q[2];
rz(0.25807992) q[2];
rz(-2.7598925) q[3];
sx q[3];
rz(-0.93770599) q[3];
sx q[3];
rz(-0.55735731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7447164) q[0];
sx q[0];
rz(-2.1602186) q[0];
sx q[0];
rz(-1.1599524) q[0];
rz(-0.57811919) q[1];
sx q[1];
rz(-1.6696397) q[1];
sx q[1];
rz(-0.32346183) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2787298) q[0];
sx q[0];
rz(-0.6375618) q[0];
sx q[0];
rz(-0.80018534) q[0];
rz(1.8238964) q[2];
sx q[2];
rz(-1.0876812) q[2];
sx q[2];
rz(-2.7716669) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66993078) q[1];
sx q[1];
rz(-2.1403229) q[1];
sx q[1];
rz(-1.3013775) q[1];
x q[2];
rz(-1.2488854) q[3];
sx q[3];
rz(-1.1309147) q[3];
sx q[3];
rz(-1.3706576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.20299992) q[2];
sx q[2];
rz(-0.77363571) q[2];
sx q[2];
rz(-0.8482376) q[2];
rz(1.2674468) q[3];
sx q[3];
rz(-1.7402612) q[3];
sx q[3];
rz(0.69376865) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0271725) q[0];
sx q[0];
rz(-0.65374756) q[0];
sx q[0];
rz(2.2681336) q[0];
rz(1.9050441) q[1];
sx q[1];
rz(-0.97573391) q[1];
sx q[1];
rz(0.083560856) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1655052) q[0];
sx q[0];
rz(-0.3437316) q[0];
sx q[0];
rz(-0.36970871) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9026743) q[2];
sx q[2];
rz(-1.0537801) q[2];
sx q[2];
rz(2.1674402) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.04490964) q[1];
sx q[1];
rz(-2.2538141) q[1];
sx q[1];
rz(-2.6234496) q[1];
rz(-pi) q[2];
rz(-0.6912937) q[3];
sx q[3];
rz(-2.6861827) q[3];
sx q[3];
rz(2.5728855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.43625912) q[2];
sx q[2];
rz(-1.6062364) q[2];
sx q[2];
rz(1.3667038) q[2];
rz(-2.8691835) q[3];
sx q[3];
rz(-2.1209769) q[3];
sx q[3];
rz(-2.3830856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90683872) q[0];
sx q[0];
rz(-0.82643569) q[0];
sx q[0];
rz(-2.2032264) q[0];
rz(0.80288184) q[1];
sx q[1];
rz(-0.96790853) q[1];
sx q[1];
rz(-2.3209007) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.11242) q[0];
sx q[0];
rz(-0.39390644) q[0];
sx q[0];
rz(1.2159416) q[0];
rz(1.9782045) q[2];
sx q[2];
rz(-0.40336868) q[2];
sx q[2];
rz(-2.5591116) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8093687) q[1];
sx q[1];
rz(-1.4648155) q[1];
sx q[1];
rz(2.1952704) q[1];
x q[2];
rz(1.9441767) q[3];
sx q[3];
rz(-2.7146517) q[3];
sx q[3];
rz(2.3641242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.20251033) q[2];
sx q[2];
rz(-2.9896917) q[2];
sx q[2];
rz(2.0738156) q[2];
rz(-2.0104525) q[3];
sx q[3];
rz(-0.83429566) q[3];
sx q[3];
rz(0.080032674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15602569) q[0];
sx q[0];
rz(-1.1556867) q[0];
sx q[0];
rz(0.40618968) q[0];
rz(1.73229) q[1];
sx q[1];
rz(-1.0180165) q[1];
sx q[1];
rz(-1.0850151) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70124088) q[0];
sx q[0];
rz(-0.95537649) q[0];
sx q[0];
rz(1.3568272) q[0];
x q[1];
rz(0.0064956587) q[2];
sx q[2];
rz(-2.5312739) q[2];
sx q[2];
rz(-1.6808866) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.93714) q[1];
sx q[1];
rz(-1.2431743) q[1];
sx q[1];
rz(0.83468584) q[1];
rz(-2.1187339) q[3];
sx q[3];
rz(-2.0798488) q[3];
sx q[3];
rz(1.7804543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5558418) q[2];
sx q[2];
rz(-2.5297574) q[2];
sx q[2];
rz(2.2108868) q[2];
rz(1.7454923) q[3];
sx q[3];
rz(-1.3627108) q[3];
sx q[3];
rz(3.0220368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68468204) q[0];
sx q[0];
rz(-0.97701183) q[0];
sx q[0];
rz(-1.3071625) q[0];
rz(0.3859418) q[1];
sx q[1];
rz(-1.394505) q[1];
sx q[1];
rz(-2.1070259) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2292896) q[0];
sx q[0];
rz(-1.467839) q[0];
sx q[0];
rz(-1.3961755) q[0];
rz(0.39888403) q[2];
sx q[2];
rz(-2.6475057) q[2];
sx q[2];
rz(-1.2436155) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.99920995) q[1];
sx q[1];
rz(-2.2015328) q[1];
sx q[1];
rz(-0.56908619) q[1];
x q[2];
rz(0.12328445) q[3];
sx q[3];
rz(-0.81466952) q[3];
sx q[3];
rz(-0.88501634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4474386) q[2];
sx q[2];
rz(-0.73221451) q[2];
sx q[2];
rz(-2.7745957) q[2];
rz(-1.1191818) q[3];
sx q[3];
rz(-2.7435591) q[3];
sx q[3];
rz(1.5252349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3601396) q[0];
sx q[0];
rz(-1.1893138) q[0];
sx q[0];
rz(-1.0839533) q[0];
rz(-3.0837334) q[1];
sx q[1];
rz(-1.2670988) q[1];
sx q[1];
rz(1.7115464) q[1];
rz(-0.51390263) q[2];
sx q[2];
rz(-0.19795098) q[2];
sx q[2];
rz(1.2951938) q[2];
rz(-0.69969768) q[3];
sx q[3];
rz(-1.2917923) q[3];
sx q[3];
rz(2.7421799) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
