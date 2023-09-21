OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3061476) q[0];
sx q[0];
rz(-2.4581576) q[0];
sx q[0];
rz(-0.47877065) q[0];
rz(-3.1105644) q[1];
sx q[1];
rz(-1.9801158) q[1];
sx q[1];
rz(-0.64136139) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35533479) q[0];
sx q[0];
rz(-1.9541652) q[0];
sx q[0];
rz(1.3134365) q[0];
rz(-pi) q[1];
x q[1];
rz(0.85513656) q[2];
sx q[2];
rz(-2.5980066) q[2];
sx q[2];
rz(0.91993514) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7176425) q[1];
sx q[1];
rz(-0.91361928) q[1];
sx q[1];
rz(-0.98492019) q[1];
x q[2];
rz(0.62946837) q[3];
sx q[3];
rz(-1.9068204) q[3];
sx q[3];
rz(2.0365086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.59149867) q[2];
sx q[2];
rz(-1.2167565) q[2];
sx q[2];
rz(0.47810289) q[2];
rz(1.6889307) q[3];
sx q[3];
rz(-2.0958459) q[3];
sx q[3];
rz(-2.5959192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1801382) q[0];
sx q[0];
rz(-0.77471662) q[0];
sx q[0];
rz(1.0189198) q[0];
rz(-1.6628751) q[1];
sx q[1];
rz(-2.5264085) q[1];
sx q[1];
rz(-0.63308024) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79257747) q[0];
sx q[0];
rz(-1.4903869) q[0];
sx q[0];
rz(-2.3809459) q[0];
x q[1];
rz(2.3929246) q[2];
sx q[2];
rz(-2.0249172) q[2];
sx q[2];
rz(-1.548896) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9620348) q[1];
sx q[1];
rz(-1.1711367) q[1];
sx q[1];
rz(0.37079294) q[1];
rz(0.31140621) q[3];
sx q[3];
rz(-0.91324556) q[3];
sx q[3];
rz(-2.1962375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7818266) q[2];
sx q[2];
rz(-0.40010139) q[2];
sx q[2];
rz(-3.0241372) q[2];
rz(-0.30101267) q[3];
sx q[3];
rz(-1.3344701) q[3];
sx q[3];
rz(1.3628179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(2.2837219) q[0];
sx q[0];
rz(-0.71325934) q[0];
sx q[0];
rz(0.088407956) q[0];
rz(-2.03405) q[1];
sx q[1];
rz(-0.91845599) q[1];
sx q[1];
rz(0.69082469) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82949626) q[0];
sx q[0];
rz(-1.4179686) q[0];
sx q[0];
rz(-1.494708) q[0];
rz(-pi) q[1];
rz(-2.8969953) q[2];
sx q[2];
rz(-2.3615712) q[2];
sx q[2];
rz(0.56842677) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9835811) q[1];
sx q[1];
rz(-2.9576655) q[1];
sx q[1];
rz(1.4278825) q[1];
rz(2.7615943) q[3];
sx q[3];
rz(-1.899154) q[3];
sx q[3];
rz(-2.2896555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.92418015) q[2];
sx q[2];
rz(-1.2574544) q[2];
sx q[2];
rz(0.10647354) q[2];
rz(1.4364093) q[3];
sx q[3];
rz(-2.7761716) q[3];
sx q[3];
rz(-1.4829372) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0728264) q[0];
sx q[0];
rz(-2.9077353) q[0];
sx q[0];
rz(-0.52247125) q[0];
rz(-0.32896313) q[1];
sx q[1];
rz(-1.6429699) q[1];
sx q[1];
rz(-0.55363384) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41933665) q[0];
sx q[0];
rz(-2.0291078) q[0];
sx q[0];
rz(2.3769747) q[0];
rz(-pi) q[1];
rz(-2.3217818) q[2];
sx q[2];
rz(-0.95633436) q[2];
sx q[2];
rz(0.43624207) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0370889) q[1];
sx q[1];
rz(-1.7040841) q[1];
sx q[1];
rz(2.3549805) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5997361) q[3];
sx q[3];
rz(-2.3674298) q[3];
sx q[3];
rz(0.72284568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.557495) q[2];
sx q[2];
rz(-0.9157052) q[2];
sx q[2];
rz(-0.49475691) q[2];
rz(-2.2385712) q[3];
sx q[3];
rz(-1.5477864) q[3];
sx q[3];
rz(1.8516699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6506127) q[0];
sx q[0];
rz(-2.3968093) q[0];
sx q[0];
rz(-2.0565128) q[0];
rz(0.96013534) q[1];
sx q[1];
rz(-1.9624058) q[1];
sx q[1];
rz(0.18403149) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2011916) q[0];
sx q[0];
rz(-1.2815676) q[0];
sx q[0];
rz(-1.4616696) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8475624) q[2];
sx q[2];
rz(-1.2185316) q[2];
sx q[2];
rz(-0.89154746) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.23136885) q[1];
sx q[1];
rz(-1.1008881) q[1];
sx q[1];
rz(0.652657) q[1];
rz(-pi) q[2];
rz(2.9156978) q[3];
sx q[3];
rz(-0.73879209) q[3];
sx q[3];
rz(-1.5238638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.90157834) q[2];
sx q[2];
rz(-0.34873909) q[2];
sx q[2];
rz(-2.0007755) q[2];
rz(0.42282894) q[3];
sx q[3];
rz(-1.6641649) q[3];
sx q[3];
rz(1.1269349) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35247701) q[0];
sx q[0];
rz(-1.1081835) q[0];
sx q[0];
rz(-2.254803) q[0];
rz(0.24041644) q[1];
sx q[1];
rz(-2.1227032) q[1];
sx q[1];
rz(-0.14850798) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0200955) q[0];
sx q[0];
rz(-1.7317061) q[0];
sx q[0];
rz(-0.73385977) q[0];
x q[1];
rz(1.4609023) q[2];
sx q[2];
rz(-1.7981148) q[2];
sx q[2];
rz(1.1306888) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2827742) q[1];
sx q[1];
rz(-0.28206477) q[1];
sx q[1];
rz(1.9726522) q[1];
rz(-pi) q[2];
rz(0.46621795) q[3];
sx q[3];
rz(-1.2122452) q[3];
sx q[3];
rz(0.31043226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.85577661) q[2];
sx q[2];
rz(-2.6591876) q[2];
sx q[2];
rz(2.6532069) q[2];
rz(-0.48940247) q[3];
sx q[3];
rz(-0.92178744) q[3];
sx q[3];
rz(1.874812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.69328904) q[0];
sx q[0];
rz(-1.8087837) q[0];
sx q[0];
rz(2.3235902) q[0];
rz(-1.8891634) q[1];
sx q[1];
rz(-0.39223448) q[1];
sx q[1];
rz(1.12524) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55114891) q[0];
sx q[0];
rz(-1.4677912) q[0];
sx q[0];
rz(1.1950905) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0411885) q[2];
sx q[2];
rz(-2.8090968) q[2];
sx q[2];
rz(1.5645129) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.352467) q[1];
sx q[1];
rz(-1.8450292) q[1];
sx q[1];
rz(-1.8791566) q[1];
rz(-pi) q[2];
rz(-0.67894499) q[3];
sx q[3];
rz(-2.3514387) q[3];
sx q[3];
rz(-1.3037579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0685048) q[2];
sx q[2];
rz(-2.1606074) q[2];
sx q[2];
rz(-2.4970064) q[2];
rz(1.6623496) q[3];
sx q[3];
rz(-0.95696604) q[3];
sx q[3];
rz(-1.7361599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1709764) q[0];
sx q[0];
rz(-3.0727486) q[0];
sx q[0];
rz(-1.53565) q[0];
rz(-1.9203141) q[1];
sx q[1];
rz(-1.6284643) q[1];
sx q[1];
rz(-2.1309526) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.772951) q[0];
sx q[0];
rz(-2.3803582) q[0];
sx q[0];
rz(-1.3377454) q[0];
rz(-pi) q[1];
rz(-1.0755195) q[2];
sx q[2];
rz(-1.8349378) q[2];
sx q[2];
rz(-1.2656982) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9606564) q[1];
sx q[1];
rz(-1.3410543) q[1];
sx q[1];
rz(0.22682637) q[1];
rz(-pi) q[2];
rz(2.741022) q[3];
sx q[3];
rz(-1.4860324) q[3];
sx q[3];
rz(-2.8469507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5236139) q[2];
sx q[2];
rz(-2.8221059) q[2];
sx q[2];
rz(-0.69407216) q[2];
rz(-2.5726035) q[3];
sx q[3];
rz(-1.4469955) q[3];
sx q[3];
rz(0.93769658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0555608) q[0];
sx q[0];
rz(-1.9984351) q[0];
sx q[0];
rz(-2.6468497) q[0];
rz(-0.61839473) q[1];
sx q[1];
rz(-1.4952375) q[1];
sx q[1];
rz(-0.075597413) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15007524) q[0];
sx q[0];
rz(-0.65438327) q[0];
sx q[0];
rz(-2.0983216) q[0];
rz(1.3952903) q[2];
sx q[2];
rz(-1.5117466) q[2];
sx q[2];
rz(0.67948558) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1765715) q[1];
sx q[1];
rz(-0.88639835) q[1];
sx q[1];
rz(1.7040764) q[1];
rz(-pi) q[2];
rz(-0.62065403) q[3];
sx q[3];
rz(-1.7194347) q[3];
sx q[3];
rz(1.312048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.33346924) q[2];
sx q[2];
rz(-1.7227017) q[2];
sx q[2];
rz(1.8010275) q[2];
rz(-0.30424413) q[3];
sx q[3];
rz(-0.86383581) q[3];
sx q[3];
rz(-0.58469599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(1.8463523) q[0];
sx q[0];
rz(-1.8363991) q[0];
sx q[0];
rz(0.17679581) q[0];
rz(-1.2416174) q[1];
sx q[1];
rz(-0.63443628) q[1];
sx q[1];
rz(0.44100824) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.037017578) q[0];
sx q[0];
rz(-1.772176) q[0];
sx q[0];
rz(1.3637533) q[0];
rz(-pi) q[1];
rz(-2.6920325) q[2];
sx q[2];
rz(-0.170378) q[2];
sx q[2];
rz(-2.5296488) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0032138) q[1];
sx q[1];
rz(-0.56204501) q[1];
sx q[1];
rz(-2.2467062) q[1];
rz(-pi) q[2];
rz(-2.7367758) q[3];
sx q[3];
rz(-2.4768618) q[3];
sx q[3];
rz(2.3354195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.56197721) q[2];
sx q[2];
rz(-2.5728971) q[2];
sx q[2];
rz(-0.30187541) q[2];
rz(-2.2484696) q[3];
sx q[3];
rz(-1.8538657) q[3];
sx q[3];
rz(-2.9368403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72398913) q[0];
sx q[0];
rz(-1.8483193) q[0];
sx q[0];
rz(1.666477) q[0];
rz(-0.026731116) q[1];
sx q[1];
rz(-1.6865128) q[1];
sx q[1];
rz(-1.7105688) q[1];
rz(-2.6731861) q[2];
sx q[2];
rz(-2.240934) q[2];
sx q[2];
rz(0.59443867) q[2];
rz(2.0704913) q[3];
sx q[3];
rz(-1.0715967) q[3];
sx q[3];
rz(-1.788492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
