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
rz(0.23213586) q[0];
sx q[0];
rz(-0.27096662) q[0];
sx q[0];
rz(-1.1722857) q[0];
rz(1.5341893) q[1];
sx q[1];
rz(-1.6521896) q[1];
sx q[1];
rz(0.54816562) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0294103) q[0];
sx q[0];
rz(-1.2971086) q[0];
sx q[0];
rz(3.0269483) q[0];
x q[1];
rz(-2.1287969) q[2];
sx q[2];
rz(-1.1269953) q[2];
sx q[2];
rz(-1.9091801) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1480268) q[1];
sx q[1];
rz(-1.2669592) q[1];
sx q[1];
rz(0.5263473) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0788623) q[3];
sx q[3];
rz(-2.0196805) q[3];
sx q[3];
rz(-2.1470438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.24070172) q[2];
sx q[2];
rz(-1.9343932) q[2];
sx q[2];
rz(-1.5473676) q[2];
rz(-2.2089925) q[3];
sx q[3];
rz(-2.2280732) q[3];
sx q[3];
rz(0.80476052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43657434) q[0];
sx q[0];
rz(-0.60282928) q[0];
sx q[0];
rz(-0.21827179) q[0];
rz(1.5087347) q[1];
sx q[1];
rz(-0.873133) q[1];
sx q[1];
rz(1.1223209) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2788804) q[0];
sx q[0];
rz(-0.76970184) q[0];
sx q[0];
rz(2.2055726) q[0];
rz(-pi) q[1];
rz(-1.7709184) q[2];
sx q[2];
rz(-0.84017838) q[2];
sx q[2];
rz(-2.6482128) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7699204) q[1];
sx q[1];
rz(-2.3974916) q[1];
sx q[1];
rz(2.4456853) q[1];
rz(-2.2436879) q[3];
sx q[3];
rz(-0.39552415) q[3];
sx q[3];
rz(0.69098587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0108769) q[2];
sx q[2];
rz(-1.9515832) q[2];
sx q[2];
rz(2.7575764) q[2];
rz(2.5859517) q[3];
sx q[3];
rz(-2.5540387) q[3];
sx q[3];
rz(-2.2491992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71690774) q[0];
sx q[0];
rz(-0.45775828) q[0];
sx q[0];
rz(-2.6440788) q[0];
rz(2.6566907) q[1];
sx q[1];
rz(-1.5496016) q[1];
sx q[1];
rz(0.44152322) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6832774) q[0];
sx q[0];
rz(-2.2234869) q[0];
sx q[0];
rz(0.21245942) q[0];
x q[1];
rz(-1.0777801) q[2];
sx q[2];
rz(-1.7642972) q[2];
sx q[2];
rz(-2.5317321) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.71601) q[1];
sx q[1];
rz(-2.0507318) q[1];
sx q[1];
rz(0.77253567) q[1];
rz(1.2526337) q[3];
sx q[3];
rz(-1.3103879) q[3];
sx q[3];
rz(-1.3459149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.32450822) q[2];
sx q[2];
rz(-1.0434979) q[2];
sx q[2];
rz(-0.95477742) q[2];
rz(-0.020141007) q[3];
sx q[3];
rz(-0.79831278) q[3];
sx q[3];
rz(-0.32625833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7666053) q[0];
sx q[0];
rz(-2.6191481) q[0];
sx q[0];
rz(0.0047542714) q[0];
rz(0.44113723) q[1];
sx q[1];
rz(-0.10103592) q[1];
sx q[1];
rz(-1.4099247) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49214464) q[0];
sx q[0];
rz(-0.36055598) q[0];
sx q[0];
rz(1.1900405) q[0];
x q[1];
rz(2.8660315) q[2];
sx q[2];
rz(-2.2803734) q[2];
sx q[2];
rz(-2.1103275) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0995304) q[1];
sx q[1];
rz(-2.5914609) q[1];
sx q[1];
rz(-1.128837) q[1];
x q[2];
rz(3.1143777) q[3];
sx q[3];
rz(-1.7593972) q[3];
sx q[3];
rz(-1.0116497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.29177523) q[2];
sx q[2];
rz(-1.3128277) q[2];
sx q[2];
rz(-2.4449091) q[2];
rz(-1.5289395) q[3];
sx q[3];
rz(-0.63642514) q[3];
sx q[3];
rz(0.67970413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-1.9686389) q[0];
sx q[0];
rz(-0.30690673) q[0];
sx q[0];
rz(2.2343743) q[0];
rz(-0.13547678) q[1];
sx q[1];
rz(-1.2716581) q[1];
sx q[1];
rz(2.4438593) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0402148) q[0];
sx q[0];
rz(-1.7003978) q[0];
sx q[0];
rz(-1.965353) q[0];
x q[1];
rz(-1.582007) q[2];
sx q[2];
rz(-0.19489842) q[2];
sx q[2];
rz(-0.29981183) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5144044) q[1];
sx q[1];
rz(-1.3980734) q[1];
sx q[1];
rz(-0.33337969) q[1];
rz(-pi) q[2];
x q[2];
rz(2.145153) q[3];
sx q[3];
rz(-3.0868885) q[3];
sx q[3];
rz(2.8131054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.1873143) q[2];
sx q[2];
rz(-2.4772187) q[2];
sx q[2];
rz(0.3453671) q[2];
rz(2.6464388) q[3];
sx q[3];
rz(-0.59585714) q[3];
sx q[3];
rz(3.0143747) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0536163) q[0];
sx q[0];
rz(-1.477366) q[0];
sx q[0];
rz(-1.9867058) q[0];
rz(2.3467973) q[1];
sx q[1];
rz(-1.844901) q[1];
sx q[1];
rz(2.6577139) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4876955) q[0];
sx q[0];
rz(-1.9761818) q[0];
sx q[0];
rz(-2.2873108) q[0];
rz(-pi) q[1];
rz(-2.0675428) q[2];
sx q[2];
rz(-1.7109461) q[2];
sx q[2];
rz(-1.4432009) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4166548) q[1];
sx q[1];
rz(-2.0339801) q[1];
sx q[1];
rz(-0.2128667) q[1];
rz(-pi) q[2];
rz(0.5838861) q[3];
sx q[3];
rz(-2.4777081) q[3];
sx q[3];
rz(-1.8107469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.38857073) q[2];
sx q[2];
rz(-1.2349671) q[2];
sx q[2];
rz(2.0335601) q[2];
rz(1.5293416) q[3];
sx q[3];
rz(-3.022091) q[3];
sx q[3];
rz(0.39335462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18688467) q[0];
sx q[0];
rz(-2.055838) q[0];
sx q[0];
rz(2.4672274) q[0];
rz(-2.2652594) q[1];
sx q[1];
rz(-0.74462157) q[1];
sx q[1];
rz(0.032329917) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6276777) q[0];
sx q[0];
rz(-1.6925294) q[0];
sx q[0];
rz(1.7072149) q[0];
x q[1];
rz(-2.9170119) q[2];
sx q[2];
rz(-1.4201284) q[2];
sx q[2];
rz(-1.9407139) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8994042) q[1];
sx q[1];
rz(-0.63540484) q[1];
sx q[1];
rz(-0.59235488) q[1];
rz(-pi) q[2];
rz(-2.1635734) q[3];
sx q[3];
rz(-1.7306574) q[3];
sx q[3];
rz(-2.7474952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.18675599) q[2];
sx q[2];
rz(-2.7080471) q[2];
sx q[2];
rz(-1.0369066) q[2];
rz(-1.2284651) q[3];
sx q[3];
rz(-0.90932536) q[3];
sx q[3];
rz(2.9508446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4036338) q[0];
sx q[0];
rz(-3.1138595) q[0];
sx q[0];
rz(-2.8926358) q[0];
rz(-0.20949334) q[1];
sx q[1];
rz(-1.8003576) q[1];
sx q[1];
rz(0.6627717) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6903721) q[0];
sx q[0];
rz(-1.5826384) q[0];
sx q[0];
rz(-1.5607587) q[0];
x q[1];
rz(-2.3784654) q[2];
sx q[2];
rz(-1.210312) q[2];
sx q[2];
rz(-0.57661001) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1472305) q[1];
sx q[1];
rz(-2.3636345) q[1];
sx q[1];
rz(-1.2176355) q[1];
rz(-pi) q[2];
rz(-2.7484077) q[3];
sx q[3];
rz(-2.4454456) q[3];
sx q[3];
rz(-0.37567155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2304307) q[2];
sx q[2];
rz(-1.5713567) q[2];
sx q[2];
rz(-2.7479808) q[2];
rz(-0.39129928) q[3];
sx q[3];
rz(-0.53520447) q[3];
sx q[3];
rz(0.75214255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1289718) q[0];
sx q[0];
rz(-2.9550645) q[0];
sx q[0];
rz(0.57998002) q[0];
rz(-2.773556) q[1];
sx q[1];
rz(-0.8152222) q[1];
sx q[1];
rz(-3.0267402) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5631487) q[0];
sx q[0];
rz(-0.73860335) q[0];
sx q[0];
rz(-1.369801) q[0];
rz(2.9909953) q[2];
sx q[2];
rz(-2.0080749) q[2];
sx q[2];
rz(2.8018746) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.29909322) q[1];
sx q[1];
rz(-2.7242492) q[1];
sx q[1];
rz(0.6562161) q[1];
rz(-pi) q[2];
rz(2.7389218) q[3];
sx q[3];
rz(-0.91710233) q[3];
sx q[3];
rz(1.5066866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4174691) q[2];
sx q[2];
rz(-0.33895156) q[2];
sx q[2];
rz(-2.4453956) q[2];
rz(2.0255069) q[3];
sx q[3];
rz(-2.3510272) q[3];
sx q[3];
rz(-0.74151403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8490132) q[0];
sx q[0];
rz(-2.1722023) q[0];
sx q[0];
rz(-0.27266362) q[0];
rz(2.4478681) q[1];
sx q[1];
rz(-2.0246181) q[1];
sx q[1];
rz(2.5781217) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7359223) q[0];
sx q[0];
rz(-1.5338681) q[0];
sx q[0];
rz(-1.6640628) q[0];
x q[1];
rz(2.5921037) q[2];
sx q[2];
rz(-0.27590431) q[2];
sx q[2];
rz(0.23978182) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8565791) q[1];
sx q[1];
rz(-1.4895298) q[1];
sx q[1];
rz(1.8829499) q[1];
rz(0.71406096) q[3];
sx q[3];
rz(-2.0787368) q[3];
sx q[3];
rz(1.270961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9844661) q[2];
sx q[2];
rz(-2.0968585) q[2];
sx q[2];
rz(-0.54894471) q[2];
rz(-0.99556154) q[3];
sx q[3];
rz(-0.82058161) q[3];
sx q[3];
rz(-2.5115749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7231049) q[0];
sx q[0];
rz(-2.1514819) q[0];
sx q[0];
rz(2.6958382) q[0];
rz(2.2743629) q[1];
sx q[1];
rz(-1.1031311) q[1];
sx q[1];
rz(1.1581609) q[1];
rz(1.340048) q[2];
sx q[2];
rz(-2.0371155) q[2];
sx q[2];
rz(-2.2178963) q[2];
rz(2.5936962) q[3];
sx q[3];
rz(-1.4461645) q[3];
sx q[3];
rz(-0.89492284) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
