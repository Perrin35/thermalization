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
rz(2.870626) q[0];
sx q[0];
rz(10.597064) q[0];
rz(1.5341893) q[1];
sx q[1];
rz(-1.6521896) q[1];
sx q[1];
rz(0.54816562) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6518657) q[0];
sx q[0];
rz(-1.681156) q[0];
sx q[0];
rz(1.8462028) q[0];
x q[1];
rz(0.83914031) q[2];
sx q[2];
rz(-0.6979894) q[2];
sx q[2];
rz(2.2006387) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.89828789) q[1];
sx q[1];
rz(-0.60051781) q[1];
sx q[1];
rz(-2.5836248) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4413766) q[3];
sx q[3];
rz(-2.6886422) q[3];
sx q[3];
rz(2.2907886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9008909) q[2];
sx q[2];
rz(-1.9343932) q[2];
sx q[2];
rz(-1.594225) q[2];
rz(-0.93260014) q[3];
sx q[3];
rz(-2.2280732) q[3];
sx q[3];
rz(-0.80476052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7050183) q[0];
sx q[0];
rz(-0.60282928) q[0];
sx q[0];
rz(0.21827179) q[0];
rz(1.5087347) q[1];
sx q[1];
rz(-0.873133) q[1];
sx q[1];
rz(-2.0192718) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9199785) q[0];
sx q[0];
rz(-1.99619) q[0];
sx q[0];
rz(2.2334188) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9232736) q[2];
sx q[2];
rz(-0.75262296) q[2];
sx q[2];
rz(2.3531331) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7699204) q[1];
sx q[1];
rz(-2.3974916) q[1];
sx q[1];
rz(2.4456853) q[1];
rz(-pi) q[2];
rz(-1.8863986) q[3];
sx q[3];
rz(-1.3282933) q[3];
sx q[3];
rz(-0.24569233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1307158) q[2];
sx q[2];
rz(-1.1900095) q[2];
sx q[2];
rz(2.7575764) q[2];
rz(2.5859517) q[3];
sx q[3];
rz(-2.5540387) q[3];
sx q[3];
rz(0.89239341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4246849) q[0];
sx q[0];
rz(-2.6838344) q[0];
sx q[0];
rz(0.49751383) q[0];
rz(-0.48490191) q[1];
sx q[1];
rz(-1.5919911) q[1];
sx q[1];
rz(2.7000694) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8988441) q[0];
sx q[0];
rz(-1.4024807) q[0];
sx q[0];
rz(0.90710137) q[0];
rz(-2.0638126) q[2];
sx q[2];
rz(-1.7642972) q[2];
sx q[2];
rz(-0.60986054) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4388112) q[1];
sx q[1];
rz(-0.88248135) q[1];
sx q[1];
rz(-2.5007894) q[1];
rz(0.86526342) q[3];
sx q[3];
rz(-2.7332716) q[3];
sx q[3];
rz(0.43864076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8170844) q[2];
sx q[2];
rz(-1.0434979) q[2];
sx q[2];
rz(-0.95477742) q[2];
rz(-3.1214516) q[3];
sx q[3];
rz(-0.79831278) q[3];
sx q[3];
rz(-2.8153343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37498736) q[0];
sx q[0];
rz(-0.52244455) q[0];
sx q[0];
rz(3.1368384) q[0];
rz(-2.7004554) q[1];
sx q[1];
rz(-0.10103592) q[1];
sx q[1];
rz(-1.4099247) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.649448) q[0];
sx q[0];
rz(-2.7810367) q[0];
sx q[0];
rz(-1.9515522) q[0];
rz(-pi) q[1];
rz(-1.2639764) q[2];
sx q[2];
rz(-0.75245082) q[2];
sx q[2];
rz(-0.6217989) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9121352) q[1];
sx q[1];
rz(-1.3452824) q[1];
sx q[1];
rz(1.0646125) q[1];
x q[2];
rz(-3.1143777) q[3];
sx q[3];
rz(-1.3821954) q[3];
sx q[3];
rz(-1.0116497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8498174) q[2];
sx q[2];
rz(-1.3128277) q[2];
sx q[2];
rz(2.4449091) q[2];
rz(1.6126532) q[3];
sx q[3];
rz(-0.63642514) q[3];
sx q[3];
rz(0.67970413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9686389) q[0];
sx q[0];
rz(-0.30690673) q[0];
sx q[0];
rz(-0.90721834) q[0];
rz(0.13547678) q[1];
sx q[1];
rz(-1.8699346) q[1];
sx q[1];
rz(-0.69773331) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22963069) q[0];
sx q[0];
rz(-2.7273587) q[0];
sx q[0];
rz(-1.8976865) q[0];
x q[1];
rz(-1.582007) q[2];
sx q[2];
rz(-2.9466942) q[2];
sx q[2];
rz(0.29981183) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4042774) q[1];
sx q[1];
rz(-2.7676146) q[1];
sx q[1];
rz(-0.48980063) q[1];
x q[2];
rz(0.99643965) q[3];
sx q[3];
rz(-3.0868885) q[3];
sx q[3];
rz(0.32848725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.1873143) q[2];
sx q[2];
rz(-2.4772187) q[2];
sx q[2];
rz(0.3453671) q[2];
rz(0.49515381) q[3];
sx q[3];
rz(-0.59585714) q[3];
sx q[3];
rz(-3.0143747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0879764) q[0];
sx q[0];
rz(-1.6642267) q[0];
sx q[0];
rz(1.1548868) q[0];
rz(2.3467973) q[1];
sx q[1];
rz(-1.2966917) q[1];
sx q[1];
rz(-2.6577139) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4876955) q[0];
sx q[0];
rz(-1.1654108) q[0];
sx q[0];
rz(0.8542819) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15911289) q[2];
sx q[2];
rz(-1.0793574) q[2];
sx q[2];
rz(0.05201498) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.7495854) q[1];
sx q[1];
rz(-1.760943) q[1];
sx q[1];
rz(2.0431678) q[1];
rz(2.5577066) q[3];
sx q[3];
rz(-0.66388452) q[3];
sx q[3];
rz(-1.8107469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.38857073) q[2];
sx q[2];
rz(-1.9066255) q[2];
sx q[2];
rz(-2.0335601) q[2];
rz(1.5293416) q[3];
sx q[3];
rz(-3.022091) q[3];
sx q[3];
rz(0.39335462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.954708) q[0];
sx q[0];
rz(-1.0857546) q[0];
sx q[0];
rz(2.4672274) q[0];
rz(-0.87633324) q[1];
sx q[1];
rz(-0.74462157) q[1];
sx q[1];
rz(-0.032329917) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5139149) q[0];
sx q[0];
rz(-1.4490632) q[0];
sx q[0];
rz(1.7072149) q[0];
rz(-pi) q[1];
rz(2.9170119) q[2];
sx q[2];
rz(-1.4201284) q[2];
sx q[2];
rz(1.9407139) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2421884) q[1];
sx q[1];
rz(-2.5061878) q[1];
sx q[1];
rz(-0.59235488) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2898212) q[3];
sx q[3];
rz(-2.5301438) q[3];
sx q[3];
rz(0.94463724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.18675599) q[2];
sx q[2];
rz(-0.43354553) q[2];
sx q[2];
rz(2.1046861) q[2];
rz(-1.9131276) q[3];
sx q[3];
rz(-2.2322673) q[3];
sx q[3];
rz(2.9508446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4036338) q[0];
sx q[0];
rz(-3.1138595) q[0];
sx q[0];
rz(2.8926358) q[0];
rz(0.20949334) q[1];
sx q[1];
rz(-1.341235) q[1];
sx q[1];
rz(0.6627717) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4512206) q[0];
sx q[0];
rz(-1.5589542) q[0];
sx q[0];
rz(1.580834) q[0];
x q[1];
rz(2.6423062) q[2];
sx q[2];
rz(-0.82816507) q[2];
sx q[2];
rz(-1.3475095) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.31964916) q[1];
sx q[1];
rz(-1.3256097) q[1];
sx q[1];
rz(2.3169434) q[1];
rz(-pi) q[2];
rz(1.2609188) q[3];
sx q[3];
rz(-0.93684465) q[3];
sx q[3];
rz(3.0216963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2304307) q[2];
sx q[2];
rz(-1.570236) q[2];
sx q[2];
rz(-0.39361185) q[2];
rz(-2.7502934) q[3];
sx q[3];
rz(-2.6063882) q[3];
sx q[3];
rz(0.75214255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0126208) q[0];
sx q[0];
rz(-0.18652815) q[0];
sx q[0];
rz(2.5616126) q[0];
rz(2.773556) q[1];
sx q[1];
rz(-0.8152222) q[1];
sx q[1];
rz(-0.11485242) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84732417) q[0];
sx q[0];
rz(-2.2912187) q[0];
sx q[0];
rz(-2.9617733) q[0];
rz(1.129135) q[2];
sx q[2];
rz(-1.7071305) q[2];
sx q[2];
rz(1.974687) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.29909322) q[1];
sx q[1];
rz(-0.4173435) q[1];
sx q[1];
rz(-0.6562161) q[1];
rz(-pi) q[2];
x q[2];
rz(0.40267085) q[3];
sx q[3];
rz(-0.91710233) q[3];
sx q[3];
rz(-1.5066866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4174691) q[2];
sx q[2];
rz(-0.33895156) q[2];
sx q[2];
rz(-0.69619703) q[2];
rz(-2.0255069) q[3];
sx q[3];
rz(-0.79056549) q[3];
sx q[3];
rz(-0.74151403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.29257947) q[0];
sx q[0];
rz(-2.1722023) q[0];
sx q[0];
rz(0.27266362) q[0];
rz(-2.4478681) q[1];
sx q[1];
rz(-1.1169746) q[1];
sx q[1];
rz(-0.56347096) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7359223) q[0];
sx q[0];
rz(-1.5338681) q[0];
sx q[0];
rz(-1.4775299) q[0];
x q[1];
rz(-2.5921037) q[2];
sx q[2];
rz(-2.8656883) q[2];
sx q[2];
rz(0.23978182) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2595926) q[1];
sx q[1];
rz(-1.881885) q[1];
sx q[1];
rz(-3.056219) q[1];
rz(0.70448168) q[3];
sx q[3];
rz(-0.84957963) q[3];
sx q[3];
rz(0.81126708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9844661) q[2];
sx q[2];
rz(-1.0447341) q[2];
sx q[2];
rz(-0.54894471) q[2];
rz(0.99556154) q[3];
sx q[3];
rz(-2.321011) q[3];
sx q[3];
rz(0.63001776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.7231049) q[0];
sx q[0];
rz(-0.99011078) q[0];
sx q[0];
rz(-0.44575442) q[0];
rz(2.2743629) q[1];
sx q[1];
rz(-1.1031311) q[1];
sx q[1];
rz(1.1581609) q[1];
rz(-0.47719279) q[2];
sx q[2];
rz(-1.365061) q[2];
sx q[2];
rz(-0.75233598) q[2];
rz(0.23602939) q[3];
sx q[3];
rz(-2.5811145) q[3];
sx q[3];
rz(0.47490912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
