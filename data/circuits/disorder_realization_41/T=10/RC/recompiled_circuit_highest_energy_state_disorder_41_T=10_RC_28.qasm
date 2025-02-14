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
rz(0.32349411) q[0];
sx q[0];
rz(2.9383724) q[0];
sx q[0];
rz(9.1246224) q[0];
rz(-0.35429859) q[1];
sx q[1];
rz(4.3759182) q[1];
sx q[1];
rz(10.195923) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3187689) q[0];
sx q[0];
rz(-1.7576801) q[0];
sx q[0];
rz(-0.19362886) q[0];
x q[1];
rz(-1.096719) q[2];
sx q[2];
rz(-1.0869622) q[2];
sx q[2];
rz(-0.7686309) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9950152) q[1];
sx q[1];
rz(-1.163401) q[1];
sx q[1];
rz(-2.7316774) q[1];
x q[2];
rz(-0.067581108) q[3];
sx q[3];
rz(-1.1376808) q[3];
sx q[3];
rz(-1.6855296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1511718) q[2];
sx q[2];
rz(-0.46151084) q[2];
sx q[2];
rz(2.0578461) q[2];
rz(-2.566973) q[3];
sx q[3];
rz(-1.5950404) q[3];
sx q[3];
rz(2.3641724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8271178) q[0];
sx q[0];
rz(-2.5125393) q[0];
sx q[0];
rz(-2.7431059) q[0];
rz(1.9972948) q[1];
sx q[1];
rz(-0.60595787) q[1];
sx q[1];
rz(-2.1872637) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1935611) q[0];
sx q[0];
rz(-1.5662417) q[0];
sx q[0];
rz(1.5713552) q[0];
rz(-pi) q[1];
rz(1.4492839) q[2];
sx q[2];
rz(-0.6612311) q[2];
sx q[2];
rz(1.1053156) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1283537) q[1];
sx q[1];
rz(-1.3283821) q[1];
sx q[1];
rz(-0.98636143) q[1];
x q[2];
rz(2.6000836) q[3];
sx q[3];
rz(-1.9810105) q[3];
sx q[3];
rz(2.4334986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2362471) q[2];
sx q[2];
rz(-0.50062847) q[2];
sx q[2];
rz(-0.11032571) q[2];
rz(-2.3162383) q[3];
sx q[3];
rz(-2.3767411) q[3];
sx q[3];
rz(0.74025214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9731307) q[0];
sx q[0];
rz(-0.2178807) q[0];
sx q[0];
rz(2.2595898) q[0];
rz(-1.0285671) q[1];
sx q[1];
rz(-2.5879637) q[1];
sx q[1];
rz(2.2291768) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7275302) q[0];
sx q[0];
rz(-1.4986218) q[0];
sx q[0];
rz(-2.0230002) q[0];
x q[1];
rz(1.4503612) q[2];
sx q[2];
rz(-2.2856973) q[2];
sx q[2];
rz(0.038106136) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0076651) q[1];
sx q[1];
rz(-1.0093736) q[1];
sx q[1];
rz(-2.7473161) q[1];
x q[2];
rz(-1.0718143) q[3];
sx q[3];
rz(-0.49919617) q[3];
sx q[3];
rz(-0.62455356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.87938157) q[2];
sx q[2];
rz(-2.8935581) q[2];
sx q[2];
rz(2.3233419) q[2];
rz(-0.36000559) q[3];
sx q[3];
rz(-1.2986978) q[3];
sx q[3];
rz(-3.103783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70571947) q[0];
sx q[0];
rz(-0.95054764) q[0];
sx q[0];
rz(-1.9807504) q[0];
rz(-2.1707161) q[1];
sx q[1];
rz(-0.91578805) q[1];
sx q[1];
rz(-0.11347778) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1176096) q[0];
sx q[0];
rz(-1.1679839) q[0];
sx q[0];
rz(-0.65448032) q[0];
rz(-1.9275437) q[2];
sx q[2];
rz(-3.103802) q[2];
sx q[2];
rz(-0.99810696) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.66333713) q[1];
sx q[1];
rz(-0.89794822) q[1];
sx q[1];
rz(0.70391432) q[1];
rz(0.63469074) q[3];
sx q[3];
rz(-1.029976) q[3];
sx q[3];
rz(-1.6010923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0299782) q[2];
sx q[2];
rz(-1.8054211) q[2];
sx q[2];
rz(0.90799904) q[2];
rz(2.8152483) q[3];
sx q[3];
rz(-0.55485266) q[3];
sx q[3];
rz(-2.7197796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.570785) q[0];
sx q[0];
rz(-2.3796005) q[0];
sx q[0];
rz(2.1245891) q[0];
rz(0.48655888) q[1];
sx q[1];
rz(-2.1162972) q[1];
sx q[1];
rz(0.09495458) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5240898) q[0];
sx q[0];
rz(-1.4727778) q[0];
sx q[0];
rz(1.6785511) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7899029) q[2];
sx q[2];
rz(-0.75828248) q[2];
sx q[2];
rz(-2.3996283) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.44014318) q[1];
sx q[1];
rz(-0.4491764) q[1];
sx q[1];
rz(2.4692332) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2566936) q[3];
sx q[3];
rz(-1.2746547) q[3];
sx q[3];
rz(1.2552798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.93593705) q[2];
sx q[2];
rz(-0.63465261) q[2];
sx q[2];
rz(-0.46979365) q[2];
rz(-0.69822407) q[3];
sx q[3];
rz(-1.9152812) q[3];
sx q[3];
rz(-2.8214112) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53223377) q[0];
sx q[0];
rz(-2.4025752) q[0];
sx q[0];
rz(-2.9267689) q[0];
rz(2.9048982) q[1];
sx q[1];
rz(-0.57699811) q[1];
sx q[1];
rz(2.7223041) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0146862) q[0];
sx q[0];
rz(-1.6492515) q[0];
sx q[0];
rz(3.0959652) q[0];
rz(2.111666) q[2];
sx q[2];
rz(-1.3020143) q[2];
sx q[2];
rz(0.93712419) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5621369) q[1];
sx q[1];
rz(-0.39135763) q[1];
sx q[1];
rz(2.1298903) q[1];
rz(-pi) q[2];
rz(-2.4287796) q[3];
sx q[3];
rz(-1.4199004) q[3];
sx q[3];
rz(-0.11980443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.99570167) q[2];
sx q[2];
rz(-0.13549165) q[2];
sx q[2];
rz(-3.0502012) q[2];
rz(2.7898096) q[3];
sx q[3];
rz(-2.3942949) q[3];
sx q[3];
rz(0.465213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2202989) q[0];
sx q[0];
rz(-1.1722246) q[0];
sx q[0];
rz(1.9660796) q[0];
rz(2.562404) q[1];
sx q[1];
rz(-2.3349031) q[1];
sx q[1];
rz(-0.18956345) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3407758) q[0];
sx q[0];
rz(-0.22332668) q[0];
sx q[0];
rz(-0.45036611) q[0];
x q[1];
rz(1.4681508) q[2];
sx q[2];
rz(-0.70404875) q[2];
sx q[2];
rz(-1.2873905) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71566641) q[1];
sx q[1];
rz(-2.3235011) q[1];
sx q[1];
rz(-0.010424213) q[1];
rz(-pi) q[2];
rz(2.0912254) q[3];
sx q[3];
rz(-1.5185906) q[3];
sx q[3];
rz(-2.5676651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.2438948) q[2];
sx q[2];
rz(-1.8420668) q[2];
sx q[2];
rz(2.6356836) q[2];
rz(-0.77887744) q[3];
sx q[3];
rz(-2.7454822) q[3];
sx q[3];
rz(1.8469384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8706354) q[0];
sx q[0];
rz(-1.0519692) q[0];
sx q[0];
rz(1.9860995) q[0];
rz(-0.0011860154) q[1];
sx q[1];
rz(-2.3540034) q[1];
sx q[1];
rz(-2.7364065) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74086607) q[0];
sx q[0];
rz(-2.4080546) q[0];
sx q[0];
rz(-2.4057968) q[0];
x q[1];
rz(-2.9878163) q[2];
sx q[2];
rz(-1.9755873) q[2];
sx q[2];
rz(-2.5396233) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2470513) q[1];
sx q[1];
rz(-2.3741268) q[1];
sx q[1];
rz(-1.6805499) q[1];
x q[2];
rz(1.6393597) q[3];
sx q[3];
rz(-1.4430832) q[3];
sx q[3];
rz(0.26940108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.81795168) q[2];
sx q[2];
rz(-1.1672856) q[2];
sx q[2];
rz(2.6200068) q[2];
rz(0.096605435) q[3];
sx q[3];
rz(-1.01869) q[3];
sx q[3];
rz(0.49083138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(2.7535962) q[0];
sx q[0];
rz(-2.2612408) q[0];
sx q[0];
rz(3.0269347) q[0];
rz(1.8278587) q[1];
sx q[1];
rz(-1.4559682) q[1];
sx q[1];
rz(-3.0065261) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11896597) q[0];
sx q[0];
rz(-0.34017902) q[0];
sx q[0];
rz(2.4426798) q[0];
rz(-pi) q[1];
rz(-0.58456771) q[2];
sx q[2];
rz(-2.6815139) q[2];
sx q[2];
rz(2.8321617) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.8986021) q[1];
sx q[1];
rz(-0.95546104) q[1];
sx q[1];
rz(2.9427285) q[1];
x q[2];
rz(-1.4536132) q[3];
sx q[3];
rz(-2.5796267) q[3];
sx q[3];
rz(-2.9148852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8257398) q[2];
sx q[2];
rz(-1.2154546) q[2];
sx q[2];
rz(-2.1960171) q[2];
rz(-0.25997508) q[3];
sx q[3];
rz(-2.1187449) q[3];
sx q[3];
rz(-1.1280577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0876227) q[0];
sx q[0];
rz(-2.051351) q[0];
sx q[0];
rz(-1.4029652) q[0];
rz(2.2258017) q[1];
sx q[1];
rz(-2.5251838) q[1];
sx q[1];
rz(2.2260407) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2820194) q[0];
sx q[0];
rz(-1.440319) q[0];
sx q[0];
rz(-0.097921485) q[0];
rz(-pi) q[1];
rz(-0.19277566) q[2];
sx q[2];
rz(-1.0571684) q[2];
sx q[2];
rz(-2.7889268) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8232223) q[1];
sx q[1];
rz(-2.530066) q[1];
sx q[1];
rz(-2.3182858) q[1];
rz(-pi) q[2];
rz(0.1696945) q[3];
sx q[3];
rz(-0.86830074) q[3];
sx q[3];
rz(2.9164721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0182858) q[2];
sx q[2];
rz(-2.4985963) q[2];
sx q[2];
rz(-2.7389738) q[2];
rz(-1.1766524) q[3];
sx q[3];
rz(-2.8549356) q[3];
sx q[3];
rz(-2.3736242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9993512) q[0];
sx q[0];
rz(-0.9724697) q[0];
sx q[0];
rz(-1.019626) q[0];
rz(-0.6877407) q[1];
sx q[1];
rz(-2.3383457) q[1];
sx q[1];
rz(1.8170423) q[1];
rz(0.68647142) q[2];
sx q[2];
rz(-2.1778637) q[2];
sx q[2];
rz(-1.7354497) q[2];
rz(-0.34360828) q[3];
sx q[3];
rz(-2.4818729) q[3];
sx q[3];
rz(-0.81036405) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
