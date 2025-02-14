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
rz(1.3951294) q[0];
sx q[0];
rz(-2.3869618) q[0];
sx q[0];
rz(-1.0838497) q[0];
rz(-0.87767449) q[1];
sx q[1];
rz(-1.2134774) q[1];
sx q[1];
rz(-0.47668639) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81945588) q[0];
sx q[0];
rz(-1.1753723) q[0];
sx q[0];
rz(2.8614392) q[0];
rz(-2.0332912) q[2];
sx q[2];
rz(-1.9841737) q[2];
sx q[2];
rz(-2.5819636) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.649316) q[1];
sx q[1];
rz(-1.5536123) q[1];
sx q[1];
rz(1.8284998) q[1];
rz(-pi) q[2];
rz(-2.7101948) q[3];
sx q[3];
rz(-1.3702624) q[3];
sx q[3];
rz(-1.3667184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1276663) q[2];
sx q[2];
rz(-1.4907336) q[2];
sx q[2];
rz(-0.15065436) q[2];
rz(1.6023747) q[3];
sx q[3];
rz(-0.37709245) q[3];
sx q[3];
rz(-0.25130513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056331228) q[0];
sx q[0];
rz(-1.769861) q[0];
sx q[0];
rz(2.7714609) q[0];
rz(0.47259694) q[1];
sx q[1];
rz(-2.8234146) q[1];
sx q[1];
rz(0.95742375) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8214398) q[0];
sx q[0];
rz(-1.6509045) q[0];
sx q[0];
rz(-1.0715241) q[0];
rz(1.9198869) q[2];
sx q[2];
rz(-2.6654976) q[2];
sx q[2];
rz(-0.27175909) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.70499252) q[1];
sx q[1];
rz(-2.3967522) q[1];
sx q[1];
rz(3.0158494) q[1];
rz(1.9022835) q[3];
sx q[3];
rz(-2.2205065) q[3];
sx q[3];
rz(-3.0420764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1364253) q[2];
sx q[2];
rz(-0.74446669) q[2];
sx q[2];
rz(1.0235419) q[2];
rz(-1.3905904) q[3];
sx q[3];
rz(-1.7460881) q[3];
sx q[3];
rz(1.3172147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35109529) q[0];
sx q[0];
rz(-1.0095162) q[0];
sx q[0];
rz(-0.56280953) q[0];
rz(-1.6273181) q[1];
sx q[1];
rz(-1.7104251) q[1];
sx q[1];
rz(-3.0624342) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3252727) q[0];
sx q[0];
rz(-2.5421341) q[0];
sx q[0];
rz(-2.5599203) q[0];
x q[1];
rz(2.6508826) q[2];
sx q[2];
rz(-1.5348002) q[2];
sx q[2];
rz(-1.3513868) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.2217949) q[1];
sx q[1];
rz(-2.256003) q[1];
sx q[1];
rz(-2.5849839) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8757365) q[3];
sx q[3];
rz(-2.1453224) q[3];
sx q[3];
rz(-2.4158626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9074273) q[2];
sx q[2];
rz(-1.5831466) q[2];
sx q[2];
rz(1.8094212) q[2];
rz(2.7857156) q[3];
sx q[3];
rz(-1.9477113) q[3];
sx q[3];
rz(-2.1607384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45977649) q[0];
sx q[0];
rz(-1.1695319) q[0];
sx q[0];
rz(2.0303149) q[0];
rz(-2.2214644) q[1];
sx q[1];
rz(-1.5524813) q[1];
sx q[1];
rz(0.2535325) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0363844) q[0];
sx q[0];
rz(-1.0691133) q[0];
sx q[0];
rz(-1.2060449) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9811834) q[2];
sx q[2];
rz(-1.6840076) q[2];
sx q[2];
rz(-0.28217523) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.92479268) q[1];
sx q[1];
rz(-0.82051045) q[1];
sx q[1];
rz(1.8245068) q[1];
rz(-1.0658137) q[3];
sx q[3];
rz(-1.8441213) q[3];
sx q[3];
rz(-1.2697288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.8186875) q[2];
sx q[2];
rz(-2.5762168) q[2];
sx q[2];
rz(-1.6820071) q[2];
rz(-2.5679576) q[3];
sx q[3];
rz(-1.2021659) q[3];
sx q[3];
rz(0.044376317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6673073) q[0];
sx q[0];
rz(-2.0720338) q[0];
sx q[0];
rz(-1.3473508) q[0];
rz(-2.5509293) q[1];
sx q[1];
rz(-1.2202411) q[1];
sx q[1];
rz(1.4264872) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9085981) q[0];
sx q[0];
rz(-0.54366606) q[0];
sx q[0];
rz(-2.3999016) q[0];
rz(-1.1828184) q[2];
sx q[2];
rz(-2.175594) q[2];
sx q[2];
rz(-1.6590349) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.062078309) q[1];
sx q[1];
rz(-0.2582363) q[1];
sx q[1];
rz(0.18194992) q[1];
x q[2];
rz(-2.7864706) q[3];
sx q[3];
rz(-1.5105845) q[3];
sx q[3];
rz(-0.2015686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.23919375) q[2];
sx q[2];
rz(-1.0603443) q[2];
sx q[2];
rz(-1.1725461) q[2];
rz(-3.0555365) q[3];
sx q[3];
rz(-0.9525758) q[3];
sx q[3];
rz(-2.9854767) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0758783) q[0];
sx q[0];
rz(-1.5471764) q[0];
sx q[0];
rz(2.2789047) q[0];
rz(-2.8942096) q[1];
sx q[1];
rz(-1.8377973) q[1];
sx q[1];
rz(-2.6695796) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.056045) q[0];
sx q[0];
rz(-1.2331095) q[0];
sx q[0];
rz(2.4187536) q[0];
rz(-pi) q[1];
rz(-1.0649933) q[2];
sx q[2];
rz(-1.7031295) q[2];
sx q[2];
rz(-1.1810609) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0942877) q[1];
sx q[1];
rz(-1.7475114) q[1];
sx q[1];
rz(-0.30728886) q[1];
rz(-pi) q[2];
rz(0.94859) q[3];
sx q[3];
rz(-1.3974117) q[3];
sx q[3];
rz(2.9297048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.90999675) q[2];
sx q[2];
rz(-0.96698499) q[2];
sx q[2];
rz(-1.3989353) q[2];
rz(1.4553962) q[3];
sx q[3];
rz(-0.78794909) q[3];
sx q[3];
rz(0.56070352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48651925) q[0];
sx q[0];
rz(-1.1460679) q[0];
sx q[0];
rz(2.6229677) q[0];
rz(-0.71867603) q[1];
sx q[1];
rz(-2.6287754) q[1];
sx q[1];
rz(1.8124883) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5868673) q[0];
sx q[0];
rz(-1.9828078) q[0];
sx q[0];
rz(-2.7464965) q[0];
rz(-1.3955922) q[2];
sx q[2];
rz(-2.8602798) q[2];
sx q[2];
rz(0.11872053) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3997201) q[1];
sx q[1];
rz(-0.73556566) q[1];
sx q[1];
rz(-1.6198141) q[1];
rz(3.1217478) q[3];
sx q[3];
rz(-0.9844508) q[3];
sx q[3];
rz(1.3278409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3507877) q[2];
sx q[2];
rz(-0.69956508) q[2];
sx q[2];
rz(0.038185509) q[2];
rz(2.3112678) q[3];
sx q[3];
rz(-1.3874976) q[3];
sx q[3];
rz(0.047688095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6787978) q[0];
sx q[0];
rz(-0.21634732) q[0];
sx q[0];
rz(1.3702673) q[0];
rz(-0.38617745) q[1];
sx q[1];
rz(-1.736707) q[1];
sx q[1];
rz(0.6699627) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7012335) q[0];
sx q[0];
rz(-2.0645077) q[0];
sx q[0];
rz(0.53892737) q[0];
rz(-pi) q[1];
rz(-0.092548056) q[2];
sx q[2];
rz(-0.71882788) q[2];
sx q[2];
rz(-2.3859442) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0112056) q[1];
sx q[1];
rz(-2.1631259) q[1];
sx q[1];
rz(-0.13641178) q[1];
rz(1.3588328) q[3];
sx q[3];
rz(-1.5342714) q[3];
sx q[3];
rz(-1.9741457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2204444) q[2];
sx q[2];
rz(-0.28382742) q[2];
sx q[2];
rz(2.0849126) q[2];
rz(1.4165261) q[3];
sx q[3];
rz(-1.8959911) q[3];
sx q[3];
rz(-1.2904803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073762745) q[0];
sx q[0];
rz(-1.0404328) q[0];
sx q[0];
rz(2.1347866) q[0];
rz(-2.6489068) q[1];
sx q[1];
rz(-1.410781) q[1];
sx q[1];
rz(2.3115092) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5093436) q[0];
sx q[0];
rz(-1.8496939) q[0];
sx q[0];
rz(-1.5362025) q[0];
rz(-pi) q[1];
rz(0.5987723) q[2];
sx q[2];
rz(-1.3922924) q[2];
sx q[2];
rz(2.0395961) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0841024) q[1];
sx q[1];
rz(-1.7508645) q[1];
sx q[1];
rz(-1.7258641) q[1];
rz(-3.1287635) q[3];
sx q[3];
rz(-2.0137824) q[3];
sx q[3];
rz(-1.6010546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3229708) q[2];
sx q[2];
rz(-1.8399723) q[2];
sx q[2];
rz(2.0446365) q[2];
rz(-1.2777626) q[3];
sx q[3];
rz(-2.4570229) q[3];
sx q[3];
rz(0.6440312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85035664) q[0];
sx q[0];
rz(-2.0391897) q[0];
sx q[0];
rz(2.7665603) q[0];
rz(-3.0229783) q[1];
sx q[1];
rz(-1.7053441) q[1];
sx q[1];
rz(-0.92432252) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7867904) q[0];
sx q[0];
rz(-1.2806007) q[0];
sx q[0];
rz(2.5470887) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0807939) q[2];
sx q[2];
rz(-0.96335232) q[2];
sx q[2];
rz(2.3404672) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.1076775) q[1];
sx q[1];
rz(-2.7795791) q[1];
sx q[1];
rz(3.0605143) q[1];
rz(-pi) q[2];
rz(2.6931346) q[3];
sx q[3];
rz(-1.8844386) q[3];
sx q[3];
rz(0.83066108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2529926) q[2];
sx q[2];
rz(-1.000095) q[2];
sx q[2];
rz(-0.89317733) q[2];
rz(1.1183974) q[3];
sx q[3];
rz(-2.1233386) q[3];
sx q[3];
rz(0.65677381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3746344) q[0];
sx q[0];
rz(-2.0603016) q[0];
sx q[0];
rz(1.7745071) q[0];
rz(-2.959666) q[1];
sx q[1];
rz(-0.55768273) q[1];
sx q[1];
rz(2.2348977) q[1];
rz(-1.295455) q[2];
sx q[2];
rz(-0.32651797) q[2];
sx q[2];
rz(1.4985195) q[2];
rz(-0.652486) q[3];
sx q[3];
rz(-0.62118821) q[3];
sx q[3];
rz(-1.4493829) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
