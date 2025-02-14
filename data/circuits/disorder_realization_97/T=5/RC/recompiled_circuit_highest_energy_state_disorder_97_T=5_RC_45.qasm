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
rz(1.3069557) q[0];
sx q[0];
rz(-2.262158) q[0];
sx q[0];
rz(0.49361324) q[0];
rz(1.8618795) q[1];
sx q[1];
rz(-0.6266098) q[1];
sx q[1];
rz(-1.6630747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8956229) q[0];
sx q[0];
rz(-0.9446656) q[0];
sx q[0];
rz(-0.62762053) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.25716146) q[2];
sx q[2];
rz(-0.92338054) q[2];
sx q[2];
rz(-2.0055298) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.93881685) q[1];
sx q[1];
rz(-2.4033053) q[1];
sx q[1];
rz(0.32719214) q[1];
rz(-pi) q[2];
rz(0.36840393) q[3];
sx q[3];
rz(-2.1737636) q[3];
sx q[3];
rz(1.4908229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0066068) q[2];
sx q[2];
rz(-2.0015494) q[2];
sx q[2];
rz(-1.1084278) q[2];
rz(1.3125575) q[3];
sx q[3];
rz(-0.24459845) q[3];
sx q[3];
rz(0.31497064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38158622) q[0];
sx q[0];
rz(-0.8257603) q[0];
sx q[0];
rz(3.1257358) q[0];
rz(2.1511757) q[1];
sx q[1];
rz(-2.3742193) q[1];
sx q[1];
rz(2.2784065) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28672934) q[0];
sx q[0];
rz(-0.83896381) q[0];
sx q[0];
rz(1.7006204) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3349971) q[2];
sx q[2];
rz(-1.2588698) q[2];
sx q[2];
rz(-0.27659135) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9941123) q[1];
sx q[1];
rz(-2.0614455) q[1];
sx q[1];
rz(1.3730241) q[1];
x q[2];
rz(-2.9921069) q[3];
sx q[3];
rz(-0.85452291) q[3];
sx q[3];
rz(0.0090741875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4414759) q[2];
sx q[2];
rz(-2.8642544) q[2];
sx q[2];
rz(-0.46879834) q[2];
rz(-0.68764728) q[3];
sx q[3];
rz(-1.629849) q[3];
sx q[3];
rz(-0.30011737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1402635) q[0];
sx q[0];
rz(-0.92324531) q[0];
sx q[0];
rz(-0.73626751) q[0];
rz(2.7983792) q[1];
sx q[1];
rz(-0.88756573) q[1];
sx q[1];
rz(-2.9578178) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60584468) q[0];
sx q[0];
rz(-0.82855201) q[0];
sx q[0];
rz(-1.6607619) q[0];
rz(1.0509346) q[2];
sx q[2];
rz(-0.94280548) q[2];
sx q[2];
rz(1.7597212) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.80189542) q[1];
sx q[1];
rz(-1.2725127) q[1];
sx q[1];
rz(2.4043179) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.36005693) q[3];
sx q[3];
rz(-1.0601794) q[3];
sx q[3];
rz(-2.1329481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3595714) q[2];
sx q[2];
rz(-0.51681334) q[2];
sx q[2];
rz(2.6256631) q[2];
rz(-1.1082209) q[3];
sx q[3];
rz(-0.55086946) q[3];
sx q[3];
rz(1.1512604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16294031) q[0];
sx q[0];
rz(-1.9300224) q[0];
sx q[0];
rz(2.6287855) q[0];
rz(2.6817952) q[1];
sx q[1];
rz(-1.5492946) q[1];
sx q[1];
rz(2.3777681) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0054454) q[0];
sx q[0];
rz(-0.17373611) q[0];
sx q[0];
rz(-0.72189118) q[0];
rz(-pi) q[1];
rz(-1.2470719) q[2];
sx q[2];
rz(-1.1737595) q[2];
sx q[2];
rz(0.70003245) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3556174) q[1];
sx q[1];
rz(-1.2316868) q[1];
sx q[1];
rz(-0.29735844) q[1];
rz(0.34448907) q[3];
sx q[3];
rz(-1.4360582) q[3];
sx q[3];
rz(-1.3400934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1226471) q[2];
sx q[2];
rz(-2.985869) q[2];
sx q[2];
rz(0.32749185) q[2];
rz(-2.859595) q[3];
sx q[3];
rz(-1.3982541) q[3];
sx q[3];
rz(-1.5544844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17664385) q[0];
sx q[0];
rz(-2.7556941) q[0];
sx q[0];
rz(-0.97736812) q[0];
rz(-0.36879677) q[1];
sx q[1];
rz(-2.2803523) q[1];
sx q[1];
rz(1.7288953) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57713014) q[0];
sx q[0];
rz(-0.53097314) q[0];
sx q[0];
rz(-1.1891548) q[0];
rz(-pi) q[1];
rz(-3.1277282) q[2];
sx q[2];
rz(-1.0540451) q[2];
sx q[2];
rz(-2.2114829) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2569132) q[1];
sx q[1];
rz(-1.2882782) q[1];
sx q[1];
rz(0.032508548) q[1];
rz(-0.47980185) q[3];
sx q[3];
rz(-0.63463026) q[3];
sx q[3];
rz(3.0632085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8480924) q[2];
sx q[2];
rz(-2.3257181) q[2];
sx q[2];
rz(-2.9917742) q[2];
rz(0.39854974) q[3];
sx q[3];
rz(-1.5236676) q[3];
sx q[3];
rz(-2.985305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3715816) q[0];
sx q[0];
rz(-0.12938975) q[0];
sx q[0];
rz(2.5883801) q[0];
rz(-2.5999056) q[1];
sx q[1];
rz(-1.0463511) q[1];
sx q[1];
rz(-0.4479301) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.124156) q[0];
sx q[0];
rz(-1.4261591) q[0];
sx q[0];
rz(0.38993024) q[0];
x q[1];
rz(2.7121353) q[2];
sx q[2];
rz(-1.1633368) q[2];
sx q[2];
rz(-2.7384659) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.74400157) q[1];
sx q[1];
rz(-2.7251456) q[1];
sx q[1];
rz(-1.9201502) q[1];
x q[2];
rz(-2.8087696) q[3];
sx q[3];
rz(-1.6245922) q[3];
sx q[3];
rz(1.4383565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.64122671) q[2];
sx q[2];
rz(-1.1643103) q[2];
sx q[2];
rz(-0.79916239) q[2];
rz(-0.3847807) q[3];
sx q[3];
rz(-0.32618263) q[3];
sx q[3];
rz(-1.0001812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.7810829) q[0];
sx q[0];
rz(-1.0436844) q[0];
sx q[0];
rz(2.5497896) q[0];
rz(-0.38161713) q[1];
sx q[1];
rz(-2.7615669) q[1];
sx q[1];
rz(0.15242481) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8417609) q[0];
sx q[0];
rz(-2.526847) q[0];
sx q[0];
rz(1.0157981) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3954787) q[2];
sx q[2];
rz(-1.0994871) q[2];
sx q[2];
rz(2.5702554) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.7970711) q[1];
sx q[1];
rz(-2.7943369) q[1];
sx q[1];
rz(3.0879435) q[1];
rz(-2.2493208) q[3];
sx q[3];
rz(-2.0990058) q[3];
sx q[3];
rz(-1.5367374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3064208) q[2];
sx q[2];
rz(-2.1495337) q[2];
sx q[2];
rz(-1.6888118) q[2];
rz(2.6460904) q[3];
sx q[3];
rz(-2.317954) q[3];
sx q[3];
rz(2.5775094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38967663) q[0];
sx q[0];
rz(-0.037411995) q[0];
sx q[0];
rz(2.9801242) q[0];
rz(3.1096733) q[1];
sx q[1];
rz(-0.64161623) q[1];
sx q[1];
rz(-1.8643103) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.204435) q[0];
sx q[0];
rz(-1.3124518) q[0];
sx q[0];
rz(-2.4169371) q[0];
x q[1];
rz(1.1487537) q[2];
sx q[2];
rz(-1.3598056) q[2];
sx q[2];
rz(2.8659897) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.016170597) q[1];
sx q[1];
rz(-1.6787663) q[1];
sx q[1];
rz(2.8362464) q[1];
x q[2];
rz(-1.2767918) q[3];
sx q[3];
rz(-1.6075896) q[3];
sx q[3];
rz(-2.5599673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6892467) q[2];
sx q[2];
rz(-2.2301058) q[2];
sx q[2];
rz(-0.39608836) q[2];
rz(0.39997697) q[3];
sx q[3];
rz(-0.59498274) q[3];
sx q[3];
rz(-2.2649435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-0.24092995) q[0];
sx q[0];
rz(-0.018095896) q[0];
sx q[0];
rz(0.15116365) q[0];
rz(2.1647272) q[1];
sx q[1];
rz(-0.38115373) q[1];
sx q[1];
rz(0.74473286) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.789294) q[0];
sx q[0];
rz(-1.2937102) q[0];
sx q[0];
rz(-0.21488551) q[0];
rz(1.9836203) q[2];
sx q[2];
rz(-1.7016439) q[2];
sx q[2];
rz(-1.3589588) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3724461) q[1];
sx q[1];
rz(-0.92248864) q[1];
sx q[1];
rz(1.5689956) q[1];
rz(2.7895097) q[3];
sx q[3];
rz(-1.3439281) q[3];
sx q[3];
rz(0.20108237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.62290827) q[2];
sx q[2];
rz(-1.9025977) q[2];
sx q[2];
rz(-2.1935479) q[2];
rz(1.0575804) q[3];
sx q[3];
rz(-1.4761997) q[3];
sx q[3];
rz(-1.074056) q[3];
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
rz(pi/2) q[0];
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
rz(-0.093216151) q[0];
sx q[0];
rz(-0.26234782) q[0];
sx q[0];
rz(-2.9076305) q[0];
rz(-1.9562862) q[1];
sx q[1];
rz(-2.2133841) q[1];
sx q[1];
rz(1.2497466) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5025185) q[0];
sx q[0];
rz(-2.5022129) q[0];
sx q[0];
rz(2.6431497) q[0];
x q[1];
rz(1.363867) q[2];
sx q[2];
rz(-2.1812058) q[2];
sx q[2];
rz(2.0774942) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0873336) q[1];
sx q[1];
rz(-2.1076084) q[1];
sx q[1];
rz(2.3786503) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4602182) q[3];
sx q[3];
rz(-1.7077966) q[3];
sx q[3];
rz(-0.43230012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.11448161) q[2];
sx q[2];
rz(-1.1610843) q[2];
sx q[2];
rz(1.7005881) q[2];
rz(-0.13127413) q[3];
sx q[3];
rz(-0.461853) q[3];
sx q[3];
rz(0.24391267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.092939) q[0];
sx q[0];
rz(-0.96554148) q[0];
sx q[0];
rz(1.1337793) q[0];
rz(2.1009905) q[1];
sx q[1];
rz(-2.2941209) q[1];
sx q[1];
rz(2.6170731) q[1];
rz(0.30529373) q[2];
sx q[2];
rz(-2.1758781) q[2];
sx q[2];
rz(0.97347409) q[2];
rz(-1.4199938) q[3];
sx q[3];
rz(-0.94298132) q[3];
sx q[3];
rz(1.0366221) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
