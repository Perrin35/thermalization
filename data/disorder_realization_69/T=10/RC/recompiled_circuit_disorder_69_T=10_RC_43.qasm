OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7944613) q[0];
sx q[0];
rz(-2.1262655) q[0];
sx q[0];
rz(-0.46749687) q[0];
rz(-0.52019083) q[1];
sx q[1];
rz(-1.3462892) q[1];
sx q[1];
rz(2.4612114) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2517393) q[0];
sx q[0];
rz(-1.4407053) q[0];
sx q[0];
rz(2.2962909) q[0];
x q[1];
rz(-0.19619588) q[2];
sx q[2];
rz(-1.8610524) q[2];
sx q[2];
rz(2.6930075) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.89741325) q[1];
sx q[1];
rz(-1.3379828) q[1];
sx q[1];
rz(1.0531055) q[1];
rz(-pi) q[2];
rz(-0.47547961) q[3];
sx q[3];
rz(-0.37215044) q[3];
sx q[3];
rz(1.9843769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.15158571) q[2];
sx q[2];
rz(-0.9881343) q[2];
sx q[2];
rz(-3.0541259) q[2];
rz(-0.72922373) q[3];
sx q[3];
rz(-2.7681523) q[3];
sx q[3];
rz(0.27174404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7049578) q[0];
sx q[0];
rz(-2.3925662) q[0];
sx q[0];
rz(1.0013642) q[0];
rz(0.17240605) q[1];
sx q[1];
rz(-2.0253851) q[1];
sx q[1];
rz(0.52406812) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.24633) q[0];
sx q[0];
rz(-1.8662013) q[0];
sx q[0];
rz(-0.88460478) q[0];
rz(-pi) q[1];
rz(-2.7153154) q[2];
sx q[2];
rz(-0.91765109) q[2];
sx q[2];
rz(0.041235812) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3327649) q[1];
sx q[1];
rz(-1.1953224) q[1];
sx q[1];
rz(-1.3198225) q[1];
rz(-1.7254132) q[3];
sx q[3];
rz(-1.039045) q[3];
sx q[3];
rz(1.0242517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9841763) q[2];
sx q[2];
rz(-2.6817862) q[2];
sx q[2];
rz(1.3400419) q[2];
rz(-0.79483461) q[3];
sx q[3];
rz(-2.0017616) q[3];
sx q[3];
rz(3.1047344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79725093) q[0];
sx q[0];
rz(-2.4448555) q[0];
sx q[0];
rz(-2.5168193) q[0];
rz(2.1773188) q[1];
sx q[1];
rz(-0.48502973) q[1];
sx q[1];
rz(0.18951167) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.097475826) q[0];
sx q[0];
rz(-2.2783845) q[0];
sx q[0];
rz(1.9301027) q[0];
rz(-pi) q[1];
rz(1.6763716) q[2];
sx q[2];
rz(-1.1974679) q[2];
sx q[2];
rz(0.8651274) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7889001) q[1];
sx q[1];
rz(-2.8955724) q[1];
sx q[1];
rz(-2.9109863) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.965851) q[3];
sx q[3];
rz(-2.1648241) q[3];
sx q[3];
rz(2.9807846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0345962) q[2];
sx q[2];
rz(-2.2980289) q[2];
sx q[2];
rz(-1.3872046) q[2];
rz(2.710279) q[3];
sx q[3];
rz(-1.286819) q[3];
sx q[3];
rz(-2.0534024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65961924) q[0];
sx q[0];
rz(-2.1932333) q[0];
sx q[0];
rz(1.5959928) q[0];
rz(2.003147) q[1];
sx q[1];
rz(-0.7754511) q[1];
sx q[1];
rz(1.0901573) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91256234) q[0];
sx q[0];
rz(-1.9711442) q[0];
sx q[0];
rz(1.9489261) q[0];
rz(-pi) q[1];
rz(2.1074739) q[2];
sx q[2];
rz(-0.82674971) q[2];
sx q[2];
rz(0.85988753) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.412147) q[1];
sx q[1];
rz(-0.73819654) q[1];
sx q[1];
rz(-0.63936887) q[1];
x q[2];
rz(0.4269883) q[3];
sx q[3];
rz(-1.3823576) q[3];
sx q[3];
rz(-2.6915336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.99889341) q[2];
sx q[2];
rz(-2.6901851) q[2];
sx q[2];
rz(2.2857655) q[2];
rz(1.9479729) q[3];
sx q[3];
rz(-1.5193628) q[3];
sx q[3];
rz(0.91317552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.49044931) q[0];
sx q[0];
rz(-0.97390807) q[0];
sx q[0];
rz(-0.14973101) q[0];
rz(-2.1504452) q[1];
sx q[1];
rz(-1.2065572) q[1];
sx q[1];
rz(-1.1700464) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8613375) q[0];
sx q[0];
rz(-2.4147408) q[0];
sx q[0];
rz(-2.0656385) q[0];
rz(2.8180426) q[2];
sx q[2];
rz(-0.91634446) q[2];
sx q[2];
rz(2.8138585) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.425881) q[1];
sx q[1];
rz(-0.60957805) q[1];
sx q[1];
rz(-2.4063935) q[1];
x q[2];
rz(0.53797651) q[3];
sx q[3];
rz(-1.7378983) q[3];
sx q[3];
rz(2.3218384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.871792) q[2];
sx q[2];
rz(-1.5912) q[2];
sx q[2];
rz(3.0043547) q[2];
rz(-1.8042701) q[3];
sx q[3];
rz(-0.58115712) q[3];
sx q[3];
rz(-1.8491245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9542434) q[0];
sx q[0];
rz(-1.3784778) q[0];
sx q[0];
rz(-1.7911918) q[0];
rz(-2.2987135) q[1];
sx q[1];
rz(-0.73892361) q[1];
sx q[1];
rz(2.4687016) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0762754) q[0];
sx q[0];
rz(-0.99940171) q[0];
sx q[0];
rz(0.12302834) q[0];
x q[1];
rz(1.5422836) q[2];
sx q[2];
rz(-0.71454222) q[2];
sx q[2];
rz(-1.5039832) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.91025483) q[1];
sx q[1];
rz(-0.64097039) q[1];
sx q[1];
rz(2.7892968) q[1];
rz(1.4937431) q[3];
sx q[3];
rz(-1.1601845) q[3];
sx q[3];
rz(-1.3239469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.792753) q[2];
sx q[2];
rz(-1.9323843) q[2];
sx q[2];
rz(-2.7065281) q[2];
rz(-1.7815636) q[3];
sx q[3];
rz(-0.74917787) q[3];
sx q[3];
rz(-0.24766651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.172794) q[0];
sx q[0];
rz(-0.89247576) q[0];
sx q[0];
rz(-2.0794179) q[0];
rz(-2.0299714) q[1];
sx q[1];
rz(-1.9042791) q[1];
sx q[1];
rz(1.7395082) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8334478) q[0];
sx q[0];
rz(-1.2045367) q[0];
sx q[0];
rz(-1.2783947) q[0];
rz(-pi) q[1];
rz(-1.3349322) q[2];
sx q[2];
rz(-2.3216322) q[2];
sx q[2];
rz(1.6279398) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.40560383) q[1];
sx q[1];
rz(-0.39499184) q[1];
sx q[1];
rz(0.40785457) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13774638) q[3];
sx q[3];
rz(-0.98689729) q[3];
sx q[3];
rz(0.13970845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7193675) q[2];
sx q[2];
rz(-0.56754595) q[2];
sx q[2];
rz(-0.68022234) q[2];
rz(0.42823544) q[3];
sx q[3];
rz(-1.2546344) q[3];
sx q[3];
rz(-1.7817106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5087886) q[0];
sx q[0];
rz(-1.7254242) q[0];
sx q[0];
rz(1.592214) q[0];
rz(-2.8920065) q[1];
sx q[1];
rz(-1.9892178) q[1];
sx q[1];
rz(-2.6002398) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4562948) q[0];
sx q[0];
rz(-2.1349499) q[0];
sx q[0];
rz(-1.5840522) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45221046) q[2];
sx q[2];
rz(-1.2621677) q[2];
sx q[2];
rz(-1.897915) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7975446) q[1];
sx q[1];
rz(-2.27564) q[1];
sx q[1];
rz(-1.0182347) q[1];
x q[2];
rz(-0.3803216) q[3];
sx q[3];
rz(-1.5690648) q[3];
sx q[3];
rz(-2.1960432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.044518746) q[2];
sx q[2];
rz(-2.7842583) q[2];
sx q[2];
rz(-2.3642335) q[2];
rz(-2.2903531) q[3];
sx q[3];
rz(-2.080353) q[3];
sx q[3];
rz(-1.8204934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38354307) q[0];
sx q[0];
rz(-1.3306916) q[0];
sx q[0];
rz(0.0099649075) q[0];
rz(1.0154356) q[1];
sx q[1];
rz(-2.3773057) q[1];
sx q[1];
rz(1.6962956) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9706668) q[0];
sx q[0];
rz(-0.88060856) q[0];
sx q[0];
rz(2.4404581) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2465835) q[2];
sx q[2];
rz(-1.0969321) q[2];
sx q[2];
rz(-0.5627788) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0547486) q[1];
sx q[1];
rz(-2.5111685) q[1];
sx q[1];
rz(-3.1151506) q[1];
x q[2];
rz(0.61049283) q[3];
sx q[3];
rz(-2.1497512) q[3];
sx q[3];
rz(-2.0905153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4426667) q[2];
sx q[2];
rz(-2.2103504) q[2];
sx q[2];
rz(2.5342069) q[2];
rz(1.4390885) q[3];
sx q[3];
rz(-1.3834229) q[3];
sx q[3];
rz(-2.3506892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8059175) q[0];
sx q[0];
rz(-1.4842002) q[0];
sx q[0];
rz(2.5277396) q[0];
rz(2.0954258) q[1];
sx q[1];
rz(-0.26509735) q[1];
sx q[1];
rz(0.39224958) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0948254) q[0];
sx q[0];
rz(-1.9691159) q[0];
sx q[0];
rz(0.70478583) q[0];
x q[1];
rz(2.4105083) q[2];
sx q[2];
rz(-0.77027551) q[2];
sx q[2];
rz(-2.0172271) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2892462) q[1];
sx q[1];
rz(-0.95818555) q[1];
sx q[1];
rz(0.5457408) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13855374) q[3];
sx q[3];
rz(-2.9346653) q[3];
sx q[3];
rz(-3.0272527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.067651) q[2];
sx q[2];
rz(-1.3839046) q[2];
sx q[2];
rz(2.5411141) q[2];
rz(1.0673808) q[3];
sx q[3];
rz(-1.821527) q[3];
sx q[3];
rz(-1.0935121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7077211) q[0];
sx q[0];
rz(-2.7086471) q[0];
sx q[0];
rz(1.4824296) q[0];
rz(-2.1451163) q[1];
sx q[1];
rz(-1.4947718) q[1];
sx q[1];
rz(1.6047118) q[1];
rz(-2.773182) q[2];
sx q[2];
rz(-0.80130063) q[2];
sx q[2];
rz(-0.038566312) q[2];
rz(2.2511803) q[3];
sx q[3];
rz(-2.012841) q[3];
sx q[3];
rz(-2.8750318) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
