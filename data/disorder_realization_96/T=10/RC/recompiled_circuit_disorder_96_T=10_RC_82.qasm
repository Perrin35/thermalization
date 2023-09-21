OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.049790073) q[0];
sx q[0];
rz(3.0135305) q[0];
sx q[0];
rz(11.749) q[0];
rz(0.983239) q[1];
sx q[1];
rz(-0.53951889) q[1];
sx q[1];
rz(-1.2004381) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.565633) q[0];
sx q[0];
rz(-0.48086777) q[0];
sx q[0];
rz(-0.54310449) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.13237662) q[2];
sx q[2];
rz(-2.1059603) q[2];
sx q[2];
rz(1.7420499) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9125036) q[1];
sx q[1];
rz(-2.2664245) q[1];
sx q[1];
rz(-1.2234664) q[1];
x q[2];
rz(-2.82253) q[3];
sx q[3];
rz(-0.82108077) q[3];
sx q[3];
rz(-1.0533489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0212705) q[2];
sx q[2];
rz(-2.5734084) q[2];
sx q[2];
rz(-1.5585287) q[2];
rz(0.99672404) q[3];
sx q[3];
rz(-2.6895027) q[3];
sx q[3];
rz(-0.42580095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5581756) q[0];
sx q[0];
rz(-0.75220627) q[0];
sx q[0];
rz(0.054071991) q[0];
rz(1.1955098) q[1];
sx q[1];
rz(-2.1046488) q[1];
sx q[1];
rz(-0.53584677) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.318583) q[0];
sx q[0];
rz(-2.5479925) q[0];
sx q[0];
rz(1.7943322) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4321062) q[2];
sx q[2];
rz(-1.3777133) q[2];
sx q[2];
rz(0.62765861) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0332196) q[1];
sx q[1];
rz(-2.6544826) q[1];
sx q[1];
rz(-0.56652041) q[1];
rz(-pi) q[2];
rz(1.7137394) q[3];
sx q[3];
rz(-1.5385475) q[3];
sx q[3];
rz(1.3168207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1198931) q[2];
sx q[2];
rz(-2.0844441) q[2];
sx q[2];
rz(-2.1832441) q[2];
rz(0.066453233) q[3];
sx q[3];
rz(-1.5575912) q[3];
sx q[3];
rz(2.691793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7217343) q[0];
sx q[0];
rz(-1.2092713) q[0];
sx q[0];
rz(-2.9911175) q[0];
rz(0.45723215) q[1];
sx q[1];
rz(-2.229264) q[1];
sx q[1];
rz(-0.025807468) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8255071) q[0];
sx q[0];
rz(-1.8206017) q[0];
sx q[0];
rz(0.020629701) q[0];
x q[1];
rz(1.2984367) q[2];
sx q[2];
rz(-2.8189427) q[2];
sx q[2];
rz(1.1622365) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7681231) q[1];
sx q[1];
rz(-2.3489531) q[1];
sx q[1];
rz(2.5440689) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13173007) q[3];
sx q[3];
rz(-1.9506491) q[3];
sx q[3];
rz(-0.47117885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0187443) q[2];
sx q[2];
rz(-0.3652502) q[2];
sx q[2];
rz(-0.5853816) q[2];
rz(-0.18150005) q[3];
sx q[3];
rz(-1.8094962) q[3];
sx q[3];
rz(1.5766778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.240775) q[0];
sx q[0];
rz(-2.5188991) q[0];
sx q[0];
rz(0.17661072) q[0];
rz(2.2606842) q[1];
sx q[1];
rz(-1.0662339) q[1];
sx q[1];
rz(2.6054629) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6484084) q[0];
sx q[0];
rz(-0.83183653) q[0];
sx q[0];
rz(2.1302845) q[0];
rz(-2.9049126) q[2];
sx q[2];
rz(-1.5152021) q[2];
sx q[2];
rz(-2.6829164) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0587479) q[1];
sx q[1];
rz(-1.9472329) q[1];
sx q[1];
rz(-2.46545) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7809308) q[3];
sx q[3];
rz(-0.70913991) q[3];
sx q[3];
rz(3.0879471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6716016) q[2];
sx q[2];
rz(-1.4208379) q[2];
sx q[2];
rz(2.0969351) q[2];
rz(-2.4345543) q[3];
sx q[3];
rz(-0.98393327) q[3];
sx q[3];
rz(-0.35693359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9064643) q[0];
sx q[0];
rz(-1.3236073) q[0];
sx q[0];
rz(-2.0902324) q[0];
rz(-1.4936739) q[1];
sx q[1];
rz(-0.58902478) q[1];
sx q[1];
rz(-3.0984745) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.608425) q[0];
sx q[0];
rz(-2.5072917) q[0];
sx q[0];
rz(-1.8932635) q[0];
rz(-pi) q[1];
rz(-1.0518603) q[2];
sx q[2];
rz(-1.8199925) q[2];
sx q[2];
rz(0.11617004) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9777898) q[1];
sx q[1];
rz(-1.8352574) q[1];
sx q[1];
rz(-0.28122854) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14857265) q[3];
sx q[3];
rz(-2.3464977) q[3];
sx q[3];
rz(0.038392301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.12895) q[2];
sx q[2];
rz(-2.6439715) q[2];
sx q[2];
rz(2.7894003) q[2];
rz(-0.59018618) q[3];
sx q[3];
rz(-0.47368172) q[3];
sx q[3];
rz(-0.56110704) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6234289) q[0];
sx q[0];
rz(-1.9264899) q[0];
sx q[0];
rz(1.9859001) q[0];
rz(-2.39134) q[1];
sx q[1];
rz(-2.2022088) q[1];
sx q[1];
rz(2.0828784) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6569865) q[0];
sx q[0];
rz(-1.3625506) q[0];
sx q[0];
rz(-0.44021846) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6201453) q[2];
sx q[2];
rz(-1.2505184) q[2];
sx q[2];
rz(0.59200586) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8844879) q[1];
sx q[1];
rz(-0.37933644) q[1];
sx q[1];
rz(-0.92909716) q[1];
rz(-2.9162354) q[3];
sx q[3];
rz(-2.4212824) q[3];
sx q[3];
rz(2.5845137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.47026149) q[2];
sx q[2];
rz(-1.7241717) q[2];
sx q[2];
rz(1.9227825) q[2];
rz(-1.1550711) q[3];
sx q[3];
rz(-2.9030436) q[3];
sx q[3];
rz(1.4412122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0999775) q[0];
sx q[0];
rz(-1.3168553) q[0];
sx q[0];
rz(-1.51145) q[0];
rz(-1.3776243) q[1];
sx q[1];
rz(-0.310985) q[1];
sx q[1];
rz(2.2999433) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11628843) q[0];
sx q[0];
rz(-1.3447273) q[0];
sx q[0];
rz(0.53209214) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5146671) q[2];
sx q[2];
rz(-1.4142087) q[2];
sx q[2];
rz(-2.2121034) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.48982606) q[1];
sx q[1];
rz(-3.0511599) q[1];
sx q[1];
rz(1.0033146) q[1];
x q[2];
rz(-1.3326725) q[3];
sx q[3];
rz(-2.3861285) q[3];
sx q[3];
rz(1.2681703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1402309) q[2];
sx q[2];
rz(-1.3683687) q[2];
sx q[2];
rz(0.049953071) q[2];
rz(2.4800381) q[3];
sx q[3];
rz(-2.619132) q[3];
sx q[3];
rz(-0.18930999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89896232) q[0];
sx q[0];
rz(-2.1147418) q[0];
sx q[0];
rz(1.7393973) q[0];
rz(0.095480355) q[1];
sx q[1];
rz(-1.9752558) q[1];
sx q[1];
rz(-2.7239674) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0324875) q[0];
sx q[0];
rz(-1.0787449) q[0];
sx q[0];
rz(2.0672654) q[0];
x q[1];
rz(-1.2049963) q[2];
sx q[2];
rz(-1.6779643) q[2];
sx q[2];
rz(1.1366213) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8199181) q[1];
sx q[1];
rz(-0.770861) q[1];
sx q[1];
rz(0.8701156) q[1];
rz(-1.1975708) q[3];
sx q[3];
rz(-1.8456568) q[3];
sx q[3];
rz(1.1653792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.79545704) q[2];
sx q[2];
rz(-1.6980349) q[2];
sx q[2];
rz(2.0987089) q[2];
rz(-2.4677094) q[3];
sx q[3];
rz(-1.6582812) q[3];
sx q[3];
rz(-2.2414482) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3089356) q[0];
sx q[0];
rz(-1.4082264) q[0];
sx q[0];
rz(2.5119264) q[0];
rz(-0.57485238) q[1];
sx q[1];
rz(-1.8300627) q[1];
sx q[1];
rz(0.94690698) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1614721) q[0];
sx q[0];
rz(-1.8959909) q[0];
sx q[0];
rz(1.8878216) q[0];
rz(2.7526555) q[2];
sx q[2];
rz(-0.10483327) q[2];
sx q[2];
rz(2.7742085) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8967594) q[1];
sx q[1];
rz(-1.6549126) q[1];
sx q[1];
rz(0.33478488) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9506748) q[3];
sx q[3];
rz(-2.0058161) q[3];
sx q[3];
rz(-2.2866979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5809014) q[2];
sx q[2];
rz(-1.1297444) q[2];
sx q[2];
rz(-1.2488731) q[2];
rz(0.71436626) q[3];
sx q[3];
rz(-1.8656105) q[3];
sx q[3];
rz(2.8760288) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0749851) q[0];
sx q[0];
rz(-0.62244901) q[0];
sx q[0];
rz(2.4865436) q[0];
rz(-0.89637268) q[1];
sx q[1];
rz(-0.90677774) q[1];
sx q[1];
rz(0.64430976) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2243758) q[0];
sx q[0];
rz(-1.7443568) q[0];
sx q[0];
rz(1.5149084) q[0];
x q[1];
rz(1.8337433) q[2];
sx q[2];
rz(-2.408228) q[2];
sx q[2];
rz(-0.76542379) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6380784) q[1];
sx q[1];
rz(-0.44736171) q[1];
sx q[1];
rz(-0.7034941) q[1];
rz(0.82018567) q[3];
sx q[3];
rz(-2.5327442) q[3];
sx q[3];
rz(-1.5370777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7908988) q[2];
sx q[2];
rz(-1.472298) q[2];
sx q[2];
rz(-2.541686) q[2];
rz(0.89896262) q[3];
sx q[3];
rz(-2.9581684) q[3];
sx q[3];
rz(1.3658587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.29466378) q[0];
sx q[0];
rz(-1.8142721) q[0];
sx q[0];
rz(-0.66692764) q[0];
rz(-2.9121493) q[1];
sx q[1];
rz(-0.89090092) q[1];
sx q[1];
rz(0.13577239) q[1];
rz(0.94639287) q[2];
sx q[2];
rz(-1.4675491) q[2];
sx q[2];
rz(1.8617392) q[2];
rz(-0.89623981) q[3];
sx q[3];
rz(-1.8825718) q[3];
sx q[3];
rz(2.752302) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];