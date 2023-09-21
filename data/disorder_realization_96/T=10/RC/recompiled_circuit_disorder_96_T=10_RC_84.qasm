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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49657208) q[0];
sx q[0];
rz(-1.3294157) q[0];
sx q[0];
rz(2.7215331) q[0];
x q[1];
rz(-3.009216) q[2];
sx q[2];
rz(-2.1059603) q[2];
sx q[2];
rz(-1.7420499) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5696722) q[1];
sx q[1];
rz(-1.306428) q[1];
sx q[1];
rz(-0.7260679) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3463763) q[3];
sx q[3];
rz(-1.3391558) q[3];
sx q[3];
rz(-0.29602805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0212705) q[2];
sx q[2];
rz(-0.56818429) q[2];
sx q[2];
rz(-1.583064) q[2];
rz(2.1448686) q[3];
sx q[3];
rz(-0.45209) q[3];
sx q[3];
rz(2.7157917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5834171) q[0];
sx q[0];
rz(-0.75220627) q[0];
sx q[0];
rz(-0.054071991) q[0];
rz(1.9460829) q[1];
sx q[1];
rz(-2.1046488) q[1];
sx q[1];
rz(-2.6057459) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0756404) q[0];
sx q[0];
rz(-1.4464805) q[0];
sx q[0];
rz(0.9888222) q[0];
x q[1];
rz(0.70948647) q[2];
sx q[2];
rz(-1.7638793) q[2];
sx q[2];
rz(-0.62765861) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.10837308) q[1];
sx q[1];
rz(-0.48711005) q[1];
sx q[1];
rz(-2.5750722) q[1];
x q[2];
rz(-1.7934947) q[3];
sx q[3];
rz(-0.14651146) q[3];
sx q[3];
rz(0.033586249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1198931) q[2];
sx q[2];
rz(-2.0844441) q[2];
sx q[2];
rz(-0.95834857) q[2];
rz(3.0751394) q[3];
sx q[3];
rz(-1.5840014) q[3];
sx q[3];
rz(2.691793) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7217343) q[0];
sx q[0];
rz(-1.2092713) q[0];
sx q[0];
rz(-0.15047519) q[0];
rz(2.6843605) q[1];
sx q[1];
rz(-0.91232863) q[1];
sx q[1];
rz(3.1157852) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8255071) q[0];
sx q[0];
rz(-1.320991) q[0];
sx q[0];
rz(0.020629701) q[0];
x q[1];
rz(-1.2984367) q[2];
sx q[2];
rz(-0.32264999) q[2];
sx q[2];
rz(-1.9793561) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3734696) q[1];
sx q[1];
rz(-0.79263955) q[1];
sx q[1];
rz(0.5975238) q[1];
rz(-0.13173007) q[3];
sx q[3];
rz(-1.1909435) q[3];
sx q[3];
rz(-2.6704138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1228483) q[2];
sx q[2];
rz(-0.3652502) q[2];
sx q[2];
rz(2.5562111) q[2];
rz(-2.9600926) q[3];
sx q[3];
rz(-1.8094962) q[3];
sx q[3];
rz(1.5649149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90081763) q[0];
sx q[0];
rz(-2.5188991) q[0];
sx q[0];
rz(-0.17661072) q[0];
rz(-0.88090849) q[1];
sx q[1];
rz(-2.0753588) q[1];
sx q[1];
rz(-2.6054629) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4767544) q[0];
sx q[0];
rz(-1.1676482) q[0];
sx q[0];
rz(2.320015) q[0];
x q[1];
rz(0.23668004) q[2];
sx q[2];
rz(-1.6263905) q[2];
sx q[2];
rz(2.6829164) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.94169468) q[1];
sx q[1];
rz(-2.3823793) q[1];
sx q[1];
rz(2.5782176) q[1];
rz(-0.36066182) q[3];
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
rz(2.7846591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9064643) q[0];
sx q[0];
rz(-1.8179853) q[0];
sx q[0];
rz(1.0513603) q[0];
rz(-1.4936739) q[1];
sx q[1];
rz(-2.5525679) q[1];
sx q[1];
rz(3.0984745) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9163141) q[0];
sx q[0];
rz(-1.7597223) q[0];
sx q[0];
rz(2.180045) q[0];
x q[1];
rz(0.28508614) q[2];
sx q[2];
rz(-2.0721772) q[2];
sx q[2];
rz(1.3146871) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4700714) q[1];
sx q[1];
rz(-0.38362353) q[1];
sx q[1];
rz(-0.77312153) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78955663) q[3];
sx q[3];
rz(-1.4649179) q[3];
sx q[3];
rz(-1.7136128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.12895) q[2];
sx q[2];
rz(-0.49762112) q[2];
sx q[2];
rz(0.35219231) q[2];
rz(0.59018618) q[3];
sx q[3];
rz(-0.47368172) q[3];
sx q[3];
rz(-2.5804856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6234289) q[0];
sx q[0];
rz(-1.9264899) q[0];
sx q[0];
rz(-1.1556926) q[0];
rz(-0.75025264) q[1];
sx q[1];
rz(-0.9393839) q[1];
sx q[1];
rz(-1.0587143) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6569865) q[0];
sx q[0];
rz(-1.779042) q[0];
sx q[0];
rz(0.44021846) q[0];
rz(-pi) q[1];
rz(-1.2054339) q[2];
sx q[2];
rz(-1.0783256) q[2];
sx q[2];
rz(-1.9838711) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7212972) q[1];
sx q[1];
rz(-1.2696206) q[1];
sx q[1];
rz(2.9073614) q[1];
rz(-pi) q[2];
rz(1.7644464) q[3];
sx q[3];
rz(-2.2691257) q[3];
sx q[3];
rz(-2.8805672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6713312) q[2];
sx q[2];
rz(-1.417421) q[2];
sx q[2];
rz(1.2188101) q[2];
rz(-1.9865215) q[3];
sx q[3];
rz(-2.9030436) q[3];
sx q[3];
rz(1.7003805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0999775) q[0];
sx q[0];
rz(-1.3168553) q[0];
sx q[0];
rz(1.6301427) q[0];
rz(-1.7639683) q[1];
sx q[1];
rz(-2.8306077) q[1];
sx q[1];
rz(2.2999433) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3232988) q[0];
sx q[0];
rz(-0.57384402) q[0];
sx q[0];
rz(0.42563514) q[0];
x q[1];
rz(-1.3782578) q[2];
sx q[2];
rz(-0.95270573) q[2];
sx q[2];
rz(-2.3877909) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6517666) q[1];
sx q[1];
rz(-0.090432743) q[1];
sx q[1];
rz(1.0033146) q[1];
rz(-pi) q[2];
rz(0.21861403) q[3];
sx q[3];
rz(-2.300005) q[3];
sx q[3];
rz(-1.5900172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1402309) q[2];
sx q[2];
rz(-1.3683687) q[2];
sx q[2];
rz(3.0916396) q[2];
rz(0.66155457) q[3];
sx q[3];
rz(-0.52246061) q[3];
sx q[3];
rz(-0.18930999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-1.4021953) q[0];
rz(-3.0461123) q[1];
sx q[1];
rz(-1.9752558) q[1];
sx q[1];
rz(0.41762525) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8861038) q[0];
sx q[0];
rz(-0.68414738) q[0];
sx q[0];
rz(2.4151001) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8629486) q[2];
sx q[2];
rz(-2.7610965) q[2];
sx q[2];
rz(0.7064864) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0954674) q[1];
sx q[1];
rz(-1.0091262) q[1];
sx q[1];
rz(-0.559505) q[1];
rz(-pi) q[2];
rz(1.9440218) q[3];
sx q[3];
rz(-1.2959359) q[3];
sx q[3];
rz(-1.1653792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.79545704) q[2];
sx q[2];
rz(-1.4435578) q[2];
sx q[2];
rz(1.0428838) q[2];
rz(-0.67388326) q[3];
sx q[3];
rz(-1.4833114) q[3];
sx q[3];
rz(0.90014443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3089356) q[0];
sx q[0];
rz(-1.4082264) q[0];
sx q[0];
rz(-2.5119264) q[0];
rz(-0.57485238) q[1];
sx q[1];
rz(-1.8300627) q[1];
sx q[1];
rz(-2.1946857) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1614721) q[0];
sx q[0];
rz(-1.8959909) q[0];
sx q[0];
rz(1.253771) q[0];
rz(-1.6106748) q[2];
sx q[2];
rz(-1.4738184) q[2];
sx q[2];
rz(-0.023488451) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0526035) q[1];
sx q[1];
rz(-2.7967884) q[1];
sx q[1];
rz(2.8903972) q[1];
rz(-2.9506748) q[3];
sx q[3];
rz(-2.0058161) q[3];
sx q[3];
rz(0.85489475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.56069121) q[2];
sx q[2];
rz(-2.0118482) q[2];
sx q[2];
rz(1.2488731) q[2];
rz(0.71436626) q[3];
sx q[3];
rz(-1.2759821) q[3];
sx q[3];
rz(0.26556382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0666075) q[0];
sx q[0];
rz(-2.5191436) q[0];
sx q[0];
rz(-0.65504909) q[0];
rz(-0.89637268) q[1];
sx q[1];
rz(-2.2348149) q[1];
sx q[1];
rz(-0.64430976) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9172168) q[0];
sx q[0];
rz(-1.7443568) q[0];
sx q[0];
rz(-1.5149084) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9115453) q[2];
sx q[2];
rz(-2.2736079) q[2];
sx q[2];
rz(0.41781296) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5035142) q[1];
sx q[1];
rz(-0.44736171) q[1];
sx q[1];
rz(-0.7034941) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82018567) q[3];
sx q[3];
rz(-0.60884848) q[3];
sx q[3];
rz(1.604515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3506938) q[2];
sx q[2];
rz(-1.6692946) q[2];
sx q[2];
rz(2.541686) q[2];
rz(0.89896262) q[3];
sx q[3];
rz(-2.9581684) q[3];
sx q[3];
rz(-1.775734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8469289) q[0];
sx q[0];
rz(-1.3273205) q[0];
sx q[0];
rz(2.474665) q[0];
rz(-0.22944336) q[1];
sx q[1];
rz(-2.2506917) q[1];
sx q[1];
rz(-3.0058203) q[1];
rz(-0.12702604) q[2];
sx q[2];
rz(-2.1913678) q[2];
sx q[2];
rz(0.36507228) q[2];
rz(1.0944081) q[3];
sx q[3];
rz(-0.73275685) q[3];
sx q[3];
rz(-2.3263596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
