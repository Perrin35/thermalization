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
rz(0.82290736) q[0];
sx q[0];
rz(-0.35879254) q[0];
sx q[0];
rz(-2.2732777) q[0];
rz(1.7262285) q[1];
sx q[1];
rz(-2.0077029) q[1];
sx q[1];
rz(0.99376065) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.438765) q[0];
sx q[0];
rz(-1.1701705) q[0];
sx q[0];
rz(1.4482657) q[0];
x q[1];
rz(2.4325718) q[2];
sx q[2];
rz(-2.5160501) q[2];
sx q[2];
rz(2.4838205) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.43683166) q[1];
sx q[1];
rz(-2.0109573) q[1];
sx q[1];
rz(-2.6635567) q[1];
rz(-0.24561974) q[3];
sx q[3];
rz(-1.1455451) q[3];
sx q[3];
rz(-2.5568145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.913784) q[2];
sx q[2];
rz(-1.8081534) q[2];
sx q[2];
rz(0.58161962) q[2];
rz(0.73451129) q[3];
sx q[3];
rz(-1.651265) q[3];
sx q[3];
rz(-0.41828004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-0.89473474) q[0];
sx q[0];
rz(-1.3792091) q[0];
sx q[0];
rz(0.49767622) q[0];
rz(2.1007762) q[1];
sx q[1];
rz(-2.7383995) q[1];
sx q[1];
rz(1.9042447) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029631296) q[0];
sx q[0];
rz(-0.8489767) q[0];
sx q[0];
rz(1.5143751) q[0];
x q[1];
rz(2.9534229) q[2];
sx q[2];
rz(-1.9261179) q[2];
sx q[2];
rz(-2.4057092) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.50857568) q[1];
sx q[1];
rz(-2.687722) q[1];
sx q[1];
rz(1.8967486) q[1];
rz(-pi) q[2];
rz(-2.8760404) q[3];
sx q[3];
rz(-1.7803811) q[3];
sx q[3];
rz(-0.095106212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.70025468) q[2];
sx q[2];
rz(-0.73107084) q[2];
sx q[2];
rz(2.8311484) q[2];
rz(1.9849518) q[3];
sx q[3];
rz(-1.1247331) q[3];
sx q[3];
rz(-1.9748851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0079086) q[0];
sx q[0];
rz(-0.5846566) q[0];
sx q[0];
rz(0.34580082) q[0];
rz(-2.8969104) q[1];
sx q[1];
rz(-2.209765) q[1];
sx q[1];
rz(-1.7049047) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37852851) q[0];
sx q[0];
rz(-0.53685704) q[0];
sx q[0];
rz(-0.92354639) q[0];
rz(-2.56836) q[2];
sx q[2];
rz(-2.4893508) q[2];
sx q[2];
rz(-2.531372) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.50592283) q[1];
sx q[1];
rz(-1.3590006) q[1];
sx q[1];
rz(2.3426272) q[1];
rz(0.64880649) q[3];
sx q[3];
rz(-2.5655364) q[3];
sx q[3];
rz(0.77133178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1338542) q[2];
sx q[2];
rz(-2.014092) q[2];
sx q[2];
rz(-0.63068843) q[2];
rz(0.33356365) q[3];
sx q[3];
rz(-1.0015229) q[3];
sx q[3];
rz(-1.001531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0665322) q[0];
sx q[0];
rz(-0.94828951) q[0];
sx q[0];
rz(-1.7370268) q[0];
rz(2.7724077) q[1];
sx q[1];
rz(-1.4063947) q[1];
sx q[1];
rz(0.085748347) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7333711) q[0];
sx q[0];
rz(-1.2722172) q[0];
sx q[0];
rz(2.7135506) q[0];
rz(-pi) q[1];
rz(-2.4253885) q[2];
sx q[2];
rz(-0.58164454) q[2];
sx q[2];
rz(3.1379791) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8814739) q[1];
sx q[1];
rz(-2.6131995) q[1];
sx q[1];
rz(2.0563988) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8885918) q[3];
sx q[3];
rz(-1.9488584) q[3];
sx q[3];
rz(-1.8658569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3021476) q[2];
sx q[2];
rz(-1.9044694) q[2];
sx q[2];
rz(2.6117924) q[2];
rz(-0.038330404) q[3];
sx q[3];
rz(-2.4119792) q[3];
sx q[3];
rz(-1.5995601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.9180651) q[0];
sx q[0];
rz(-2.6730838) q[0];
sx q[0];
rz(1.9388306) q[0];
rz(-0.27944061) q[1];
sx q[1];
rz(-1.025082) q[1];
sx q[1];
rz(0.7739982) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27085486) q[0];
sx q[0];
rz(-1.2574728) q[0];
sx q[0];
rz(2.3148651) q[0];
rz(-pi) q[1];
rz(-0.33542893) q[2];
sx q[2];
rz(-0.24615363) q[2];
sx q[2];
rz(0.85573643) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7420089) q[1];
sx q[1];
rz(-3.0019433) q[1];
sx q[1];
rz(1.4781471) q[1];
rz(1.0353885) q[3];
sx q[3];
rz(-1.5683953) q[3];
sx q[3];
rz(-0.46214275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1191795) q[2];
sx q[2];
rz(-3.000562) q[2];
sx q[2];
rz(-0.5640344) q[2];
rz(2.4885079) q[3];
sx q[3];
rz(-2.0061195) q[3];
sx q[3];
rz(0.45026067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.7804724) q[0];
sx q[0];
rz(-2.5887964) q[0];
sx q[0];
rz(0.64315382) q[0];
rz(1.9505352) q[1];
sx q[1];
rz(-1.6845614) q[1];
sx q[1];
rz(2.4868884) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9438671) q[0];
sx q[0];
rz(-1.4023313) q[0];
sx q[0];
rz(-1.6351624) q[0];
rz(1.5327318) q[2];
sx q[2];
rz(-0.88721472) q[2];
sx q[2];
rz(2.3823007) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1015849) q[1];
sx q[1];
rz(-1.6363278) q[1];
sx q[1];
rz(2.1145909) q[1];
x q[2];
rz(-1.9611465) q[3];
sx q[3];
rz(-1.6359513) q[3];
sx q[3];
rz(-0.81126838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5074629) q[2];
sx q[2];
rz(-2.7644988) q[2];
sx q[2];
rz(1.6667574) q[2];
rz(-0.48464388) q[3];
sx q[3];
rz(-2.1438997) q[3];
sx q[3];
rz(-1.8528329) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5019048) q[0];
sx q[0];
rz(-2.8785093) q[0];
sx q[0];
rz(2.4581773) q[0];
rz(3.0112093) q[1];
sx q[1];
rz(-1.5905292) q[1];
sx q[1];
rz(-3.108976) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5033753) q[0];
sx q[0];
rz(-2.1299358) q[0];
sx q[0];
rz(0.20737623) q[0];
rz(-0.24222272) q[2];
sx q[2];
rz(-0.70021473) q[2];
sx q[2];
rz(-1.2716573) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.72717818) q[1];
sx q[1];
rz(-1.7373996) q[1];
sx q[1];
rz(-1.3018621) q[1];
x q[2];
rz(-1.873305) q[3];
sx q[3];
rz(-2.4166346) q[3];
sx q[3];
rz(-2.8444949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.42328295) q[2];
sx q[2];
rz(-1.1085359) q[2];
sx q[2];
rz(-0.56524593) q[2];
rz(-1.7806753) q[3];
sx q[3];
rz(-2.8748685) q[3];
sx q[3];
rz(-2.3274073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77124202) q[0];
sx q[0];
rz(-0.057436198) q[0];
sx q[0];
rz(-2.8420319) q[0];
rz(1.6948505) q[1];
sx q[1];
rz(-2.2694777) q[1];
sx q[1];
rz(-2.1902693) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2593237) q[0];
sx q[0];
rz(-0.20765064) q[0];
sx q[0];
rz(1.728968) q[0];
x q[1];
rz(1.6137684) q[2];
sx q[2];
rz(-2.9298721) q[2];
sx q[2];
rz(-0.73390244) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.658706) q[1];
sx q[1];
rz(-1.1466007) q[1];
sx q[1];
rz(1.5931604) q[1];
rz(0.087835066) q[3];
sx q[3];
rz(-1.3099652) q[3];
sx q[3];
rz(-1.5276599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.99913725) q[2];
sx q[2];
rz(-0.74644011) q[2];
sx q[2];
rz(2.8738521) q[2];
rz(2.1382051) q[3];
sx q[3];
rz(-2.181874) q[3];
sx q[3];
rz(2.3333534) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7830911) q[0];
sx q[0];
rz(-0.49810228) q[0];
sx q[0];
rz(0.95712334) q[0];
rz(-2.762291) q[1];
sx q[1];
rz(-2.7298268) q[1];
sx q[1];
rz(-2.9764825) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6388164) q[0];
sx q[0];
rz(-2.2383225) q[0];
sx q[0];
rz(1.1291885) q[0];
x q[1];
rz(3.0909507) q[2];
sx q[2];
rz(-1.0989185) q[2];
sx q[2];
rz(0.34220055) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.86805008) q[1];
sx q[1];
rz(-2.1434577) q[1];
sx q[1];
rz(-0.21577253) q[1];
rz(-0.96747193) q[3];
sx q[3];
rz(-1.3730197) q[3];
sx q[3];
rz(2.8654049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2740606) q[2];
sx q[2];
rz(-1.2313077) q[2];
sx q[2];
rz(0.77862281) q[2];
rz(-2.8042931) q[3];
sx q[3];
rz(-1.0236579) q[3];
sx q[3];
rz(-1.5571099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.857665) q[0];
sx q[0];
rz(-1.090467) q[0];
sx q[0];
rz(2.0528059) q[0];
rz(-0.42824832) q[1];
sx q[1];
rz(-1.0232404) q[1];
sx q[1];
rz(0.64839378) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8729432) q[0];
sx q[0];
rz(-2.6336484) q[0];
sx q[0];
rz(0.56468876) q[0];
x q[1];
rz(1.6971385) q[2];
sx q[2];
rz(-0.34046945) q[2];
sx q[2];
rz(-0.44035092) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8191145) q[1];
sx q[1];
rz(-1.5286865) q[1];
sx q[1];
rz(-2.4933715) q[1];
rz(-pi) q[2];
rz(-2.5779453) q[3];
sx q[3];
rz(-1.4741338) q[3];
sx q[3];
rz(3.0954297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.26583656) q[2];
sx q[2];
rz(-2.0529604) q[2];
sx q[2];
rz(2.9618373) q[2];
rz(1.1579375) q[3];
sx q[3];
rz(-1.4561184) q[3];
sx q[3];
rz(-1.9013083) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4614048) q[0];
sx q[0];
rz(-2.9869933) q[0];
sx q[0];
rz(-2.3416478) q[0];
rz(1.9752621) q[1];
sx q[1];
rz(-0.98465289) q[1];
sx q[1];
rz(-2.226895) q[1];
rz(-2.2565319) q[2];
sx q[2];
rz(-0.37715465) q[2];
sx q[2];
rz(-2.5436795) q[2];
rz(0.55462454) q[3];
sx q[3];
rz(-2.7477874) q[3];
sx q[3];
rz(1.2870233) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
