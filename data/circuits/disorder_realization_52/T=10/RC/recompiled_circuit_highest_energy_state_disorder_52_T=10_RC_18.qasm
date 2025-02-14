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
rz(-2.4177457) q[0];
sx q[0];
rz(-1.9404193) q[0];
sx q[0];
rz(2.6475651) q[0];
rz(-1.5863034) q[1];
sx q[1];
rz(-2.2145693) q[1];
sx q[1];
rz(2.291099) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9181831) q[0];
sx q[0];
rz(-2.4505058) q[0];
sx q[0];
rz(-1.3084433) q[0];
rz(-pi) q[1];
rz(-2.0319967) q[2];
sx q[2];
rz(-2.1716433) q[2];
sx q[2];
rz(-0.33723649) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1728178) q[1];
sx q[1];
rz(-0.39645312) q[1];
sx q[1];
rz(-0.60703599) q[1];
rz(-pi) q[2];
rz(-1.4981573) q[3];
sx q[3];
rz(-1.9359255) q[3];
sx q[3];
rz(-2.3453494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4832619) q[2];
sx q[2];
rz(-0.61507812) q[2];
sx q[2];
rz(2.1966546) q[2];
rz(3.1350709) q[3];
sx q[3];
rz(-2.3801453) q[3];
sx q[3];
rz(2.884088) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1061123) q[0];
sx q[0];
rz(-0.98863125) q[0];
sx q[0];
rz(-0.60580564) q[0];
rz(0.93217355) q[1];
sx q[1];
rz(-1.4235539) q[1];
sx q[1];
rz(3.0359643) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8307997) q[0];
sx q[0];
rz(-0.89185878) q[0];
sx q[0];
rz(-1.7693158) q[0];
rz(-pi) q[1];
rz(-2.058624) q[2];
sx q[2];
rz(-1.2529071) q[2];
sx q[2];
rz(-2.6018104) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.88813284) q[1];
sx q[1];
rz(-1.2769298) q[1];
sx q[1];
rz(2.6685358) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9808045) q[3];
sx q[3];
rz(-1.2224397) q[3];
sx q[3];
rz(-2.1961371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2848844) q[2];
sx q[2];
rz(-1.363089) q[2];
sx q[2];
rz(-1.523783) q[2];
rz(2.3273322) q[3];
sx q[3];
rz(-1.3596478) q[3];
sx q[3];
rz(-0.8849357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.098323671) q[0];
sx q[0];
rz(-0.79541484) q[0];
sx q[0];
rz(1.5650308) q[0];
rz(-0.99524975) q[1];
sx q[1];
rz(-2.1627656) q[1];
sx q[1];
rz(2.3547122) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3580619) q[0];
sx q[0];
rz(-1.3710877) q[0];
sx q[0];
rz(3.003503) q[0];
x q[1];
rz(-1.6549395) q[2];
sx q[2];
rz(-1.9168789) q[2];
sx q[2];
rz(3.0402407) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7461024) q[1];
sx q[1];
rz(-1.8999892) q[1];
sx q[1];
rz(-1.2389061) q[1];
rz(-pi) q[2];
rz(1.6305109) q[3];
sx q[3];
rz(-1.6303326) q[3];
sx q[3];
rz(-1.2950031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.83871049) q[2];
sx q[2];
rz(-1.4883214) q[2];
sx q[2];
rz(-0.6380471) q[2];
rz(-2.6436515) q[3];
sx q[3];
rz(-1.0718071) q[3];
sx q[3];
rz(2.0518484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2447253) q[0];
sx q[0];
rz(-1.3571955) q[0];
sx q[0];
rz(-3.0778399) q[0];
rz(1.6925192) q[1];
sx q[1];
rz(-1.3056825) q[1];
sx q[1];
rz(-1.3099028) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5777912) q[0];
sx q[0];
rz(-0.088106958) q[0];
sx q[0];
rz(2.7035575) q[0];
rz(-0.74964995) q[2];
sx q[2];
rz(-2.8733265) q[2];
sx q[2];
rz(1.1531354) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3514362) q[1];
sx q[1];
rz(-2.3871941) q[1];
sx q[1];
rz(-0.18138563) q[1];
rz(-2.2579262) q[3];
sx q[3];
rz(-1.2901297) q[3];
sx q[3];
rz(-1.3786045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.070907585) q[2];
sx q[2];
rz(-2.0347774) q[2];
sx q[2];
rz(-0.10565383) q[2];
rz(1.9581155) q[3];
sx q[3];
rz(-2.2453997) q[3];
sx q[3];
rz(-0.67160523) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79528177) q[0];
sx q[0];
rz(-2.2310937) q[0];
sx q[0];
rz(1.7972535) q[0];
rz(-0.67289871) q[1];
sx q[1];
rz(-1.3105323) q[1];
sx q[1];
rz(0.013462822) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56820368) q[0];
sx q[0];
rz(-2.4512198) q[0];
sx q[0];
rz(-1.5296049) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6031262) q[2];
sx q[2];
rz(-0.41439498) q[2];
sx q[2];
rz(0.76782819) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.33121651) q[1];
sx q[1];
rz(-1.5840285) q[1];
sx q[1];
rz(0.41535901) q[1];
x q[2];
rz(-0.3810639) q[3];
sx q[3];
rz(-2.237202) q[3];
sx q[3];
rz(-1.8997418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5992735) q[2];
sx q[2];
rz(-0.66212526) q[2];
sx q[2];
rz(1.8947961) q[2];
rz(0.072619297) q[3];
sx q[3];
rz(-0.19425546) q[3];
sx q[3];
rz(-2.8323925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70960629) q[0];
sx q[0];
rz(-2.015634) q[0];
sx q[0];
rz(-0.11269888) q[0];
rz(0.20507774) q[1];
sx q[1];
rz(-1.3321184) q[1];
sx q[1];
rz(-2.233706) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4604707) q[0];
sx q[0];
rz(-2.2325071) q[0];
sx q[0];
rz(-0.76028334) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5285883) q[2];
sx q[2];
rz(-1.2132436) q[2];
sx q[2];
rz(1.2956217) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.0081387719) q[1];
sx q[1];
rz(-0.63888237) q[1];
sx q[1];
rz(2.9879346) q[1];
x q[2];
rz(0.6648074) q[3];
sx q[3];
rz(-1.7403354) q[3];
sx q[3];
rz(0.72532082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1418566) q[2];
sx q[2];
rz(-2.8047968) q[2];
sx q[2];
rz(3.0736308) q[2];
rz(1.147602) q[3];
sx q[3];
rz(-1.2679029) q[3];
sx q[3];
rz(-1.8806774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1934018) q[0];
sx q[0];
rz(-1.3113439) q[0];
sx q[0];
rz(-0.46352682) q[0];
rz(0.14687982) q[1];
sx q[1];
rz(-1.905922) q[1];
sx q[1];
rz(0.99259496) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79560876) q[0];
sx q[0];
rz(-2.3844686) q[0];
sx q[0];
rz(-0.36619314) q[0];
x q[1];
rz(1.2847177) q[2];
sx q[2];
rz(-1.6353288) q[2];
sx q[2];
rz(-2.4308824) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.22659644) q[1];
sx q[1];
rz(-2.2961756) q[1];
sx q[1];
rz(-2.5195049) q[1];
rz(1.1168043) q[3];
sx q[3];
rz(-2.2510248) q[3];
sx q[3];
rz(-0.46453634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2831882) q[2];
sx q[2];
rz(-1.7121199) q[2];
sx q[2];
rz(2.6962213) q[2];
rz(1.3238268) q[3];
sx q[3];
rz(-2.0711074) q[3];
sx q[3];
rz(-1.0858735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0367947) q[0];
sx q[0];
rz(-0.0093655149) q[0];
sx q[0];
rz(-0.15583663) q[0];
rz(1.2771295) q[1];
sx q[1];
rz(-1.3469478) q[1];
sx q[1];
rz(-3.1256622) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.984337) q[0];
sx q[0];
rz(-0.80828342) q[0];
sx q[0];
rz(1.2356204) q[0];
rz(-pi) q[1];
rz(-2.3492947) q[2];
sx q[2];
rz(-1.3092666) q[2];
sx q[2];
rz(0.025394414) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6813035) q[1];
sx q[1];
rz(-0.88424129) q[1];
sx q[1];
rz(0.89214561) q[1];
x q[2];
rz(2.5941284) q[3];
sx q[3];
rz(-0.71540912) q[3];
sx q[3];
rz(-0.61111952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.62898034) q[2];
sx q[2];
rz(-0.11233687) q[2];
sx q[2];
rz(-3.0158499) q[2];
rz(-0.9564774) q[3];
sx q[3];
rz(-1.4230909) q[3];
sx q[3];
rz(-2.8023348) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9823031) q[0];
sx q[0];
rz(-2.90137) q[0];
sx q[0];
rz(-2.6851728) q[0];
rz(1.5688815) q[1];
sx q[1];
rz(-2.1379037) q[1];
sx q[1];
rz(-0.68663866) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4149041) q[0];
sx q[0];
rz(-2.3285638) q[0];
sx q[0];
rz(-0.70720478) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1988611) q[2];
sx q[2];
rz(-0.21272993) q[2];
sx q[2];
rz(2.5036734) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9514044) q[1];
sx q[1];
rz(-1.4238289) q[1];
sx q[1];
rz(-1.6920056) q[1];
x q[2];
rz(1.6410646) q[3];
sx q[3];
rz(-1.3965551) q[3];
sx q[3];
rz(-1.4605923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.75866428) q[2];
sx q[2];
rz(-2.0186581) q[2];
sx q[2];
rz(-2.0056966) q[2];
rz(-2.1618333) q[3];
sx q[3];
rz(-1.8085248) q[3];
sx q[3];
rz(0.69825828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45282388) q[0];
sx q[0];
rz(-1.8166421) q[0];
sx q[0];
rz(-2.685637) q[0];
rz(0.042757209) q[1];
sx q[1];
rz(-1.9330934) q[1];
sx q[1];
rz(-0.6932238) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1826349) q[0];
sx q[0];
rz(-1.4719677) q[0];
sx q[0];
rz(1.1846428) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1666388) q[2];
sx q[2];
rz(-1.1956788) q[2];
sx q[2];
rz(-0.98304316) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.4148414) q[1];
sx q[1];
rz(-1.5637239) q[1];
sx q[1];
rz(-1.6688136) q[1];
rz(-pi) q[2];
rz(0.45155489) q[3];
sx q[3];
rz(-1.8579036) q[3];
sx q[3];
rz(-2.5352458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.91528714) q[2];
sx q[2];
rz(-1.2682356) q[2];
sx q[2];
rz(0.66429663) q[2];
rz(-1.213446) q[3];
sx q[3];
rz(-2.4197141) q[3];
sx q[3];
rz(-0.87219316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(1.6912457) q[0];
sx q[0];
rz(-2.3007614) q[0];
sx q[0];
rz(1.7578516) q[0];
rz(-1.6830403) q[1];
sx q[1];
rz(-0.28266193) q[1];
sx q[1];
rz(-2.2255486) q[1];
rz(-0.37049313) q[2];
sx q[2];
rz(-1.5975614) q[2];
sx q[2];
rz(2.3536828) q[2];
rz(-2.4356859) q[3];
sx q[3];
rz(-2.0173895) q[3];
sx q[3];
rz(2.8593393) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
