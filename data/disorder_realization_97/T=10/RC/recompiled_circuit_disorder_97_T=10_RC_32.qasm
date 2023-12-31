OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72397435) q[0];
sx q[0];
rz(-1.6516049) q[0];
sx q[0];
rz(0.93044257) q[0];
rz(0.62970495) q[1];
sx q[1];
rz(4.2760744) q[1];
sx q[1];
rz(8.3174336) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6247834) q[0];
sx q[0];
rz(-2.2342355) q[0];
sx q[0];
rz(-1.1532564) q[0];
rz(-pi) q[1];
x q[1];
rz(0.8503352) q[2];
sx q[2];
rz(-1.2002266) q[2];
sx q[2];
rz(-0.71760273) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.90802898) q[1];
sx q[1];
rz(-1.6828013) q[1];
sx q[1];
rz(1.2909375) q[1];
rz(-pi) q[2];
rz(-0.9355448) q[3];
sx q[3];
rz(-1.2860635) q[3];
sx q[3];
rz(-0.92393827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0779695) q[2];
sx q[2];
rz(-0.72903967) q[2];
sx q[2];
rz(1.8135653) q[2];
rz(-2.8207181) q[3];
sx q[3];
rz(-0.98595536) q[3];
sx q[3];
rz(-0.13197556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(0.48224738) q[0];
sx q[0];
rz(-3.0292065) q[0];
sx q[0];
rz(-2.2609718) q[0];
rz(-1.2940787) q[1];
sx q[1];
rz(-2.7236415) q[1];
sx q[1];
rz(-2.3243288) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5215223) q[0];
sx q[0];
rz(-1.2261454) q[0];
sx q[0];
rz(-2.9917813) q[0];
rz(-pi) q[1];
rz(-1.9643289) q[2];
sx q[2];
rz(-2.6173008) q[2];
sx q[2];
rz(0.92698586) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4546928) q[1];
sx q[1];
rz(-0.94140879) q[1];
sx q[1];
rz(-1.8886186) q[1];
x q[2];
rz(-1.0975295) q[3];
sx q[3];
rz(-1.739199) q[3];
sx q[3];
rz(2.4646204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.91784224) q[2];
sx q[2];
rz(-2.4607401) q[2];
sx q[2];
rz(-2.7775653) q[2];
rz(-0.98637995) q[3];
sx q[3];
rz(-1.4168408) q[3];
sx q[3];
rz(-1.6769489) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7746975) q[0];
sx q[0];
rz(-2.3092473) q[0];
sx q[0];
rz(0.96631518) q[0];
rz(2.9486588) q[1];
sx q[1];
rz(-1.0886334) q[1];
sx q[1];
rz(-1.6945217) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15196249) q[0];
sx q[0];
rz(-1.785779) q[0];
sx q[0];
rz(-1.5528029) q[0];
rz(-pi) q[1];
rz(1.3005199) q[2];
sx q[2];
rz(-0.30105653) q[2];
sx q[2];
rz(1.9384055) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2599064) q[1];
sx q[1];
rz(-0.73598624) q[1];
sx q[1];
rz(2.5658539) q[1];
rz(-1.0560016) q[3];
sx q[3];
rz(-0.47009531) q[3];
sx q[3];
rz(-1.6197636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.43859279) q[2];
sx q[2];
rz(-0.48406988) q[2];
sx q[2];
rz(2.036371) q[2];
rz(0.74622074) q[3];
sx q[3];
rz(-1.4893702) q[3];
sx q[3];
rz(2.1658649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(2.1352017) q[0];
sx q[0];
rz(-2.6171896) q[0];
sx q[0];
rz(1.4659457) q[0];
rz(-0.28494596) q[1];
sx q[1];
rz(-2.0712712) q[1];
sx q[1];
rz(-2.7526061) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31127351) q[0];
sx q[0];
rz(-1.8000326) q[0];
sx q[0];
rz(-0.29005187) q[0];
rz(-pi) q[1];
rz(2.0573425) q[2];
sx q[2];
rz(-1.7440737) q[2];
sx q[2];
rz(-1.7692406) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.092408808) q[1];
sx q[1];
rz(-1.0256983) q[1];
sx q[1];
rz(2.0712907) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3936478) q[3];
sx q[3];
rz(-1.7224632) q[3];
sx q[3];
rz(-0.15350728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4185562) q[2];
sx q[2];
rz(-1.9779466) q[2];
sx q[2];
rz(2.6848865) q[2];
rz(1.6263973) q[3];
sx q[3];
rz(-2.1925192) q[3];
sx q[3];
rz(2.8592498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91963768) q[0];
sx q[0];
rz(-1.1397521) q[0];
sx q[0];
rz(-2.7815681) q[0];
rz(-0.64741627) q[1];
sx q[1];
rz(-1.5292239) q[1];
sx q[1];
rz(-2.6470851) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9279328) q[0];
sx q[0];
rz(-1.4466009) q[0];
sx q[0];
rz(2.859982) q[0];
rz(-3.102166) q[2];
sx q[2];
rz(-1.4354424) q[2];
sx q[2];
rz(-0.53265041) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.33489409) q[1];
sx q[1];
rz(-0.94167275) q[1];
sx q[1];
rz(1.4966399) q[1];
x q[2];
rz(-0.4440998) q[3];
sx q[3];
rz(-1.5937623) q[3];
sx q[3];
rz(2.3820153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.71022025) q[2];
sx q[2];
rz(-0.65588313) q[2];
sx q[2];
rz(0.29850706) q[2];
rz(-2.8295637) q[3];
sx q[3];
rz(-1.3307064) q[3];
sx q[3];
rz(-0.63849866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.9962149) q[0];
sx q[0];
rz(-2.4185116) q[0];
sx q[0];
rz(2.9456855) q[0];
rz(-0.021082489) q[1];
sx q[1];
rz(-1.3985876) q[1];
sx q[1];
rz(1.9063937) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7346749) q[0];
sx q[0];
rz(-1.2376681) q[0];
sx q[0];
rz(1.245265) q[0];
rz(-pi) q[1];
rz(1.9094798) q[2];
sx q[2];
rz(-0.89502305) q[2];
sx q[2];
rz(2.0896926) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0683806) q[1];
sx q[1];
rz(-0.74406032) q[1];
sx q[1];
rz(0.63853227) q[1];
x q[2];
rz(-2.9052832) q[3];
sx q[3];
rz(-1.4893388) q[3];
sx q[3];
rz(2.8636275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.49332508) q[2];
sx q[2];
rz(-2.2154634) q[2];
sx q[2];
rz(-2.7098999) q[2];
rz(1.4124983) q[3];
sx q[3];
rz(-0.72237152) q[3];
sx q[3];
rz(-0.13599642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39111185) q[0];
sx q[0];
rz(-1.9727805) q[0];
sx q[0];
rz(3.0294763) q[0];
rz(0.21513367) q[1];
sx q[1];
rz(-1.5810177) q[1];
sx q[1];
rz(-1.1134061) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68590251) q[0];
sx q[0];
rz(-0.4710353) q[0];
sx q[0];
rz(2.5409565) q[0];
rz(-pi) q[1];
rz(1.8298803) q[2];
sx q[2];
rz(-1.8763181) q[2];
sx q[2];
rz(-2.6079026) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0472764) q[1];
sx q[1];
rz(-2.2502406) q[1];
sx q[1];
rz(-1.1761155) q[1];
x q[2];
rz(0.21344276) q[3];
sx q[3];
rz(-2.0091972) q[3];
sx q[3];
rz(1.9813117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1404184) q[2];
sx q[2];
rz(-0.73264709) q[2];
sx q[2];
rz(0.12602885) q[2];
rz(1.0472939) q[3];
sx q[3];
rz(-1.3207366) q[3];
sx q[3];
rz(2.4333911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4814608) q[0];
sx q[0];
rz(-0.75868693) q[0];
sx q[0];
rz(1.460176) q[0];
rz(-1.2449645) q[1];
sx q[1];
rz(-2.0472725) q[1];
sx q[1];
rz(-1.9326899) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93223244) q[0];
sx q[0];
rz(-2.8673842) q[0];
sx q[0];
rz(-1.2742395) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2751341) q[2];
sx q[2];
rz(-2.2668215) q[2];
sx q[2];
rz(2.8528086) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3371256) q[1];
sx q[1];
rz(-0.85782385) q[1];
sx q[1];
rz(1.1285524) q[1];
rz(2.9221411) q[3];
sx q[3];
rz(-2.1900574) q[3];
sx q[3];
rz(2.0451562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7631491) q[2];
sx q[2];
rz(-1.903879) q[2];
sx q[2];
rz(-1.8219927) q[2];
rz(2.54946) q[3];
sx q[3];
rz(-1.416128) q[3];
sx q[3];
rz(-0.035141703) q[3];
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
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2329344) q[0];
sx q[0];
rz(-2.6180551) q[0];
sx q[0];
rz(-1.3611025) q[0];
rz(-1.9305485) q[1];
sx q[1];
rz(-2.2369604) q[1];
sx q[1];
rz(-2.7499054) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64205326) q[0];
sx q[0];
rz(-1.5649438) q[0];
sx q[0];
rz(-0.43453479) q[0];
x q[1];
rz(0.39536706) q[2];
sx q[2];
rz(-2.6569416) q[2];
sx q[2];
rz(-1.5026827) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.45410941) q[1];
sx q[1];
rz(-0.81067639) q[1];
sx q[1];
rz(2.0337385) q[1];
rz(-pi) q[2];
rz(3.016032) q[3];
sx q[3];
rz(-2.6124622) q[3];
sx q[3];
rz(-0.76847968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6212375) q[2];
sx q[2];
rz(-1.3427799) q[2];
sx q[2];
rz(-1.3367782) q[2];
rz(0.76198602) q[3];
sx q[3];
rz(-2.8218994) q[3];
sx q[3];
rz(2.3412162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(15/(14*pi)) q[0];
sx q[0];
rz(-0.30277345) q[0];
sx q[0];
rz(0.57089943) q[0];
rz(-1.4292498) q[1];
sx q[1];
rz(-2.0776904) q[1];
sx q[1];
rz(-0.16194078) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3085092) q[0];
sx q[0];
rz(-0.97531318) q[0];
sx q[0];
rz(-2.7916629) q[0];
rz(-pi) q[1];
rz(-1.5230721) q[2];
sx q[2];
rz(-0.18880162) q[2];
sx q[2];
rz(1.4023086) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.83878126) q[1];
sx q[1];
rz(-1.4421842) q[1];
sx q[1];
rz(1.7437115) q[1];
rz(-0.52311388) q[3];
sx q[3];
rz(-1.599708) q[3];
sx q[3];
rz(1.6541964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1841715) q[2];
sx q[2];
rz(-2.0364169) q[2];
sx q[2];
rz(-1.7133678) q[2];
rz(1.9421633) q[3];
sx q[3];
rz(-1.0933484) q[3];
sx q[3];
rz(-1.7709581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7006871) q[0];
sx q[0];
rz(-2.6661243) q[0];
sx q[0];
rz(2.0211924) q[0];
rz(-1.7715001) q[1];
sx q[1];
rz(-2.1961828) q[1];
sx q[1];
rz(-0.97074769) q[1];
rz(1.4233521) q[2];
sx q[2];
rz(-1.6885919) q[2];
sx q[2];
rz(-0.13292776) q[2];
rz(-0.45392848) q[3];
sx q[3];
rz(-2.6700927) q[3];
sx q[3];
rz(2.3910458) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
