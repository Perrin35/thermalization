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
rz(-0.55627745) q[0];
sx q[0];
rz(-0.039160691) q[0];
sx q[0];
rz(2.477159) q[0];
rz(-0.18222624) q[1];
sx q[1];
rz(-1.4540949) q[1];
sx q[1];
rz(1.7323642) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.389654) q[0];
sx q[0];
rz(-2.1028127) q[0];
sx q[0];
rz(1.7641032) q[0];
x q[1];
rz(-1.9033236) q[2];
sx q[2];
rz(-1.9938139) q[2];
sx q[2];
rz(0.11655434) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4454882) q[1];
sx q[1];
rz(-1.0173807) q[1];
sx q[1];
rz(-1.9457818) q[1];
rz(-0.18376155) q[3];
sx q[3];
rz(-1.9270883) q[3];
sx q[3];
rz(-0.7687062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0613784) q[2];
sx q[2];
rz(-0.55157101) q[2];
sx q[2];
rz(-0.75458327) q[2];
rz(-1.4590229) q[3];
sx q[3];
rz(-0.85586923) q[3];
sx q[3];
rz(-1.7350908) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9262806) q[0];
sx q[0];
rz(-2.5975241) q[0];
sx q[0];
rz(1.9790443) q[0];
rz(2.4303719) q[1];
sx q[1];
rz(-2.4237207) q[1];
sx q[1];
rz(-0.1444764) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4220336) q[0];
sx q[0];
rz(-1.8633909) q[0];
sx q[0];
rz(-1.4841311) q[0];
x q[1];
rz(-2.2644193) q[2];
sx q[2];
rz(-2.2899004) q[2];
sx q[2];
rz(-0.70110496) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.45119845) q[1];
sx q[1];
rz(-1.2022126) q[1];
sx q[1];
rz(0.75779961) q[1];
rz(0.30938332) q[3];
sx q[3];
rz(-2.2094634) q[3];
sx q[3];
rz(2.2974599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7084536) q[2];
sx q[2];
rz(-1.837156) q[2];
sx q[2];
rz(-0.30678314) q[2];
rz(2.3661546) q[3];
sx q[3];
rz(-2.756835) q[3];
sx q[3];
rz(0.14498372) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48258346) q[0];
sx q[0];
rz(-2.6982396) q[0];
sx q[0];
rz(2.3495667) q[0];
rz(-1.5576942) q[1];
sx q[1];
rz(-1.3521103) q[1];
sx q[1];
rz(0.99383324) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4980443) q[0];
sx q[0];
rz(-0.99731612) q[0];
sx q[0];
rz(-1.9924966) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4491399) q[2];
sx q[2];
rz(-0.98441974) q[2];
sx q[2];
rz(2.3687378) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6043338) q[1];
sx q[1];
rz(-0.52642614) q[1];
sx q[1];
rz(-0.98186214) q[1];
rz(-pi) q[2];
x q[2];
rz(2.951521) q[3];
sx q[3];
rz(-1.1255923) q[3];
sx q[3];
rz(-2.0855188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9516248) q[2];
sx q[2];
rz(-0.8364532) q[2];
sx q[2];
rz(2.9659029) q[2];
rz(0.22819337) q[3];
sx q[3];
rz(-0.56932813) q[3];
sx q[3];
rz(0.5051676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93155414) q[0];
sx q[0];
rz(-0.84489548) q[0];
sx q[0];
rz(-1.5616052) q[0];
rz(-2.2700894) q[1];
sx q[1];
rz(-1.611064) q[1];
sx q[1];
rz(-2.9759488) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4145395) q[0];
sx q[0];
rz(-3.139608) q[0];
sx q[0];
rz(2.5059047) q[0];
x q[1];
rz(0.0087426337) q[2];
sx q[2];
rz(-2.0164818) q[2];
sx q[2];
rz(-2.2042556) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1042479) q[1];
sx q[1];
rz(-2.3583721) q[1];
sx q[1];
rz(1.6820289) q[1];
x q[2];
rz(1.8481726) q[3];
sx q[3];
rz(-1.0078537) q[3];
sx q[3];
rz(1.790188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.44276253) q[2];
sx q[2];
rz(-0.58219588) q[2];
sx q[2];
rz(-2.3637135) q[2];
rz(2.2569518) q[3];
sx q[3];
rz(-1.9886465) q[3];
sx q[3];
rz(-0.9790498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.00064656249) q[0];
sx q[0];
rz(-2.1053173) q[0];
sx q[0];
rz(0.86877862) q[0];
rz(-2.2390305) q[1];
sx q[1];
rz(-1.9156009) q[1];
sx q[1];
rz(2.5696519) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40232274) q[0];
sx q[0];
rz(-0.86193854) q[0];
sx q[0];
rz(-2.2174382) q[0];
rz(2.6380013) q[2];
sx q[2];
rz(-1.9707754) q[2];
sx q[2];
rz(1.7970038) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9801367) q[1];
sx q[1];
rz(-1.2677578) q[1];
sx q[1];
rz(-2.0281467) q[1];
x q[2];
rz(2.1359613) q[3];
sx q[3];
rz(-0.4400357) q[3];
sx q[3];
rz(2.7081851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.14704412) q[2];
sx q[2];
rz(-2.2037555) q[2];
sx q[2];
rz(-0.63329548) q[2];
rz(-0.58081943) q[3];
sx q[3];
rz(-1.3588901) q[3];
sx q[3];
rz(0.66494989) q[3];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0112515) q[0];
sx q[0];
rz(-2.7524152) q[0];
sx q[0];
rz(-1.3379541) q[0];
rz(-1.7300216) q[1];
sx q[1];
rz(-2.7346225) q[1];
sx q[1];
rz(0.79708797) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4914843) q[0];
sx q[0];
rz(-0.31768018) q[0];
sx q[0];
rz(-1.0450715) q[0];
rz(1.6930137) q[2];
sx q[2];
rz(-0.56899348) q[2];
sx q[2];
rz(-0.48556604) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.63070272) q[1];
sx q[1];
rz(-1.7308917) q[1];
sx q[1];
rz(-1.2994205) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0510526) q[3];
sx q[3];
rz(-1.7281579) q[3];
sx q[3];
rz(-1.5998154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7756614) q[2];
sx q[2];
rz(-0.66880995) q[2];
sx q[2];
rz(-0.68525165) q[2];
rz(0.25964409) q[3];
sx q[3];
rz(-2.1681867) q[3];
sx q[3];
rz(1.2118305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24945666) q[0];
sx q[0];
rz(-2.9849755) q[0];
sx q[0];
rz(-2.6020965) q[0];
rz(2.9501713) q[1];
sx q[1];
rz(-1.6554662) q[1];
sx q[1];
rz(2.6991381) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6689597) q[0];
sx q[0];
rz(-1.5711391) q[0];
sx q[0];
rz(1.5616722) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0422969) q[2];
sx q[2];
rz(-2.4467154) q[2];
sx q[2];
rz(2.881584) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6574508) q[1];
sx q[1];
rz(-1.2353578) q[1];
sx q[1];
rz(-2.2141333) q[1];
rz(0.78667647) q[3];
sx q[3];
rz(-2.1128251) q[3];
sx q[3];
rz(-1.6539128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.29360867) q[2];
sx q[2];
rz(-2.831735) q[2];
sx q[2];
rz(-2.4388745) q[2];
rz(2.8637049) q[3];
sx q[3];
rz(-2.0746456) q[3];
sx q[3];
rz(1.3387298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.24669312) q[0];
sx q[0];
rz(-2.2791857) q[0];
sx q[0];
rz(2.605751) q[0];
rz(-2.4193173) q[1];
sx q[1];
rz(-1.7012137) q[1];
sx q[1];
rz(2.3586418) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0412126) q[0];
sx q[0];
rz(-1.6579275) q[0];
sx q[0];
rz(1.9308596) q[0];
x q[1];
rz(-1.6269685) q[2];
sx q[2];
rz(-1.1740285) q[2];
sx q[2];
rz(2.6068991) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.49753518) q[1];
sx q[1];
rz(-2.3573207) q[1];
sx q[1];
rz(-0.80401827) q[1];
rz(0.43704982) q[3];
sx q[3];
rz(-1.0659704) q[3];
sx q[3];
rz(-2.7849891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4202412) q[2];
sx q[2];
rz(-2.6768117) q[2];
sx q[2];
rz(-0.50344023) q[2];
rz(-0.33468801) q[3];
sx q[3];
rz(-1.3379593) q[3];
sx q[3];
rz(0.23505841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6675785) q[0];
sx q[0];
rz(-2.773371) q[0];
sx q[0];
rz(-1.0598805) q[0];
rz(-2.8267951) q[1];
sx q[1];
rz(-1.7960725) q[1];
sx q[1];
rz(1.5816636) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2602291) q[0];
sx q[0];
rz(-0.90103982) q[0];
sx q[0];
rz(-1.2854281) q[0];
rz(1.8354206) q[2];
sx q[2];
rz(-1.6034063) q[2];
sx q[2];
rz(2.5196241) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8851354) q[1];
sx q[1];
rz(-0.68084913) q[1];
sx q[1];
rz(1.225994) q[1];
rz(-1.6208036) q[3];
sx q[3];
rz(-2.7259318) q[3];
sx q[3];
rz(-0.60189825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9554837) q[2];
sx q[2];
rz(-1.1649705) q[2];
sx q[2];
rz(-2.1060409) q[2];
rz(-1.1459076) q[3];
sx q[3];
rz(-1.112554) q[3];
sx q[3];
rz(2.9884393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.938852) q[0];
sx q[0];
rz(-0.11883141) q[0];
sx q[0];
rz(-2.3760997) q[0];
rz(-1.0368404) q[1];
sx q[1];
rz(-1.8427269) q[1];
sx q[1];
rz(3.0293005) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85690391) q[0];
sx q[0];
rz(-2.5717126) q[0];
sx q[0];
rz(2.5068552) q[0];
rz(-2.5300171) q[2];
sx q[2];
rz(-2.2627352) q[2];
sx q[2];
rz(-2.126978) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4203305) q[1];
sx q[1];
rz(-0.76618505) q[1];
sx q[1];
rz(1.2797194) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3577317) q[3];
sx q[3];
rz(-1.5376159) q[3];
sx q[3];
rz(0.23058321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3969193) q[2];
sx q[2];
rz(-1.2770709) q[2];
sx q[2];
rz(-3.0199158) q[2];
rz(0.35950288) q[3];
sx q[3];
rz(-0.72327852) q[3];
sx q[3];
rz(-0.017688964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.551238) q[0];
sx q[0];
rz(-1.4689162) q[0];
sx q[0];
rz(1.3015251) q[0];
rz(-1.3941258) q[1];
sx q[1];
rz(-1.9448517) q[1];
sx q[1];
rz(1.3762884) q[1];
rz(-2.5689498) q[2];
sx q[2];
rz(-0.70664684) q[2];
sx q[2];
rz(1.6942506) q[2];
rz(1.2944503) q[3];
sx q[3];
rz(-2.69097) q[3];
sx q[3];
rz(0.16926058) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
