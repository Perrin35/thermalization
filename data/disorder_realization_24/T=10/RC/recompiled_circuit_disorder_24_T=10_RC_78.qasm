OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.55968094) q[0];
sx q[0];
rz(-1.0868602) q[0];
sx q[0];
rz(-1.342919) q[0];
rz(1.5496594) q[1];
sx q[1];
rz(-0.078443371) q[1];
sx q[1];
rz(2.4989541) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064518236) q[0];
sx q[0];
rz(-2.1097578) q[0];
sx q[0];
rz(-1.5383188) q[0];
rz(-pi) q[1];
rz(0.3354934) q[2];
sx q[2];
rz(-0.53288424) q[2];
sx q[2];
rz(-0.83836183) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5832311) q[1];
sx q[1];
rz(-1.6915295) q[1];
sx q[1];
rz(-1.4503149) q[1];
rz(-pi) q[2];
x q[2];
rz(1.60728) q[3];
sx q[3];
rz(-1.3368703) q[3];
sx q[3];
rz(0.28947383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6671483) q[2];
sx q[2];
rz(-0.64626226) q[2];
sx q[2];
rz(1.8623964) q[2];
rz(0.71875087) q[3];
sx q[3];
rz(-1.5712534) q[3];
sx q[3];
rz(0.32354245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-0.46368018) q[0];
sx q[0];
rz(-1.3543411) q[0];
sx q[0];
rz(2.1610778) q[0];
rz(-0.15788831) q[1];
sx q[1];
rz(-1.4973367) q[1];
sx q[1];
rz(-2.352879) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6533587) q[0];
sx q[0];
rz(-1.4569067) q[0];
sx q[0];
rz(-1.6949523) q[0];
rz(-pi) q[1];
rz(-0.2561432) q[2];
sx q[2];
rz(-1.6015341) q[2];
sx q[2];
rz(2.5275633) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.84856725) q[1];
sx q[1];
rz(-1.2663406) q[1];
sx q[1];
rz(2.5679563) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4103568) q[3];
sx q[3];
rz(-2.2255219) q[3];
sx q[3];
rz(-2.5996641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.366189) q[2];
sx q[2];
rz(-1.0752233) q[2];
sx q[2];
rz(-0.186084) q[2];
rz(-2.4880593) q[3];
sx q[3];
rz(-1.7553522) q[3];
sx q[3];
rz(-3.0236566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2114975) q[0];
sx q[0];
rz(-2.1989172) q[0];
sx q[0];
rz(-2.511456) q[0];
rz(-0.12763003) q[1];
sx q[1];
rz(-2.4928513) q[1];
sx q[1];
rz(-0.72174597) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8080224) q[0];
sx q[0];
rz(-2.434242) q[0];
sx q[0];
rz(-2.7471514) q[0];
x q[1];
rz(2.6038405) q[2];
sx q[2];
rz(-2.3719412) q[2];
sx q[2];
rz(0.69307454) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0112146) q[1];
sx q[1];
rz(-1.5323258) q[1];
sx q[1];
rz(-2.6973666) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33629041) q[3];
sx q[3];
rz(-0.97067562) q[3];
sx q[3];
rz(1.0322514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3815986) q[2];
sx q[2];
rz(-2.1888032) q[2];
sx q[2];
rz(-0.90908137) q[2];
rz(-0.51820731) q[3];
sx q[3];
rz(-2.0776694) q[3];
sx q[3];
rz(-2.853493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58406126) q[0];
sx q[0];
rz(-2.4052305) q[0];
sx q[0];
rz(-0.042536143) q[0];
rz(-2.361239) q[1];
sx q[1];
rz(-2.6413554) q[1];
sx q[1];
rz(-0.76400486) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1799058) q[0];
sx q[0];
rz(-0.67877239) q[0];
sx q[0];
rz(-1.5508482) q[0];
rz(-pi) q[1];
rz(1.4444703) q[2];
sx q[2];
rz(-0.9300803) q[2];
sx q[2];
rz(1.7510406) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2149787) q[1];
sx q[1];
rz(-1.1264631) q[1];
sx q[1];
rz(-0.39919969) q[1];
rz(-pi) q[2];
rz(0.56960168) q[3];
sx q[3];
rz(-2.3929425) q[3];
sx q[3];
rz(-0.99018712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2970695) q[2];
sx q[2];
rz(-1.917118) q[2];
sx q[2];
rz(-0.53304535) q[2];
rz(2.879203) q[3];
sx q[3];
rz(-0.636594) q[3];
sx q[3];
rz(-2.6791402) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2247291) q[0];
sx q[0];
rz(-2.841195) q[0];
sx q[0];
rz(-1.8977144) q[0];
rz(0.22661701) q[1];
sx q[1];
rz(-0.77426568) q[1];
sx q[1];
rz(-0.40333834) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82542244) q[0];
sx q[0];
rz(-0.24992019) q[0];
sx q[0];
rz(0.73096801) q[0];
rz(-pi) q[1];
rz(1.2071768) q[2];
sx q[2];
rz(-2.223613) q[2];
sx q[2];
rz(-1.2448685) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9215645) q[1];
sx q[1];
rz(-1.0462531) q[1];
sx q[1];
rz(1.3925874) q[1];
rz(-pi) q[2];
rz(-1.0395398) q[3];
sx q[3];
rz(-0.96635339) q[3];
sx q[3];
rz(1.2951375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.32101813) q[2];
sx q[2];
rz(-2.0901188) q[2];
sx q[2];
rz(-2.3590951) q[2];
rz(1.1123505) q[3];
sx q[3];
rz(-2.7045043) q[3];
sx q[3];
rz(-2.1330244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.430442) q[0];
sx q[0];
rz(-2.7395881) q[0];
sx q[0];
rz(-0.45561403) q[0];
rz(3.0420711) q[1];
sx q[1];
rz(-1.1385304) q[1];
sx q[1];
rz(-0.0064370357) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59505263) q[0];
sx q[0];
rz(-1.8227302) q[0];
sx q[0];
rz(0.76058723) q[0];
rz(-pi) q[1];
rz(-1.3077277) q[2];
sx q[2];
rz(-2.7527713) q[2];
sx q[2];
rz(0.013465492) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.16033123) q[1];
sx q[1];
rz(-0.48109522) q[1];
sx q[1];
rz(-1.0465924) q[1];
rz(2.3859343) q[3];
sx q[3];
rz(-1.0010127) q[3];
sx q[3];
rz(-0.33982402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7727938) q[2];
sx q[2];
rz(-1.6364748) q[2];
sx q[2];
rz(2.2951365) q[2];
rz(-2.1438697) q[3];
sx q[3];
rz(-2.3322451) q[3];
sx q[3];
rz(-1.6830106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.098175123) q[0];
sx q[0];
rz(-2.257453) q[0];
sx q[0];
rz(0.055710677) q[0];
rz(-0.78272351) q[1];
sx q[1];
rz(-1.1306154) q[1];
sx q[1];
rz(-1.1605211) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6877277) q[0];
sx q[0];
rz(-1.9181607) q[0];
sx q[0];
rz(-2.0485282) q[0];
rz(2.0918526) q[2];
sx q[2];
rz(-1.9869291) q[2];
sx q[2];
rz(-3.1396438) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.22735587) q[1];
sx q[1];
rz(-1.0096692) q[1];
sx q[1];
rz(0.44740541) q[1];
rz(0.79302391) q[3];
sx q[3];
rz(-0.51699713) q[3];
sx q[3];
rz(-1.5152064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0760076) q[2];
sx q[2];
rz(-2.2202754) q[2];
sx q[2];
rz(-2.356142) q[2];
rz(-0.75585946) q[3];
sx q[3];
rz(-1.8463585) q[3];
sx q[3];
rz(2.7261962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8287559) q[0];
sx q[0];
rz(-0.090933196) q[0];
sx q[0];
rz(-1.998741) q[0];
rz(-1.3061334) q[1];
sx q[1];
rz(-1.9530692) q[1];
sx q[1];
rz(-2.725504) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7284262) q[0];
sx q[0];
rz(-0.51983716) q[0];
sx q[0];
rz(-1.1632989) q[0];
x q[1];
rz(2.0390688) q[2];
sx q[2];
rz(-1.6981914) q[2];
sx q[2];
rz(2.1786736) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9302952) q[1];
sx q[1];
rz(-0.43356178) q[1];
sx q[1];
rz(-0.3890721) q[1];
x q[2];
rz(3.0307426) q[3];
sx q[3];
rz(-1.1006315) q[3];
sx q[3];
rz(1.449031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0351506) q[2];
sx q[2];
rz(-1.687655) q[2];
sx q[2];
rz(2.8184334) q[2];
rz(-2.9390826) q[3];
sx q[3];
rz(-1.2812307) q[3];
sx q[3];
rz(2.5642853) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37339661) q[0];
sx q[0];
rz(-2.51238) q[0];
sx q[0];
rz(-1.1449822) q[0];
rz(1.1960944) q[1];
sx q[1];
rz(-2.9856666) q[1];
sx q[1];
rz(2.6224565) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25871823) q[0];
sx q[0];
rz(-1.9472194) q[0];
sx q[0];
rz(-0.14278485) q[0];
x q[1];
rz(-0.22186188) q[2];
sx q[2];
rz(-1.3405181) q[2];
sx q[2];
rz(0.78267539) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5362894) q[1];
sx q[1];
rz(-1.481206) q[1];
sx q[1];
rz(0.78219608) q[1];
rz(-1.7464697) q[3];
sx q[3];
rz(-1.8050977) q[3];
sx q[3];
rz(-2.7845886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8273932) q[2];
sx q[2];
rz(-1.9284356) q[2];
sx q[2];
rz(-0.040977565) q[2];
rz(-0.86769062) q[3];
sx q[3];
rz(-2.4882081) q[3];
sx q[3];
rz(-2.6303671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.7460019) q[0];
sx q[0];
rz(-0.87930167) q[0];
sx q[0];
rz(-1.6145153) q[0];
rz(1.7136259) q[1];
sx q[1];
rz(-2.7797996) q[1];
sx q[1];
rz(-1.6428927) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41560995) q[0];
sx q[0];
rz(-1.590953) q[0];
sx q[0];
rz(1.495342) q[0];
rz(1.7259049) q[2];
sx q[2];
rz(-1.6181706) q[2];
sx q[2];
rz(2.2021289) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0000671) q[1];
sx q[1];
rz(-2.2728517) q[1];
sx q[1];
rz(-2.6932004) q[1];
x q[2];
rz(0.54159553) q[3];
sx q[3];
rz(-2.9062727) q[3];
sx q[3];
rz(0.02863392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.75858086) q[2];
sx q[2];
rz(-1.0772394) q[2];
sx q[2];
rz(-2.995058) q[2];
rz(0.81418973) q[3];
sx q[3];
rz(-0.79629961) q[3];
sx q[3];
rz(0.41671419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.9183337) q[0];
sx q[0];
rz(-1.8948566) q[0];
sx q[0];
rz(-2.1444453) q[0];
rz(1.2659484) q[1];
sx q[1];
rz(-2.1389778) q[1];
sx q[1];
rz(-1.9139342) q[1];
rz(-1.1351893) q[2];
sx q[2];
rz(-2.8786567) q[2];
sx q[2];
rz(-1.565956) q[2];
rz(-2.4155865) q[3];
sx q[3];
rz(-2.2719759) q[3];
sx q[3];
rz(-2.0235973) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];