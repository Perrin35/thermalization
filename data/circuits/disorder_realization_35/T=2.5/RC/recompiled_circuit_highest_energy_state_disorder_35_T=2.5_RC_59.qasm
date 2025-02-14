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
rz(2.2284722) q[0];
sx q[0];
rz(-2.1051536) q[0];
sx q[0];
rz(1.5358465) q[0];
rz(-2.4802471) q[1];
sx q[1];
rz(-1.9547434) q[1];
sx q[1];
rz(-0.43651906) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52659863) q[0];
sx q[0];
rz(-0.33373555) q[0];
sx q[0];
rz(-1.8769391) q[0];
x q[1];
rz(-0.5646602) q[2];
sx q[2];
rz(-1.9604948) q[2];
sx q[2];
rz(-1.0124568) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4492717) q[1];
sx q[1];
rz(-0.8682478) q[1];
sx q[1];
rz(-0.0535123) q[1];
rz(2.9645676) q[3];
sx q[3];
rz(-0.89460556) q[3];
sx q[3];
rz(-2.2023622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.46301699) q[2];
sx q[2];
rz(-0.43872321) q[2];
sx q[2];
rz(0.92301816) q[2];
rz(3.0758514) q[3];
sx q[3];
rz(-1.3275423) q[3];
sx q[3];
rz(-1.8008697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.0098669212) q[0];
sx q[0];
rz(-2.0781524) q[0];
sx q[0];
rz(-0.03751066) q[0];
rz(2.5610979) q[1];
sx q[1];
rz(-0.66941222) q[1];
sx q[1];
rz(-2.4494749) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1624533) q[0];
sx q[0];
rz(-1.7458436) q[0];
sx q[0];
rz(1.9003504) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0522444) q[2];
sx q[2];
rz(-2.5855067) q[2];
sx q[2];
rz(-2.5646445) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2381535) q[1];
sx q[1];
rz(-2.3336556) q[1];
sx q[1];
rz(-2.417516) q[1];
x q[2];
rz(-1.397527) q[3];
sx q[3];
rz(-2.3152346) q[3];
sx q[3];
rz(2.2196774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6633501) q[2];
sx q[2];
rz(-1.9448091) q[2];
sx q[2];
rz(-1.3956068) q[2];
rz(0.1869525) q[3];
sx q[3];
rz(-1.9713277) q[3];
sx q[3];
rz(-0.54005867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66889399) q[0];
sx q[0];
rz(-0.63303328) q[0];
sx q[0];
rz(-2.582666) q[0];
rz(0.82866296) q[1];
sx q[1];
rz(-1.5507973) q[1];
sx q[1];
rz(-2.8899946) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54788816) q[0];
sx q[0];
rz(-1.3798837) q[0];
sx q[0];
rz(0.44661354) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0992804) q[2];
sx q[2];
rz(-0.61868514) q[2];
sx q[2];
rz(1.9306192) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4797552) q[1];
sx q[1];
rz(-2.2244456) q[1];
sx q[1];
rz(0.33262555) q[1];
x q[2];
rz(-2.3763856) q[3];
sx q[3];
rz(-0.83929378) q[3];
sx q[3];
rz(-2.8347903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.805213) q[2];
sx q[2];
rz(-0.42625913) q[2];
sx q[2];
rz(-1.4317929) q[2];
rz(-0.21670565) q[3];
sx q[3];
rz(-1.5313989) q[3];
sx q[3];
rz(1.7823035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8828204) q[0];
sx q[0];
rz(-0.50676218) q[0];
sx q[0];
rz(1.2910507) q[0];
rz(-0.82921118) q[1];
sx q[1];
rz(-2.0334838) q[1];
sx q[1];
rz(1.0145899) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050796789) q[0];
sx q[0];
rz(-1.8136588) q[0];
sx q[0];
rz(-0.54625752) q[0];
rz(-pi) q[1];
rz(-1.9564187) q[2];
sx q[2];
rz(-2.1004538) q[2];
sx q[2];
rz(2.5593064) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9660849) q[1];
sx q[1];
rz(-2.1455742) q[1];
sx q[1];
rz(3.1360583) q[1];
rz(-0.75506702) q[3];
sx q[3];
rz(-1.6757378) q[3];
sx q[3];
rz(1.5800832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1318876) q[2];
sx q[2];
rz(-1.9405126) q[2];
sx q[2];
rz(0.76060549) q[2];
rz(-0.46918121) q[3];
sx q[3];
rz(-1.6719336) q[3];
sx q[3];
rz(0.73452264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(2.0932662) q[0];
sx q[0];
rz(-2.6272197) q[0];
sx q[0];
rz(2.560428) q[0];
rz(-1.5515074) q[1];
sx q[1];
rz(-2.4947512) q[1];
sx q[1];
rz(0.69492984) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.061175) q[0];
sx q[0];
rz(-1.8616849) q[0];
sx q[0];
rz(0.24334749) q[0];
rz(-pi) q[1];
rz(-2.8492691) q[2];
sx q[2];
rz(-1.6609123) q[2];
sx q[2];
rz(-2.9499049) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6617115) q[1];
sx q[1];
rz(-1.5269566) q[1];
sx q[1];
rz(2.5978894) q[1];
rz(-pi) q[2];
rz(-1.3881237) q[3];
sx q[3];
rz(-1.338306) q[3];
sx q[3];
rz(-0.35216613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.86357) q[2];
sx q[2];
rz(-1.5333971) q[2];
sx q[2];
rz(-0.30280534) q[2];
rz(-1.4263724) q[3];
sx q[3];
rz(-2.7808166) q[3];
sx q[3];
rz(0.58437955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67702174) q[0];
sx q[0];
rz(-2.7320221) q[0];
sx q[0];
rz(0.69497481) q[0];
rz(2.8547844) q[1];
sx q[1];
rz(-0.60568714) q[1];
sx q[1];
rz(-1.8668176) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54166486) q[0];
sx q[0];
rz(-0.58690161) q[0];
sx q[0];
rz(-0.98595263) q[0];
x q[1];
rz(-2.5996501) q[2];
sx q[2];
rz(-1.9795609) q[2];
sx q[2];
rz(-0.33318502) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.50837738) q[1];
sx q[1];
rz(-1.9267836) q[1];
sx q[1];
rz(-2.2683925) q[1];
x q[2];
rz(-2.5596095) q[3];
sx q[3];
rz(-1.3701539) q[3];
sx q[3];
rz(-2.0394005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0605165) q[2];
sx q[2];
rz(-0.21848564) q[2];
sx q[2];
rz(-2.1742353) q[2];
rz(-2.4216912) q[3];
sx q[3];
rz(-0.94341174) q[3];
sx q[3];
rz(-1.2758183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8571781) q[0];
sx q[0];
rz(-3.1107749) q[0];
sx q[0];
rz(1.6946633) q[0];
rz(-2.6743496) q[1];
sx q[1];
rz(-1.8485565) q[1];
sx q[1];
rz(-1.9460868) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23654437) q[0];
sx q[0];
rz(-2.0285188) q[0];
sx q[0];
rz(2.6603856) q[0];
x q[1];
rz(0.81584064) q[2];
sx q[2];
rz(-1.8058449) q[2];
sx q[2];
rz(1.5624969) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.48796921) q[1];
sx q[1];
rz(-1.5088827) q[1];
sx q[1];
rz(1.0714404) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9105114) q[3];
sx q[3];
rz(-3.1034443) q[3];
sx q[3];
rz(1.3919786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82470992) q[2];
sx q[2];
rz(-2.3066545) q[2];
sx q[2];
rz(-0.52428025) q[2];
rz(-0.25238642) q[3];
sx q[3];
rz(-0.10656825) q[3];
sx q[3];
rz(3.0805123) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97279945) q[0];
sx q[0];
rz(-2.0357098) q[0];
sx q[0];
rz(-1.9148781) q[0];
rz(-0.75617689) q[1];
sx q[1];
rz(-2.1935479) q[1];
sx q[1];
rz(2.7511168) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3813864) q[0];
sx q[0];
rz(-2.0122271) q[0];
sx q[0];
rz(-1.2797829) q[0];
x q[1];
rz(3.0608589) q[2];
sx q[2];
rz(-2.6962523) q[2];
sx q[2];
rz(-0.75917086) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9041216) q[1];
sx q[1];
rz(-1.7842222) q[1];
sx q[1];
rz(-1.2280812) q[1];
x q[2];
rz(2.195695) q[3];
sx q[3];
rz(-0.74453738) q[3];
sx q[3];
rz(2.8685399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8172424) q[2];
sx q[2];
rz(-2.4411026) q[2];
sx q[2];
rz(-1.3520757) q[2];
rz(0.17524854) q[3];
sx q[3];
rz(-1.8027571) q[3];
sx q[3];
rz(1.9372743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4125724) q[0];
sx q[0];
rz(-2.6798798) q[0];
sx q[0];
rz(2.4793258) q[0];
rz(-2.5114255) q[1];
sx q[1];
rz(-1.8616734) q[1];
sx q[1];
rz(0.23552775) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2199687) q[0];
sx q[0];
rz(-1.3597836) q[0];
sx q[0];
rz(-1.7485662) q[0];
rz(-pi) q[1];
rz(-0.56925242) q[2];
sx q[2];
rz(-1.6496837) q[2];
sx q[2];
rz(1.9663262) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3429361) q[1];
sx q[1];
rz(-1.6321811) q[1];
sx q[1];
rz(-3.117901) q[1];
rz(-pi) q[2];
rz(-0.8742378) q[3];
sx q[3];
rz(-1.4676508) q[3];
sx q[3];
rz(0.49498765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0528637) q[2];
sx q[2];
rz(-0.95053089) q[2];
sx q[2];
rz(0.55727422) q[2];
rz(2.3142464) q[3];
sx q[3];
rz(-1.5118303) q[3];
sx q[3];
rz(-1.5000337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92852229) q[0];
sx q[0];
rz(-0.7970354) q[0];
sx q[0];
rz(0.56761566) q[0];
rz(0.40191832) q[1];
sx q[1];
rz(-2.4979976) q[1];
sx q[1];
rz(-1.3622805) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4618191) q[0];
sx q[0];
rz(-1.1704233) q[0];
sx q[0];
rz(-2.7072565) q[0];
x q[1];
rz(2.882233) q[2];
sx q[2];
rz(-0.96732124) q[2];
sx q[2];
rz(-0.063613907) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6580672) q[1];
sx q[1];
rz(-0.9416675) q[1];
sx q[1];
rz(-3.0444798) q[1];
x q[2];
rz(-2.3990223) q[3];
sx q[3];
rz(-0.97108632) q[3];
sx q[3];
rz(-2.7903583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7207429) q[2];
sx q[2];
rz(-1.2030615) q[2];
sx q[2];
rz(2.8813349) q[2];
rz(0.69089729) q[3];
sx q[3];
rz(-0.8553718) q[3];
sx q[3];
rz(-1.6183841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.8061168) q[0];
sx q[0];
rz(-1.5609043) q[0];
sx q[0];
rz(1.5326395) q[0];
rz(2.7083022) q[1];
sx q[1];
rz(-1.8186124) q[1];
sx q[1];
rz(-1.1846452) q[1];
rz(-0.70788371) q[2];
sx q[2];
rz(-1.0546726) q[2];
sx q[2];
rz(2.2451014) q[2];
rz(-2.3512202) q[3];
sx q[3];
rz(-1.4740667) q[3];
sx q[3];
rz(2.4355751) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
