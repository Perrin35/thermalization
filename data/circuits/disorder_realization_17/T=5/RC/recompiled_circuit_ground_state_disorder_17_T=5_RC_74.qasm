OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.021304) q[0];
sx q[0];
rz(2.6074183) q[0];
sx q[0];
rz(8.3745126) q[0];
rz(1.8771111) q[1];
sx q[1];
rz(-0.85327947) q[1];
sx q[1];
rz(0.54860151) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7954338) q[0];
sx q[0];
rz(-1.8654658) q[0];
sx q[0];
rz(0.3064038) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5919551) q[2];
sx q[2];
rz(-0.97351626) q[2];
sx q[2];
rz(-0.51366546) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7214515) q[1];
sx q[1];
rz(-2.4597617) q[1];
sx q[1];
rz(-0.89116606) q[1];
rz(0.24031432) q[3];
sx q[3];
rz(-1.263947) q[3];
sx q[3];
rz(0.14123973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7301664) q[2];
sx q[2];
rz(-0.41894087) q[2];
sx q[2];
rz(1.0116928) q[2];
rz(-2.2569979) q[3];
sx q[3];
rz(-1.9832289) q[3];
sx q[3];
rz(1.8959034) q[3];
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
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0737792) q[0];
sx q[0];
rz(-2.1429006) q[0];
sx q[0];
rz(-0.78773898) q[0];
rz(0.9785606) q[1];
sx q[1];
rz(-1.1509044) q[1];
sx q[1];
rz(-1.2501134) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38215559) q[0];
sx q[0];
rz(-0.50001745) q[0];
sx q[0];
rz(1.3445205) q[0];
rz(-1.9735632) q[2];
sx q[2];
rz(-1.8284697) q[2];
sx q[2];
rz(2.7855025) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.82455222) q[1];
sx q[1];
rz(-2.3352156) q[1];
sx q[1];
rz(3.1081852) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9216209) q[3];
sx q[3];
rz(-1.7860054) q[3];
sx q[3];
rz(2.6107091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0222187) q[2];
sx q[2];
rz(-1.4639857) q[2];
sx q[2];
rz(0.12164965) q[2];
rz(0.71074784) q[3];
sx q[3];
rz(-0.23356479) q[3];
sx q[3];
rz(2.5392883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1043333) q[0];
sx q[0];
rz(-2.4132044) q[0];
sx q[0];
rz(-2.4816568) q[0];
rz(-0.51689369) q[1];
sx q[1];
rz(-0.70534244) q[1];
sx q[1];
rz(2.0491811) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3583864) q[0];
sx q[0];
rz(-2.6899245) q[0];
sx q[0];
rz(-0.93675254) q[0];
rz(-pi) q[1];
rz(2.8404929) q[2];
sx q[2];
rz(-2.2458255) q[2];
sx q[2];
rz(-1.1464034) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9398492) q[1];
sx q[1];
rz(-0.99786109) q[1];
sx q[1];
rz(-2.4948289) q[1];
rz(-pi) q[2];
rz(-0.33340065) q[3];
sx q[3];
rz(-0.67615792) q[3];
sx q[3];
rz(-2.859883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2163781) q[2];
sx q[2];
rz(-2.7984012) q[2];
sx q[2];
rz(-1.421831) q[2];
rz(-0.4839932) q[3];
sx q[3];
rz(-1.657594) q[3];
sx q[3];
rz(1.7234195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5928818) q[0];
sx q[0];
rz(-3.0573248) q[0];
sx q[0];
rz(-1.3053869) q[0];
rz(0.0072172324) q[1];
sx q[1];
rz(-2.9199298) q[1];
sx q[1];
rz(-2.4086187) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3696267) q[0];
sx q[0];
rz(-1.7031755) q[0];
sx q[0];
rz(-3.0218829) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4568366) q[2];
sx q[2];
rz(-1.2780398) q[2];
sx q[2];
rz(1.2639015) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0187015) q[1];
sx q[1];
rz(-0.7801434) q[1];
sx q[1];
rz(-2.200526) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6812612) q[3];
sx q[3];
rz(-2.2909431) q[3];
sx q[3];
rz(0.11321774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2923773) q[2];
sx q[2];
rz(-1.7376309) q[2];
sx q[2];
rz(-2.2461069) q[2];
rz(-2.605947) q[3];
sx q[3];
rz(-1.5629385) q[3];
sx q[3];
rz(-0.061804684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2484922) q[0];
sx q[0];
rz(-2.94815) q[0];
sx q[0];
rz(-2.6336811) q[0];
rz(-1.6912564) q[1];
sx q[1];
rz(-1.0117057) q[1];
sx q[1];
rz(-1.3571665) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88801182) q[0];
sx q[0];
rz(-1.5216899) q[0];
sx q[0];
rz(-0.11195575) q[0];
rz(-2.6772887) q[2];
sx q[2];
rz(-1.1912734) q[2];
sx q[2];
rz(-0.10019856) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8352767) q[1];
sx q[1];
rz(-2.7359254) q[1];
sx q[1];
rz(-1.2951281) q[1];
x q[2];
rz(0.37973865) q[3];
sx q[3];
rz(-1.0061227) q[3];
sx q[3];
rz(0.60864757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.42673972) q[2];
sx q[2];
rz(-1.4740976) q[2];
sx q[2];
rz(-1.1023785) q[2];
rz(2.0057996) q[3];
sx q[3];
rz(-1.9250684) q[3];
sx q[3];
rz(-0.77146012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.256846) q[0];
sx q[0];
rz(-2.1717635) q[0];
sx q[0];
rz(0.92887512) q[0];
rz(-0.38201395) q[1];
sx q[1];
rz(-2.1451352) q[1];
sx q[1];
rz(1.4487723) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0576768) q[0];
sx q[0];
rz(-1.2681715) q[0];
sx q[0];
rz(0.50047366) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.63183) q[2];
sx q[2];
rz(-2.054913) q[2];
sx q[2];
rz(-0.33554892) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0121921) q[1];
sx q[1];
rz(-1.4621648) q[1];
sx q[1];
rz(1.9700178) q[1];
x q[2];
rz(1.153451) q[3];
sx q[3];
rz(-2.7355237) q[3];
sx q[3];
rz(1.690133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6021619) q[2];
sx q[2];
rz(-2.7950037) q[2];
sx q[2];
rz(0.76751417) q[2];
rz(-1.8481988) q[3];
sx q[3];
rz(-2.4496205) q[3];
sx q[3];
rz(-2.9366711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-0.9354189) q[0];
sx q[0];
rz(-2.967301) q[0];
sx q[0];
rz(2.8906004) q[0];
rz(-2.9639066) q[1];
sx q[1];
rz(-1.6975941) q[1];
sx q[1];
rz(0.65690717) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2816683) q[0];
sx q[0];
rz(-0.31552464) q[0];
sx q[0];
rz(0.84976999) q[0];
x q[1];
rz(0.20505242) q[2];
sx q[2];
rz(-2.4775213) q[2];
sx q[2];
rz(-0.98552209) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4272139) q[1];
sx q[1];
rz(-1.8698815) q[1];
sx q[1];
rz(-2.7370791) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4961996) q[3];
sx q[3];
rz(-1.2155967) q[3];
sx q[3];
rz(1.3076289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3333007) q[2];
sx q[2];
rz(-0.57922816) q[2];
sx q[2];
rz(-2.2283238) q[2];
rz(0.67780668) q[3];
sx q[3];
rz(-1.6401688) q[3];
sx q[3];
rz(1.0031797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.102757) q[0];
sx q[0];
rz(-1.1433733) q[0];
sx q[0];
rz(0.58297771) q[0];
rz(1.673117) q[1];
sx q[1];
rz(-0.99248326) q[1];
sx q[1];
rz(0.34117064) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7852029) q[0];
sx q[0];
rz(-1.613253) q[0];
sx q[0];
rz(-1.7701985) q[0];
rz(-1.2874574) q[2];
sx q[2];
rz(-1.3379127) q[2];
sx q[2];
rz(-2.8849059) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0778596) q[1];
sx q[1];
rz(-1.9394411) q[1];
sx q[1];
rz(2.7474982) q[1];
rz(-pi) q[2];
rz(-2.2489585) q[3];
sx q[3];
rz(-1.8934819) q[3];
sx q[3];
rz(1.762515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.94656452) q[2];
sx q[2];
rz(-2.3174758) q[2];
sx q[2];
rz(0.80671802) q[2];
rz(0.45754704) q[3];
sx q[3];
rz(-2.1541607) q[3];
sx q[3];
rz(0.88361067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6954527) q[0];
sx q[0];
rz(-1.0574295) q[0];
sx q[0];
rz(2.9651508) q[0];
rz(-2.8314619) q[1];
sx q[1];
rz(-1.6519494) q[1];
sx q[1];
rz(2.2055221) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53545685) q[0];
sx q[0];
rz(-1.274489) q[0];
sx q[0];
rz(1.7264051) q[0];
rz(2.5199982) q[2];
sx q[2];
rz(-0.69658579) q[2];
sx q[2];
rz(0.36795771) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8011978) q[1];
sx q[1];
rz(-2.7615393) q[1];
sx q[1];
rz(1.317522) q[1];
rz(1.3388757) q[3];
sx q[3];
rz(-2.6565564) q[3];
sx q[3];
rz(0.95279653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.04756847) q[2];
sx q[2];
rz(-1.8393643) q[2];
sx q[2];
rz(0.0085208323) q[2];
rz(-0.73355567) q[3];
sx q[3];
rz(-2.3365884) q[3];
sx q[3];
rz(-3.0146397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2215288) q[0];
sx q[0];
rz(-2.7286752) q[0];
sx q[0];
rz(0.76989663) q[0];
rz(0.13433111) q[1];
sx q[1];
rz(-2.3513992) q[1];
sx q[1];
rz(-1.2199527) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7766321) q[0];
sx q[0];
rz(-1.3377995) q[0];
sx q[0];
rz(2.8316336) q[0];
rz(-1.8656857) q[2];
sx q[2];
rz(-2.1572626) q[2];
sx q[2];
rz(2.6457654) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1618859) q[1];
sx q[1];
rz(-2.5570013) q[1];
sx q[1];
rz(2.4120055) q[1];
x q[2];
rz(1.72206) q[3];
sx q[3];
rz(-1.6409887) q[3];
sx q[3];
rz(-2.9541525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8613345) q[2];
sx q[2];
rz(-0.65797776) q[2];
sx q[2];
rz(0.78557837) q[2];
rz(-2.806459) q[3];
sx q[3];
rz(-0.98888713) q[3];
sx q[3];
rz(0.82211632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.817374) q[0];
sx q[0];
rz(-0.86295177) q[0];
sx q[0];
rz(-1.2865768) q[0];
rz(-2.4556976) q[1];
sx q[1];
rz(-0.58882014) q[1];
sx q[1];
rz(-1.3480766) q[1];
rz(-3.1361572) q[2];
sx q[2];
rz(-2.1519244) q[2];
sx q[2];
rz(-1.5627418) q[2];
rz(1.2780581) q[3];
sx q[3];
rz(-1.5354029) q[3];
sx q[3];
rz(-1.7913747) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
