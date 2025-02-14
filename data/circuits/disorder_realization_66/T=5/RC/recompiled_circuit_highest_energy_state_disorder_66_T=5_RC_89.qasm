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
rz(3.003886) q[0];
sx q[0];
rz(-2.6258111) q[0];
sx q[0];
rz(3.1006324) q[0];
rz(0.75086683) q[1];
sx q[1];
rz(2.6584396) q[1];
sx q[1];
rz(9.2106342) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1233323) q[0];
sx q[0];
rz(-1.4054448) q[0];
sx q[0];
rz(1.6613217) q[0];
x q[1];
rz(-3.0532367) q[2];
sx q[2];
rz(-1.7680829) q[2];
sx q[2];
rz(-0.046869761) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.029672) q[1];
sx q[1];
rz(-0.73603928) q[1];
sx q[1];
rz(1.8007832) q[1];
rz(2.3365306) q[3];
sx q[3];
rz(-0.65316155) q[3];
sx q[3];
rz(1.6860608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.21138771) q[2];
sx q[2];
rz(-1.7519506) q[2];
sx q[2];
rz(2.1432121) q[2];
rz(-1.8006648) q[3];
sx q[3];
rz(-1.0705592) q[3];
sx q[3];
rz(0.65795952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(2.9407161) q[0];
sx q[0];
rz(-2.2287892) q[0];
sx q[0];
rz(-0.87769133) q[0];
rz(2.8857723) q[1];
sx q[1];
rz(-2.1820549) q[1];
sx q[1];
rz(-2.6928601) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.701544) q[0];
sx q[0];
rz(-2.8075135) q[0];
sx q[0];
rz(0.17091093) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0340292) q[2];
sx q[2];
rz(-1.7834181) q[2];
sx q[2];
rz(2.4580301) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4325368) q[1];
sx q[1];
rz(-2.2655055) q[1];
sx q[1];
rz(1.1251775) q[1];
rz(-0.93495448) q[3];
sx q[3];
rz(-1.9391914) q[3];
sx q[3];
rz(1.4810497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9933219) q[2];
sx q[2];
rz(-1.4485056) q[2];
sx q[2];
rz(-1.5987965) q[2];
rz(-0.79105061) q[3];
sx q[3];
rz(-1.6796716) q[3];
sx q[3];
rz(0.79234523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1087082) q[0];
sx q[0];
rz(-0.59297639) q[0];
sx q[0];
rz(1.2805043) q[0];
rz(2.9325824) q[1];
sx q[1];
rz(-2.0031395) q[1];
sx q[1];
rz(1.7144405) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7775734) q[0];
sx q[0];
rz(-0.74467842) q[0];
sx q[0];
rz(0.68666057) q[0];
rz(3.0328906) q[2];
sx q[2];
rz(-1.4042028) q[2];
sx q[2];
rz(0.76736375) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.64066891) q[1];
sx q[1];
rz(-0.66447778) q[1];
sx q[1];
rz(-0.44698235) q[1];
x q[2];
rz(2.9455037) q[3];
sx q[3];
rz(-0.59299378) q[3];
sx q[3];
rz(-0.65341572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7872494) q[2];
sx q[2];
rz(-2.3086583) q[2];
sx q[2];
rz(-0.99861097) q[2];
rz(0.35501114) q[3];
sx q[3];
rz(-1.7576926) q[3];
sx q[3];
rz(-1.9173054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76096475) q[0];
sx q[0];
rz(-0.17292085) q[0];
sx q[0];
rz(-0.91249102) q[0];
rz(-1.4997222) q[1];
sx q[1];
rz(-1.8955756) q[1];
sx q[1];
rz(1.46371) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74418236) q[0];
sx q[0];
rz(-1.46755) q[0];
sx q[0];
rz(-2.849177) q[0];
rz(-pi) q[1];
rz(0.60419146) q[2];
sx q[2];
rz(-1.6407654) q[2];
sx q[2];
rz(1.7330012) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50747516) q[1];
sx q[1];
rz(-1.9582796) q[1];
sx q[1];
rz(-1.3204367) q[1];
rz(0.2777956) q[3];
sx q[3];
rz(-2.2738289) q[3];
sx q[3];
rz(-2.4206116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0616167) q[2];
sx q[2];
rz(-1.4392263) q[2];
sx q[2];
rz(1.7495135) q[2];
rz(-1.1213087) q[3];
sx q[3];
rz(-2.0544572) q[3];
sx q[3];
rz(0.054718941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012080972) q[0];
sx q[0];
rz(-1.5665781) q[0];
sx q[0];
rz(2.6636301) q[0];
rz(-1.6732008) q[1];
sx q[1];
rz(-1.5719527) q[1];
sx q[1];
rz(-2.0782616) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2713094) q[0];
sx q[0];
rz(-2.0905604) q[0];
sx q[0];
rz(0.55488806) q[0];
rz(-pi) q[1];
rz(2.9929286) q[2];
sx q[2];
rz(-1.0188661) q[2];
sx q[2];
rz(-2.8665598) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8420452) q[1];
sx q[1];
rz(-2.2631524) q[1];
sx q[1];
rz(-2.1238219) q[1];
rz(-pi) q[2];
rz(2.8394978) q[3];
sx q[3];
rz(-1.8595795) q[3];
sx q[3];
rz(2.0022415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.31074506) q[2];
sx q[2];
rz(-1.1771392) q[2];
sx q[2];
rz(-1.2486521) q[2];
rz(-2.2511075) q[3];
sx q[3];
rz(-1.9200446) q[3];
sx q[3];
rz(-0.42731592) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.86522) q[0];
sx q[0];
rz(-0.38580147) q[0];
sx q[0];
rz(1.9355829) q[0];
rz(1.6845866) q[1];
sx q[1];
rz(-2.147069) q[1];
sx q[1];
rz(2.1280033) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.550279) q[0];
sx q[0];
rz(-3.075789) q[0];
sx q[0];
rz(0.71424152) q[0];
x q[1];
rz(-0.59047575) q[2];
sx q[2];
rz(-1.9370859) q[2];
sx q[2];
rz(0.72792887) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9070753) q[1];
sx q[1];
rz(-2.5401118) q[1];
sx q[1];
rz(0.3882577) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7238377) q[3];
sx q[3];
rz(-0.42649999) q[3];
sx q[3];
rz(-0.41446009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4509585) q[2];
sx q[2];
rz(-2.3834507) q[2];
sx q[2];
rz(-3.0699733) q[2];
rz(0.13876638) q[3];
sx q[3];
rz(-1.3774739) q[3];
sx q[3];
rz(2.5630786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-0.11107681) q[0];
sx q[0];
rz(-1.0826305) q[0];
sx q[0];
rz(-1.6336596) q[0];
rz(1.1208447) q[1];
sx q[1];
rz(-1.0034674) q[1];
sx q[1];
rz(2.9337163) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5198398) q[0];
sx q[0];
rz(-0.43448453) q[0];
sx q[0];
rz(1.012201) q[0];
rz(2.6034749) q[2];
sx q[2];
rz(-2.8413515) q[2];
sx q[2];
rz(3.0150692) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.058931303) q[1];
sx q[1];
rz(-1.3724695) q[1];
sx q[1];
rz(1.9160509) q[1];
x q[2];
rz(3.0886544) q[3];
sx q[3];
rz(-1.8013305) q[3];
sx q[3];
rz(-2.9028877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8230744) q[2];
sx q[2];
rz(-1.069331) q[2];
sx q[2];
rz(-2.6413909) q[2];
rz(2.9253166) q[3];
sx q[3];
rz(-1.6612771) q[3];
sx q[3];
rz(-1.1094619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08865393) q[0];
sx q[0];
rz(-1.5733938) q[0];
sx q[0];
rz(-1.1608359) q[0];
rz(-2.898518) q[1];
sx q[1];
rz(-1.4832152) q[1];
sx q[1];
rz(1.9879139) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5192272) q[0];
sx q[0];
rz(-1.3524242) q[0];
sx q[0];
rz(-1.1535991) q[0];
x q[1];
rz(-0.52277182) q[2];
sx q[2];
rz(-0.47975329) q[2];
sx q[2];
rz(2.4186277) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.28330559) q[1];
sx q[1];
rz(-1.518462) q[1];
sx q[1];
rz(-0.01776617) q[1];
rz(-pi) q[2];
rz(-1.6457993) q[3];
sx q[3];
rz(-0.60216367) q[3];
sx q[3];
rz(1.259089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.941075) q[2];
sx q[2];
rz(-0.87824559) q[2];
sx q[2];
rz(-0.78594691) q[2];
rz(-3.1089697) q[3];
sx q[3];
rz(-2.2226108) q[3];
sx q[3];
rz(2.6526764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0068479) q[0];
sx q[0];
rz(-2.5866046) q[0];
sx q[0];
rz(-2.3299589) q[0];
rz(2.6875878) q[1];
sx q[1];
rz(-1.3521399) q[1];
sx q[1];
rz(2.6754191) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3562856) q[0];
sx q[0];
rz(-1.476383) q[0];
sx q[0];
rz(1.0862907) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2978372) q[2];
sx q[2];
rz(-2.5409219) q[2];
sx q[2];
rz(0.39988437) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.3442163) q[1];
sx q[1];
rz(-2.600311) q[1];
sx q[1];
rz(0.75030542) q[1];
x q[2];
rz(-0.85316633) q[3];
sx q[3];
rz(-0.93255842) q[3];
sx q[3];
rz(-0.555942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0780645) q[2];
sx q[2];
rz(-1.177265) q[2];
sx q[2];
rz(-1.9130116) q[2];
rz(-0.56033963) q[3];
sx q[3];
rz(-1.8599583) q[3];
sx q[3];
rz(1.3250215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.233577) q[0];
sx q[0];
rz(-2.9486739) q[0];
sx q[0];
rz(0.31780258) q[0];
rz(0.69215149) q[1];
sx q[1];
rz(-2.0798123) q[1];
sx q[1];
rz(2.1923776) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5458706) q[0];
sx q[0];
rz(-0.4341653) q[0];
sx q[0];
rz(-1.4465141) q[0];
rz(-pi) q[1];
x q[1];
rz(0.91429488) q[2];
sx q[2];
rz(-1.8021447) q[2];
sx q[2];
rz(1.0745688) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.80133841) q[1];
sx q[1];
rz(-2.0726554) q[1];
sx q[1];
rz(-1.5141634) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0824954) q[3];
sx q[3];
rz(-1.4708843) q[3];
sx q[3];
rz(2.8076415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.074241) q[2];
sx q[2];
rz(-0.94506741) q[2];
sx q[2];
rz(-2.2828339) q[2];
rz(0.37352118) q[3];
sx q[3];
rz(-1.3770603) q[3];
sx q[3];
rz(3.0612194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9141948) q[0];
sx q[0];
rz(-2.5136431) q[0];
sx q[0];
rz(0.11970438) q[0];
rz(-1.1788728) q[1];
sx q[1];
rz(-2.1736455) q[1];
sx q[1];
rz(-1.4812352) q[1];
rz(0.79201067) q[2];
sx q[2];
rz(-0.7274193) q[2];
sx q[2];
rz(-0.37972478) q[2];
rz(-1.6681485) q[3];
sx q[3];
rz(-2.23776) q[3];
sx q[3];
rz(0.8791578) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
