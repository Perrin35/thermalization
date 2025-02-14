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
rz(4.1780317) q[0];
sx q[0];
rz(10.960624) q[0];
rz(0.66134557) q[1];
sx q[1];
rz(-1.1868492) q[1];
sx q[1];
rz(0.43651906) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8494328) q[0];
sx q[0];
rz(-1.8884553) q[0];
sx q[0];
rz(-3.0374797) q[0];
rz(-pi) q[1];
rz(1.1182734) q[2];
sx q[2];
rz(-2.088639) q[2];
sx q[2];
rz(2.81942) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.6923209) q[1];
sx q[1];
rz(-0.8682478) q[1];
sx q[1];
rz(-3.0880804) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17702509) q[3];
sx q[3];
rz(-2.2469871) q[3];
sx q[3];
rz(0.9392304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.46301699) q[2];
sx q[2];
rz(-0.43872321) q[2];
sx q[2];
rz(0.92301816) q[2];
rz(-0.065741278) q[3];
sx q[3];
rz(-1.8140503) q[3];
sx q[3];
rz(-1.340723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0098669212) q[0];
sx q[0];
rz(-1.0634402) q[0];
sx q[0];
rz(0.03751066) q[0];
rz(-2.5610979) q[1];
sx q[1];
rz(-0.66941222) q[1];
sx q[1];
rz(2.4494749) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6094309) q[0];
sx q[0];
rz(-1.2464644) q[0];
sx q[0];
rz(-2.956809) q[0];
rz(-pi) q[1];
rz(1.0672928) q[2];
sx q[2];
rz(-1.3238591) q[2];
sx q[2];
rz(-1.7301138) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.0040652635) q[1];
sx q[1];
rz(-2.1430227) q[1];
sx q[1];
rz(2.1768084) q[1];
x q[2];
rz(-1.7440657) q[3];
sx q[3];
rz(-2.3152346) q[3];
sx q[3];
rz(0.92191523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.47824255) q[2];
sx q[2];
rz(-1.9448091) q[2];
sx q[2];
rz(-1.3956068) q[2];
rz(-0.1869525) q[3];
sx q[3];
rz(-1.9713277) q[3];
sx q[3];
rz(0.54005867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4726987) q[0];
sx q[0];
rz(-2.5085594) q[0];
sx q[0];
rz(2.582666) q[0];
rz(0.82866296) q[1];
sx q[1];
rz(-1.5507973) q[1];
sx q[1];
rz(-2.8899946) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5937045) q[0];
sx q[0];
rz(-1.3798837) q[0];
sx q[0];
rz(2.6949791) q[0];
x q[1];
rz(0.34464406) q[2];
sx q[2];
rz(-2.0953669) q[2];
sx q[2];
rz(1.3087496) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.66183749) q[1];
sx q[1];
rz(-2.2244456) q[1];
sx q[1];
rz(0.33262555) q[1];
rz(-0.67686848) q[3];
sx q[3];
rz(-2.1123611) q[3];
sx q[3];
rz(-1.8344017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.805213) q[2];
sx q[2];
rz(-0.42625913) q[2];
sx q[2];
rz(1.4317929) q[2];
rz(-0.21670565) q[3];
sx q[3];
rz(-1.6101937) q[3];
sx q[3];
rz(1.3592892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8828204) q[0];
sx q[0];
rz(-0.50676218) q[0];
sx q[0];
rz(-1.850542) q[0];
rz(-2.3123815) q[1];
sx q[1];
rz(-2.0334838) q[1];
sx q[1];
rz(2.1270027) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0907959) q[0];
sx q[0];
rz(-1.3279339) q[0];
sx q[0];
rz(-0.54625752) q[0];
rz(-pi) q[1];
x q[1];
rz(1.185174) q[2];
sx q[2];
rz(-2.1004538) q[2];
sx q[2];
rz(2.5593064) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.17550771) q[1];
sx q[1];
rz(-2.1455742) q[1];
sx q[1];
rz(3.1360583) q[1];
x q[2];
rz(1.7144373) q[3];
sx q[3];
rz(-2.3206986) q[3];
sx q[3];
rz(-0.088975541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1318876) q[2];
sx q[2];
rz(-1.9405126) q[2];
sx q[2];
rz(-0.76060549) q[2];
rz(2.6724114) q[3];
sx q[3];
rz(-1.6719336) q[3];
sx q[3];
rz(0.73452264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0932662) q[0];
sx q[0];
rz(-0.51437298) q[0];
sx q[0];
rz(0.58116466) q[0];
rz(1.5515074) q[1];
sx q[1];
rz(-0.64684144) q[1];
sx q[1];
rz(0.69492984) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7030554) q[0];
sx q[0];
rz(-1.8037272) q[0];
sx q[0];
rz(-1.2715879) q[0];
x q[1];
rz(1.476711) q[2];
sx q[2];
rz(-1.2796937) q[2];
sx q[2];
rz(-1.7895607) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4798812) q[1];
sx q[1];
rz(-1.5269566) q[1];
sx q[1];
rz(-0.54370329) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7534689) q[3];
sx q[3];
rz(-1.8032866) q[3];
sx q[3];
rz(-0.35216613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.86357) q[2];
sx q[2];
rz(-1.5333971) q[2];
sx q[2];
rz(-0.30280534) q[2];
rz(-1.4263724) q[3];
sx q[3];
rz(-0.36077603) q[3];
sx q[3];
rz(-0.58437955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4645709) q[0];
sx q[0];
rz(-0.40957054) q[0];
sx q[0];
rz(-2.4466178) q[0];
rz(-0.28680828) q[1];
sx q[1];
rz(-2.5359055) q[1];
sx q[1];
rz(-1.274775) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5329842) q[0];
sx q[0];
rz(-1.8814981) q[0];
sx q[0];
rz(2.0771189) q[0];
x q[1];
rz(-2.4430577) q[2];
sx q[2];
rz(-2.4753127) q[2];
sx q[2];
rz(0.65411386) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3466323) q[1];
sx q[1];
rz(-2.2169276) q[1];
sx q[1];
rz(2.6898795) q[1];
rz(-pi) q[2];
rz(1.3319837) q[3];
sx q[3];
rz(-1.0019571) q[3];
sx q[3];
rz(2.5425926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0605165) q[2];
sx q[2];
rz(-0.21848564) q[2];
sx q[2];
rz(2.1742353) q[2];
rz(0.71990144) q[3];
sx q[3];
rz(-0.94341174) q[3];
sx q[3];
rz(-1.2758183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28441456) q[0];
sx q[0];
rz(-3.1107749) q[0];
sx q[0];
rz(-1.6946633) q[0];
rz(-2.6743496) q[1];
sx q[1];
rz(-1.8485565) q[1];
sx q[1];
rz(1.1955059) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23654437) q[0];
sx q[0];
rz(-2.0285188) q[0];
sx q[0];
rz(0.48120705) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31767665) q[2];
sx q[2];
rz(-0.84140771) q[2];
sx q[2];
rz(-2.9175772) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0250281) q[1];
sx q[1];
rz(-2.0691074) q[1];
sx q[1];
rz(-3.0710941) q[1];
rz(1.2310813) q[3];
sx q[3];
rz(-3.1034443) q[3];
sx q[3];
rz(1.7496141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.82470992) q[2];
sx q[2];
rz(-2.3066545) q[2];
sx q[2];
rz(-2.6173124) q[2];
rz(0.25238642) q[3];
sx q[3];
rz(-0.10656825) q[3];
sx q[3];
rz(0.061080385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97279945) q[0];
sx q[0];
rz(-1.1058829) q[0];
sx q[0];
rz(1.2267145) q[0];
rz(-2.3854158) q[1];
sx q[1];
rz(-2.1935479) q[1];
sx q[1];
rz(0.39047584) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9928241) q[0];
sx q[0];
rz(-2.618194) q[0];
sx q[0];
rz(0.54570178) q[0];
rz(-3.0608589) q[2];
sx q[2];
rz(-2.6962523) q[2];
sx q[2];
rz(-2.3824218) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9041216) q[1];
sx q[1];
rz(-1.3573705) q[1];
sx q[1];
rz(-1.2280812) q[1];
x q[2];
rz(-0.94589767) q[3];
sx q[3];
rz(-2.3970553) q[3];
sx q[3];
rz(-2.8685399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.32435027) q[2];
sx q[2];
rz(-0.70049006) q[2];
sx q[2];
rz(-1.7895169) q[2];
rz(2.9663441) q[3];
sx q[3];
rz(-1.8027571) q[3];
sx q[3];
rz(1.2043183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4125724) q[0];
sx q[0];
rz(-2.6798798) q[0];
sx q[0];
rz(0.66226688) q[0];
rz(0.63016713) q[1];
sx q[1];
rz(-1.2799193) q[1];
sx q[1];
rz(-0.23552775) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92162392) q[0];
sx q[0];
rz(-1.7818091) q[0];
sx q[0];
rz(-1.3930265) q[0];
rz(1.6643737) q[2];
sx q[2];
rz(-2.1380599) q[2];
sx q[2];
rz(-2.7964489) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3429361) q[1];
sx q[1];
rz(-1.6321811) q[1];
sx q[1];
rz(3.117901) q[1];
x q[2];
rz(0.8742378) q[3];
sx q[3];
rz(-1.4676508) q[3];
sx q[3];
rz(2.646605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0887289) q[2];
sx q[2];
rz(-0.95053089) q[2];
sx q[2];
rz(-2.5843184) q[2];
rz(-2.3142464) q[3];
sx q[3];
rz(-1.5118303) q[3];
sx q[3];
rz(1.5000337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92852229) q[0];
sx q[0];
rz(-2.3445573) q[0];
sx q[0];
rz(0.56761566) q[0];
rz(2.7396743) q[1];
sx q[1];
rz(-2.4979976) q[1];
sx q[1];
rz(-1.7793122) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4618191) q[0];
sx q[0];
rz(-1.1704233) q[0];
sx q[0];
rz(2.7072565) q[0];
rz(-2.1902607) q[2];
sx q[2];
rz(-1.7835604) q[2];
sx q[2];
rz(-1.484953) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1115845) q[1];
sx q[1];
rz(-1.4923195) q[1];
sx q[1];
rz(-2.2021738) q[1];
rz(-pi) q[2];
rz(2.3990223) q[3];
sx q[3];
rz(-0.97108632) q[3];
sx q[3];
rz(-0.35123435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7207429) q[2];
sx q[2];
rz(-1.2030615) q[2];
sx q[2];
rz(-0.26025772) q[2];
rz(2.4506954) q[3];
sx q[3];
rz(-0.8553718) q[3];
sx q[3];
rz(1.6183841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33547587) q[0];
sx q[0];
rz(-1.5806883) q[0];
sx q[0];
rz(-1.6089532) q[0];
rz(2.7083022) q[1];
sx q[1];
rz(-1.8186124) q[1];
sx q[1];
rz(-1.1846452) q[1];
rz(2.2122907) q[2];
sx q[2];
rz(-0.96972307) q[2];
sx q[2];
rz(-2.0676421) q[2];
rz(-0.79037249) q[3];
sx q[3];
rz(-1.6675259) q[3];
sx q[3];
rz(-0.70601757) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
