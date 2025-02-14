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
rz(1.16601) q[0];
sx q[0];
rz(-0.87870413) q[0];
sx q[0];
rz(-2.9455844) q[0];
rz(-1.4930383) q[1];
sx q[1];
rz(-0.9382481) q[1];
sx q[1];
rz(-2.2641163) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8931819) q[0];
sx q[0];
rz(-2.0201004) q[0];
sx q[0];
rz(-2.9025159) q[0];
x q[1];
rz(-0.67812596) q[2];
sx q[2];
rz(-2.6161751) q[2];
sx q[2];
rz(3.0209288) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1371378) q[1];
sx q[1];
rz(-0.52912213) q[1];
sx q[1];
rz(2.1251388) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0095413) q[3];
sx q[3];
rz(-0.51902229) q[3];
sx q[3];
rz(-1.2281017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8087372) q[2];
sx q[2];
rz(-0.70781195) q[2];
sx q[2];
rz(-0.037192496) q[2];
rz(2.228179) q[3];
sx q[3];
rz(-2.1593058) q[3];
sx q[3];
rz(1.5090401) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7971802) q[0];
sx q[0];
rz(-2.0672815) q[0];
sx q[0];
rz(0.2221701) q[0];
rz(-0.19368681) q[1];
sx q[1];
rz(-1.2682468) q[1];
sx q[1];
rz(-2.8817315) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3701137) q[0];
sx q[0];
rz(-0.51379097) q[0];
sx q[0];
rz(-2.0822078) q[0];
rz(-pi) q[1];
rz(1.3919984) q[2];
sx q[2];
rz(-1.6769209) q[2];
sx q[2];
rz(0.44580844) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9495595) q[1];
sx q[1];
rz(-1.38654) q[1];
sx q[1];
rz(2.9263587) q[1];
rz(-pi) q[2];
rz(-0.50825676) q[3];
sx q[3];
rz(-2.314869) q[3];
sx q[3];
rz(1.9342157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.76693177) q[2];
sx q[2];
rz(-2.5587475) q[2];
sx q[2];
rz(0.15677162) q[2];
rz(-2.3356656) q[3];
sx q[3];
rz(-1.0904049) q[3];
sx q[3];
rz(-2.6167615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.2369775) q[0];
sx q[0];
rz(-0.1122864) q[0];
sx q[0];
rz(1.1861381) q[0];
rz(3.0778432) q[1];
sx q[1];
rz(-1.6155764) q[1];
sx q[1];
rz(2.729111) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8500687) q[0];
sx q[0];
rz(-2.0578034) q[0];
sx q[0];
rz(-0.97536941) q[0];
rz(2.2507452) q[2];
sx q[2];
rz(-1.4483223) q[2];
sx q[2];
rz(-2.9735931) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.30668614) q[1];
sx q[1];
rz(-1.2962771) q[1];
sx q[1];
rz(-1.7815018) q[1];
rz(-pi) q[2];
rz(0.90610151) q[3];
sx q[3];
rz(-2.0482488) q[3];
sx q[3];
rz(0.64732823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.56846642) q[2];
sx q[2];
rz(-0.44034475) q[2];
sx q[2];
rz(-1.1661412) q[2];
rz(0.79683534) q[3];
sx q[3];
rz(-1.7723869) q[3];
sx q[3];
rz(-2.3603175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.9154938) q[0];
sx q[0];
rz(-1.508536) q[0];
sx q[0];
rz(-1.728212) q[0];
rz(0.49890292) q[1];
sx q[1];
rz(-1.7585124) q[1];
sx q[1];
rz(1.2938719) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8167433) q[0];
sx q[0];
rz(-1.6559282) q[0];
sx q[0];
rz(3.1146634) q[0];
rz(-pi) q[1];
rz(1.8982688) q[2];
sx q[2];
rz(-1.6427543) q[2];
sx q[2];
rz(1.5775934) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1509959) q[1];
sx q[1];
rz(-1.5739723) q[1];
sx q[1];
rz(3.0795829) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3975375) q[3];
sx q[3];
rz(-3.0500055) q[3];
sx q[3];
rz(2.0360144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52900195) q[2];
sx q[2];
rz(-0.68098536) q[2];
sx q[2];
rz(0.16981086) q[2];
rz(-0.42810193) q[3];
sx q[3];
rz(-1.7120818) q[3];
sx q[3];
rz(2.9812109) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47996461) q[0];
sx q[0];
rz(-0.39230883) q[0];
sx q[0];
rz(2.3062134) q[0];
rz(2.5021878) q[1];
sx q[1];
rz(-0.6650005) q[1];
sx q[1];
rz(1.0909874) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8919286) q[0];
sx q[0];
rz(-0.87452379) q[0];
sx q[0];
rz(1.554717) q[0];
rz(-pi) q[1];
rz(-1.7382938) q[2];
sx q[2];
rz(-1.8665627) q[2];
sx q[2];
rz(-2.5633321) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9613232) q[1];
sx q[1];
rz(-1.1774447) q[1];
sx q[1];
rz(0.8913656) q[1];
rz(1.0476607) q[3];
sx q[3];
rz(-1.2854648) q[3];
sx q[3];
rz(-0.77013515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.368448) q[2];
sx q[2];
rz(-1.1870563) q[2];
sx q[2];
rz(0.0040357987) q[2];
rz(-1.076738) q[3];
sx q[3];
rz(-0.69759798) q[3];
sx q[3];
rz(2.9546886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3183597) q[0];
sx q[0];
rz(-1.3613181) q[0];
sx q[0];
rz(-0.63364345) q[0];
rz(0.40114316) q[1];
sx q[1];
rz(-0.70924962) q[1];
sx q[1];
rz(-0.69222442) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.322583) q[0];
sx q[0];
rz(-0.87912175) q[0];
sx q[0];
rz(2.5965967) q[0];
rz(-pi) q[1];
rz(1.879331) q[2];
sx q[2];
rz(-2.3156328) q[2];
sx q[2];
rz(-2.8070297) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.11897381) q[1];
sx q[1];
rz(-1.3217719) q[1];
sx q[1];
rz(-0.98414661) q[1];
x q[2];
rz(-2.9561437) q[3];
sx q[3];
rz(-1.4223961) q[3];
sx q[3];
rz(0.68391358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9358518) q[2];
sx q[2];
rz(-2.3571099) q[2];
sx q[2];
rz(1.9007696) q[2];
rz(-1.1848909) q[3];
sx q[3];
rz(-1.3013867) q[3];
sx q[3];
rz(1.5968116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.009509) q[0];
sx q[0];
rz(-1.8385831) q[0];
sx q[0];
rz(1.3707772) q[0];
rz(2.7235183) q[1];
sx q[1];
rz(-1.4825753) q[1];
sx q[1];
rz(0.48666993) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.226017) q[0];
sx q[0];
rz(-0.56779438) q[0];
sx q[0];
rz(-1.3894807) q[0];
rz(2.1407632) q[2];
sx q[2];
rz(-2.5116133) q[2];
sx q[2];
rz(-2.2905245) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8312953) q[1];
sx q[1];
rz(-2.8716209) q[1];
sx q[1];
rz(0.84862535) q[1];
x q[2];
rz(-3.0492856) q[3];
sx q[3];
rz(-2.669791) q[3];
sx q[3];
rz(-0.86154729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.70249867) q[2];
sx q[2];
rz(-1.0489901) q[2];
sx q[2];
rz(-0.19719633) q[2];
rz(2.6321453) q[3];
sx q[3];
rz(-2.8809437) q[3];
sx q[3];
rz(0.82836866) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8615123) q[0];
sx q[0];
rz(-2.0674288) q[0];
sx q[0];
rz(1.7751088) q[0];
rz(-0.81661433) q[1];
sx q[1];
rz(-2.0520703) q[1];
sx q[1];
rz(2.8588967) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1408843) q[0];
sx q[0];
rz(-1.3907897) q[0];
sx q[0];
rz(-0.96830826) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6655695) q[2];
sx q[2];
rz(-1.2605209) q[2];
sx q[2];
rz(2.5233334) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3779439) q[1];
sx q[1];
rz(-0.70893439) q[1];
sx q[1];
rz(-1.8590742) q[1];
x q[2];
rz(2.5211996) q[3];
sx q[3];
rz(-0.61479688) q[3];
sx q[3];
rz(2.7063587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1064328) q[2];
sx q[2];
rz(-2.1356434) q[2];
sx q[2];
rz(0.98313588) q[2];
rz(-1.2784917) q[3];
sx q[3];
rz(-2.216414) q[3];
sx q[3];
rz(-0.92888752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4814257) q[0];
sx q[0];
rz(-1.7919414) q[0];
sx q[0];
rz(0.31326681) q[0];
rz(0.52472862) q[1];
sx q[1];
rz(-2.8786761) q[1];
sx q[1];
rz(1.252334) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1850644) q[0];
sx q[0];
rz(-1.6913497) q[0];
sx q[0];
rz(-0.53114364) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1519442) q[2];
sx q[2];
rz(-2.2454442) q[2];
sx q[2];
rz(1.0824114) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.58432594) q[1];
sx q[1];
rz(-1.0385822) q[1];
sx q[1];
rz(-0.6289215) q[1];
rz(-pi) q[2];
rz(-3.0602686) q[3];
sx q[3];
rz(-2.3235851) q[3];
sx q[3];
rz(0.52018702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.19920443) q[2];
sx q[2];
rz(-2.51666) q[2];
sx q[2];
rz(-1.3814629) q[2];
rz(2.840461) q[3];
sx q[3];
rz(-0.98265177) q[3];
sx q[3];
rz(-2.7605831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23406601) q[0];
sx q[0];
rz(-1.0028239) q[0];
sx q[0];
rz(1.3475077) q[0];
rz(-2.2566336) q[1];
sx q[1];
rz(-2.5159409) q[1];
sx q[1];
rz(1.7431097) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8162154) q[0];
sx q[0];
rz(-2.7669816) q[0];
sx q[0];
rz(1.5999567) q[0];
rz(2.276097) q[2];
sx q[2];
rz(-1.7112321) q[2];
sx q[2];
rz(0.27117929) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8902379) q[1];
sx q[1];
rz(-1.4904463) q[1];
sx q[1];
rz(1.608196) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1330256) q[3];
sx q[3];
rz(-1.3763714) q[3];
sx q[3];
rz(-1.0510707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0884023) q[2];
sx q[2];
rz(-1.8209499) q[2];
sx q[2];
rz(-2.1957446) q[2];
rz(2.451402) q[3];
sx q[3];
rz(-1.918957) q[3];
sx q[3];
rz(1.3010196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.239554) q[0];
sx q[0];
rz(-1.6033462) q[0];
sx q[0];
rz(-0.70815804) q[0];
rz(-0.061307727) q[1];
sx q[1];
rz(-2.5262482) q[1];
sx q[1];
rz(1.6940438) q[1];
rz(2.8715677) q[2];
sx q[2];
rz(-2.8403828) q[2];
sx q[2];
rz(-2.3595911) q[2];
rz(-1.3907018) q[3];
sx q[3];
rz(-1.4994499) q[3];
sx q[3];
rz(-2.3105619) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
