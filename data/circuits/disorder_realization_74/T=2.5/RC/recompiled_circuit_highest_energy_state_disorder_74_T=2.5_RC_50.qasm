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
rz(0.81646252) q[0];
sx q[0];
rz(3.2433885) q[0];
sx q[0];
rz(9.9643702) q[0];
rz(0.49630961) q[1];
sx q[1];
rz(-0.30975431) q[1];
sx q[1];
rz(-2.6024979) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64083034) q[0];
sx q[0];
rz(-1.1295365) q[0];
sx q[0];
rz(1.5426969) q[0];
rz(-pi) q[1];
rz(-0.35959682) q[2];
sx q[2];
rz(-1.7126377) q[2];
sx q[2];
rz(-1.3027018) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2869563) q[1];
sx q[1];
rz(-2.7883734) q[1];
sx q[1];
rz(-1.9853206) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92181262) q[3];
sx q[3];
rz(-1.505449) q[3];
sx q[3];
rz(1.1999038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.76217905) q[2];
sx q[2];
rz(-2.2654686) q[2];
sx q[2];
rz(1.1890821) q[2];
rz(1.1842229) q[3];
sx q[3];
rz(-0.90457478) q[3];
sx q[3];
rz(-1.5466461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9963843) q[0];
sx q[0];
rz(-1.5897911) q[0];
sx q[0];
rz(1.3231963) q[0];
rz(-2.6595751) q[1];
sx q[1];
rz(-2.2299485) q[1];
sx q[1];
rz(0.97420305) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0486748) q[0];
sx q[0];
rz(-0.67103025) q[0];
sx q[0];
rz(1.8604408) q[0];
rz(-pi) q[1];
rz(-0.85136885) q[2];
sx q[2];
rz(-0.15002827) q[2];
sx q[2];
rz(1.139037) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7108591) q[1];
sx q[1];
rz(-1.8054334) q[1];
sx q[1];
rz(-2.5905142) q[1];
rz(-pi) q[2];
rz(1.810435) q[3];
sx q[3];
rz(-1.0175034) q[3];
sx q[3];
rz(-0.026913337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1366068) q[2];
sx q[2];
rz(-2.5544781) q[2];
sx q[2];
rz(-0.74964398) q[2];
rz(-2.7096115) q[3];
sx q[3];
rz(-2.0260729) q[3];
sx q[3];
rz(1.3031134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62464803) q[0];
sx q[0];
rz(-0.2722781) q[0];
sx q[0];
rz(-2.5352449) q[0];
rz(-1.3336522) q[1];
sx q[1];
rz(-1.9275815) q[1];
sx q[1];
rz(2.8820754) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1071651) q[0];
sx q[0];
rz(-1.8176862) q[0];
sx q[0];
rz(0.16040032) q[0];
rz(-pi) q[1];
rz(2.9664842) q[2];
sx q[2];
rz(-1.7705743) q[2];
sx q[2];
rz(1.4005043) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3644581) q[1];
sx q[1];
rz(-1.196047) q[1];
sx q[1];
rz(1.115429) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0084148) q[3];
sx q[3];
rz(-1.9451688) q[3];
sx q[3];
rz(1.238747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2567265) q[2];
sx q[2];
rz(-1.3639516) q[2];
sx q[2];
rz(-1.5474896) q[2];
rz(2.3571842) q[3];
sx q[3];
rz(-1.6268077) q[3];
sx q[3];
rz(-1.5659531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6467658) q[0];
sx q[0];
rz(-2.0443199) q[0];
sx q[0];
rz(0.25949091) q[0];
rz(-1.9208113) q[1];
sx q[1];
rz(-2.6916598) q[1];
sx q[1];
rz(1.4422013) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15029112) q[0];
sx q[0];
rz(-1.1027005) q[0];
sx q[0];
rz(0.97458124) q[0];
rz(-pi) q[1];
rz(1.3791612) q[2];
sx q[2];
rz(-0.60376061) q[2];
sx q[2];
rz(1.1549283) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3842786) q[1];
sx q[1];
rz(-0.43506778) q[1];
sx q[1];
rz(2.5237094) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6569225) q[3];
sx q[3];
rz(-2.4575876) q[3];
sx q[3];
rz(0.67953426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0356902) q[2];
sx q[2];
rz(-1.3373969) q[2];
sx q[2];
rz(-0.40994677) q[2];
rz(1.2116872) q[3];
sx q[3];
rz(-0.68710059) q[3];
sx q[3];
rz(-1.80779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.425151) q[0];
sx q[0];
rz(-0.71617675) q[0];
sx q[0];
rz(-2.8675365) q[0];
rz(-0.43824276) q[1];
sx q[1];
rz(-1.848105) q[1];
sx q[1];
rz(-0.73572198) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0276707) q[0];
sx q[0];
rz(-0.7589853) q[0];
sx q[0];
rz(-0.12667292) q[0];
rz(-pi) q[1];
rz(2.9922036) q[2];
sx q[2];
rz(-1.7138897) q[2];
sx q[2];
rz(2.9516475) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2949038) q[1];
sx q[1];
rz(-1.3219993) q[1];
sx q[1];
rz(1.873068) q[1];
rz(-pi) q[2];
rz(-2.208136) q[3];
sx q[3];
rz(-2.3846755) q[3];
sx q[3];
rz(-0.33801038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2030486) q[2];
sx q[2];
rz(-1.6837348) q[2];
sx q[2];
rz(2.5277444) q[2];
rz(-0.40890536) q[3];
sx q[3];
rz(-2.1730065) q[3];
sx q[3];
rz(-2.3992505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78855377) q[0];
sx q[0];
rz(-2.5034294) q[0];
sx q[0];
rz(1.9045389) q[0];
rz(-1.6963814) q[1];
sx q[1];
rz(-2.684869) q[1];
sx q[1];
rz(2.9663185) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0171368) q[0];
sx q[0];
rz(-1.5457898) q[0];
sx q[0];
rz(-1.3978005) q[0];
x q[1];
rz(2.484283) q[2];
sx q[2];
rz(-2.016531) q[2];
sx q[2];
rz(0.64195871) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3289017) q[1];
sx q[1];
rz(-1.1393424) q[1];
sx q[1];
rz(-2.575909) q[1];
x q[2];
rz(-0.57717429) q[3];
sx q[3];
rz(-1.0374616) q[3];
sx q[3];
rz(-1.8624901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1767629) q[2];
sx q[2];
rz(-1.8076597) q[2];
sx q[2];
rz(0.97243398) q[2];
rz(-1.6644299) q[3];
sx q[3];
rz(-1.1494613) q[3];
sx q[3];
rz(-0.76178637) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85457388) q[0];
sx q[0];
rz(-0.022495689) q[0];
sx q[0];
rz(2.7884685) q[0];
rz(-1.0278206) q[1];
sx q[1];
rz(-1.2007583) q[1];
sx q[1];
rz(-1.0110528) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42732692) q[0];
sx q[0];
rz(-1.2135226) q[0];
sx q[0];
rz(-1.5554122) q[0];
rz(-pi) q[1];
rz(-0.89247668) q[2];
sx q[2];
rz(-0.8304285) q[2];
sx q[2];
rz(2.3811946) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.29652913) q[1];
sx q[1];
rz(-1.4140714) q[1];
sx q[1];
rz(-1.2178376) q[1];
rz(2.3156741) q[3];
sx q[3];
rz(-0.41676918) q[3];
sx q[3];
rz(0.59405223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.31356835) q[2];
sx q[2];
rz(-0.32005969) q[2];
sx q[2];
rz(1.3671406) q[2];
rz(1.8732871) q[3];
sx q[3];
rz(-1.4788078) q[3];
sx q[3];
rz(2.2936599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1139514) q[0];
sx q[0];
rz(-2.5340762) q[0];
sx q[0];
rz(2.0116346) q[0];
rz(-2.9226411) q[1];
sx q[1];
rz(-1.4066701) q[1];
sx q[1];
rz(2.2311282) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2856871) q[0];
sx q[0];
rz(-0.98074161) q[0];
sx q[0];
rz(1.9737712) q[0];
rz(-1.6565645) q[2];
sx q[2];
rz(-1.0940486) q[2];
sx q[2];
rz(-1.3439182) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.76669756) q[1];
sx q[1];
rz(-1.8378705) q[1];
sx q[1];
rz(1.0407991) q[1];
x q[2];
rz(1.9764141) q[3];
sx q[3];
rz(-2.2542103) q[3];
sx q[3];
rz(0.77714099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3160481) q[2];
sx q[2];
rz(-2.3393708) q[2];
sx q[2];
rz(0.82553378) q[2];
rz(2.6801706) q[3];
sx q[3];
rz(-1.8749219) q[3];
sx q[3];
rz(-2.6958444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6398741) q[0];
sx q[0];
rz(-1.8032782) q[0];
sx q[0];
rz(-0.83183944) q[0];
rz(2.4141451) q[1];
sx q[1];
rz(-1.4201545) q[1];
sx q[1];
rz(-0.096253455) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.03503) q[0];
sx q[0];
rz(-2.0756531) q[0];
sx q[0];
rz(-2.6120606) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6321858) q[2];
sx q[2];
rz(-1.824531) q[2];
sx q[2];
rz(-2.369427) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5288594) q[1];
sx q[1];
rz(-1.3599571) q[1];
sx q[1];
rz(-0.092540578) q[1];
rz(0.15877896) q[3];
sx q[3];
rz(-1.3769866) q[3];
sx q[3];
rz(2.4474582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40621296) q[2];
sx q[2];
rz(-1.7071743) q[2];
sx q[2];
rz(1.5926788) q[2];
rz(-2.5088572) q[3];
sx q[3];
rz(-1.2318719) q[3];
sx q[3];
rz(2.8922141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(11/(9*pi)) q[0];
sx q[0];
rz(-1.0405552) q[0];
sx q[0];
rz(-2.7885875) q[0];
rz(0.32265916) q[1];
sx q[1];
rz(-2.8162075) q[1];
sx q[1];
rz(0.37193146) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73814017) q[0];
sx q[0];
rz(-2.7785025) q[0];
sx q[0];
rz(2.1573261) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7499128) q[2];
sx q[2];
rz(-1.8955911) q[2];
sx q[2];
rz(1.3472404) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.10807762) q[1];
sx q[1];
rz(-1.4107499) q[1];
sx q[1];
rz(-2.8041502) q[1];
x q[2];
rz(2.2953565) q[3];
sx q[3];
rz(-1.5414943) q[3];
sx q[3];
rz(-1.3361564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9670664) q[2];
sx q[2];
rz(-1.1897949) q[2];
sx q[2];
rz(-1.8444427) q[2];
rz(0.8477115) q[3];
sx q[3];
rz(-1.1573557) q[3];
sx q[3];
rz(-1.004809) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90947718) q[0];
sx q[0];
rz(-1.3790601) q[0];
sx q[0];
rz(-2.7098304) q[0];
rz(1.3879981) q[1];
sx q[1];
rz(-1.9276062) q[1];
sx q[1];
rz(3.066317) q[1];
rz(2.3222011) q[2];
sx q[2];
rz(-1.0181622) q[2];
sx q[2];
rz(-0.43140585) q[2];
rz(2.7305215) q[3];
sx q[3];
rz(-0.84079327) q[3];
sx q[3];
rz(-0.83415915) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
