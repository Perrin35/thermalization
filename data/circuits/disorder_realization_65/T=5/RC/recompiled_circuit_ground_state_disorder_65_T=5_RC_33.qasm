OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.4545257) q[0];
sx q[0];
rz(5.0205686) q[0];
sx q[0];
rz(10.159259) q[0];
rz(-1.4990516) q[1];
sx q[1];
rz(-1.8416815) q[1];
sx q[1];
rz(-0.97839657) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9578732) q[0];
sx q[0];
rz(-1.5786202) q[0];
sx q[0];
rz(-1.70723) q[0];
rz(-2.1026272) q[2];
sx q[2];
rz(-2.2072993) q[2];
sx q[2];
rz(0.34953618) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.96127993) q[1];
sx q[1];
rz(-0.35113564) q[1];
sx q[1];
rz(1.1457972) q[1];
rz(-pi) q[2];
rz(0.079732883) q[3];
sx q[3];
rz(-2.0679722) q[3];
sx q[3];
rz(2.2261208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3413099) q[2];
sx q[2];
rz(-1.492123) q[2];
sx q[2];
rz(-0.71551234) q[2];
rz(-1.4095151) q[3];
sx q[3];
rz(-0.87604299) q[3];
sx q[3];
rz(1.6375665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3187123) q[0];
sx q[0];
rz(-0.57682288) q[0];
sx q[0];
rz(-0.034828287) q[0];
rz(2.4233129) q[1];
sx q[1];
rz(-1.7363345) q[1];
sx q[1];
rz(1.1911596) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6189335) q[0];
sx q[0];
rz(-0.95574035) q[0];
sx q[0];
rz(0.19490029) q[0];
rz(-pi) q[1];
rz(0.68231742) q[2];
sx q[2];
rz(-1.6860262) q[2];
sx q[2];
rz(0.29733411) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4264448) q[1];
sx q[1];
rz(-2.1472125) q[1];
sx q[1];
rz(1.2395241) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1805405) q[3];
sx q[3];
rz(-2.6724216) q[3];
sx q[3];
rz(2.495595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8027773) q[2];
sx q[2];
rz(-1.85314) q[2];
sx q[2];
rz(-0.046099376) q[2];
rz(0.80373803) q[3];
sx q[3];
rz(-1.9019144) q[3];
sx q[3];
rz(2.5900335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7130724) q[0];
sx q[0];
rz(-1.5621194) q[0];
sx q[0];
rz(0.61306104) q[0];
rz(-0.26519457) q[1];
sx q[1];
rz(-1.3399597) q[1];
sx q[1];
rz(-0.048096098) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6702794) q[0];
sx q[0];
rz(-1.0889772) q[0];
sx q[0];
rz(0.0028208931) q[0];
x q[1];
rz(1.0499642) q[2];
sx q[2];
rz(-2.0893059) q[2];
sx q[2];
rz(1.1074378) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9768757) q[1];
sx q[1];
rz(-1.5111458) q[1];
sx q[1];
rz(-0.05257923) q[1];
rz(2.7907324) q[3];
sx q[3];
rz(-1.8985629) q[3];
sx q[3];
rz(-1.3223437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.6951822) q[2];
sx q[2];
rz(-1.1626652) q[2];
sx q[2];
rz(-2.7282558) q[2];
rz(2.4731596) q[3];
sx q[3];
rz(-1.3276525) q[3];
sx q[3];
rz(1.3774762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17871019) q[0];
sx q[0];
rz(-2.9229735) q[0];
sx q[0];
rz(-0.14337732) q[0];
rz(-2.7168221) q[1];
sx q[1];
rz(-1.2062585) q[1];
sx q[1];
rz(-0.43407789) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.089325) q[0];
sx q[0];
rz(-0.10073034) q[0];
sx q[0];
rz(-1.9585376) q[0];
x q[1];
rz(0.1527359) q[2];
sx q[2];
rz(-0.65780168) q[2];
sx q[2];
rz(-1.7694103) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7460355) q[1];
sx q[1];
rz(-1.0537123) q[1];
sx q[1];
rz(0.16090572) q[1];
rz(-pi) q[2];
rz(0.8114154) q[3];
sx q[3];
rz(-1.4196809) q[3];
sx q[3];
rz(2.9600157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.8726337) q[2];
sx q[2];
rz(-2.3838398) q[2];
sx q[2];
rz(2.7244275) q[2];
rz(2.4560691) q[3];
sx q[3];
rz(-1.439582) q[3];
sx q[3];
rz(-0.050366966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3454469) q[0];
sx q[0];
rz(-2.8247927) q[0];
sx q[0];
rz(1.5078804) q[0];
rz(-0.66868526) q[1];
sx q[1];
rz(-1.1983127) q[1];
sx q[1];
rz(-1.1315469) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2270344) q[0];
sx q[0];
rz(-2.5426604) q[0];
sx q[0];
rz(1.5127439) q[0];
rz(-pi) q[1];
rz(-1.065992) q[2];
sx q[2];
rz(-2.5194296) q[2];
sx q[2];
rz(1.0823859) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.9100489) q[1];
sx q[1];
rz(-2.0374415) q[1];
sx q[1];
rz(1.779516) q[1];
rz(-pi) q[2];
rz(1.442836) q[3];
sx q[3];
rz(-1.6887512) q[3];
sx q[3];
rz(-1.02533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9097627) q[2];
sx q[2];
rz(-0.54002395) q[2];
sx q[2];
rz(-0.46169272) q[2];
rz(-0.26990226) q[3];
sx q[3];
rz(-1.2609127) q[3];
sx q[3];
rz(0.65129483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5676024) q[0];
sx q[0];
rz(-2.7157768) q[0];
sx q[0];
rz(2.0800166) q[0];
rz(-0.65879446) q[1];
sx q[1];
rz(-1.4219204) q[1];
sx q[1];
rz(-0.40491358) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8773738) q[0];
sx q[0];
rz(-1.4658017) q[0];
sx q[0];
rz(-2.4355302) q[0];
x q[1];
rz(-1.986386) q[2];
sx q[2];
rz(-1.0672369) q[2];
sx q[2];
rz(-1.0114618) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9820653) q[1];
sx q[1];
rz(-1.7100466) q[1];
sx q[1];
rz(-1.0604481) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.177321) q[3];
sx q[3];
rz(-2.4654536) q[3];
sx q[3];
rz(-0.87489389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.095526516) q[2];
sx q[2];
rz(-1.0045071) q[2];
sx q[2];
rz(1.3859762) q[2];
rz(-0.69685495) q[3];
sx q[3];
rz(-2.160852) q[3];
sx q[3];
rz(0.81010404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5974469) q[0];
sx q[0];
rz(-2.5283041) q[0];
sx q[0];
rz(2.5469653) q[0];
rz(-1.1320629) q[1];
sx q[1];
rz(-1.7375528) q[1];
sx q[1];
rz(-2.6304257) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60301946) q[0];
sx q[0];
rz(-1.5165189) q[0];
sx q[0];
rz(0.57247573) q[0];
rz(-pi) q[1];
rz(1.2434741) q[2];
sx q[2];
rz(-1.3053857) q[2];
sx q[2];
rz(-2.3922269) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3222924) q[1];
sx q[1];
rz(-0.27931133) q[1];
sx q[1];
rz(-1.929024) q[1];
rz(-pi) q[2];
rz(2.6066512) q[3];
sx q[3];
rz(-0.80748338) q[3];
sx q[3];
rz(0.51715467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1749997) q[2];
sx q[2];
rz(-0.24093691) q[2];
sx q[2];
rz(-0.28755507) q[2];
rz(2.0368841) q[3];
sx q[3];
rz(-1.6738626) q[3];
sx q[3];
rz(3.0159359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48968807) q[0];
sx q[0];
rz(-2.7629485) q[0];
sx q[0];
rz(-2.4878159) q[0];
rz(0.50103029) q[1];
sx q[1];
rz(-2.4153695) q[1];
sx q[1];
rz(2.3562866) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0806431) q[0];
sx q[0];
rz(-1.7082346) q[0];
sx q[0];
rz(2.5713443) q[0];
x q[1];
rz(-0.64898934) q[2];
sx q[2];
rz(-2.4675998) q[2];
sx q[2];
rz(-1.4907995) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.874215) q[1];
sx q[1];
rz(-1.3611882) q[1];
sx q[1];
rz(-2.4506635) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6527965) q[3];
sx q[3];
rz(-1.3424216) q[3];
sx q[3];
rz(-0.60199753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1886957) q[2];
sx q[2];
rz(-1.6152103) q[2];
sx q[2];
rz(2.7719899) q[2];
rz(0.26436198) q[3];
sx q[3];
rz(-1.2814458) q[3];
sx q[3];
rz(-0.00082357426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0109176) q[0];
sx q[0];
rz(-2.4322746) q[0];
sx q[0];
rz(1.5189019) q[0];
rz(-1.9021289) q[1];
sx q[1];
rz(-1.112097) q[1];
sx q[1];
rz(2.1942203) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1993123) q[0];
sx q[0];
rz(-0.40422842) q[0];
sx q[0];
rz(-1.72615) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5079673) q[2];
sx q[2];
rz(-1.9312973) q[2];
sx q[2];
rz(0.45180991) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4650146) q[1];
sx q[1];
rz(-1.4099219) q[1];
sx q[1];
rz(-1.4032928) q[1];
rz(-pi) q[2];
rz(-0.059809879) q[3];
sx q[3];
rz(-1.734077) q[3];
sx q[3];
rz(-2.9684396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.47306791) q[2];
sx q[2];
rz(-2.1418085) q[2];
sx q[2];
rz(1.1119941) q[2];
rz(-2.9396577) q[3];
sx q[3];
rz(-2.5578121) q[3];
sx q[3];
rz(-2.7659265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9428381) q[0];
sx q[0];
rz(-2.9534979) q[0];
sx q[0];
rz(-1.4659708) q[0];
rz(1.6000043) q[1];
sx q[1];
rz(-1.3009289) q[1];
sx q[1];
rz(-2.8864554) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45642325) q[0];
sx q[0];
rz(-0.13482929) q[0];
sx q[0];
rz(1.6998057) q[0];
rz(2.3681247) q[2];
sx q[2];
rz(-1.3313401) q[2];
sx q[2];
rz(-0.6780034) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0855116) q[1];
sx q[1];
rz(-1.2930808) q[1];
sx q[1];
rz(0.081760689) q[1];
rz(-pi) q[2];
rz(0.76948036) q[3];
sx q[3];
rz(-1.7090685) q[3];
sx q[3];
rz(1.6618123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5280891) q[2];
sx q[2];
rz(-1.1897503) q[2];
sx q[2];
rz(-0.36821723) q[2];
rz(-2.8585989) q[3];
sx q[3];
rz(-2.8734983) q[3];
sx q[3];
rz(2.8158269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.014864347) q[0];
sx q[0];
rz(-0.8664425) q[0];
sx q[0];
rz(2.0436825) q[0];
rz(-2.6425843) q[1];
sx q[1];
rz(-1.2393163) q[1];
sx q[1];
rz(-2.7543482) q[1];
rz(-2.784619) q[2];
sx q[2];
rz(-0.35335159) q[2];
sx q[2];
rz(-0.54749827) q[2];
rz(-1.8659429) q[3];
sx q[3];
rz(-2.9316918) q[3];
sx q[3];
rz(1.5422921) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
