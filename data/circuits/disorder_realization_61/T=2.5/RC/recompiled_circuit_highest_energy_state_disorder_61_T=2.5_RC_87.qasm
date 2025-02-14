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
rz(0.50245589) q[0];
sx q[0];
rz(-2.034019) q[0];
sx q[0];
rz(-1.7662319) q[0];
rz(-0.31887588) q[1];
sx q[1];
rz(-1.640929) q[1];
sx q[1];
rz(0.73135102) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36424822) q[0];
sx q[0];
rz(-2.4349182) q[0];
sx q[0];
rz(2.0752491) q[0];
rz(3.1362304) q[2];
sx q[2];
rz(-0.66197936) q[2];
sx q[2];
rz(-0.84214166) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3870942) q[1];
sx q[1];
rz(-1.6162786) q[1];
sx q[1];
rz(0.78630739) q[1];
rz(-1.4764686) q[3];
sx q[3];
rz(-1.7324311) q[3];
sx q[3];
rz(1.9548655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2607164) q[2];
sx q[2];
rz(-1.9274351) q[2];
sx q[2];
rz(0.94240776) q[2];
rz(2.2830394) q[3];
sx q[3];
rz(-0.61045727) q[3];
sx q[3];
rz(2.4220991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.061676625) q[0];
sx q[0];
rz(-2.583857) q[0];
sx q[0];
rz(-0.058636531) q[0];
rz(-0.05642852) q[1];
sx q[1];
rz(-1.7516878) q[1];
sx q[1];
rz(-2.0243534) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6031044) q[0];
sx q[0];
rz(-1.177801) q[0];
sx q[0];
rz(-0.80727838) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2710167) q[2];
sx q[2];
rz(-1.6938871) q[2];
sx q[2];
rz(-1.276615) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.27858116) q[1];
sx q[1];
rz(-1.3323235) q[1];
sx q[1];
rz(3.0568074) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9100203) q[3];
sx q[3];
rz(-2.6548214) q[3];
sx q[3];
rz(-2.6269905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.77659082) q[2];
sx q[2];
rz(-0.33010179) q[2];
sx q[2];
rz(-2.7351725) q[2];
rz(-0.31600076) q[3];
sx q[3];
rz(-1.4603115) q[3];
sx q[3];
rz(0.78280226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2773748) q[0];
sx q[0];
rz(-0.45806956) q[0];
sx q[0];
rz(-1.3360485) q[0];
rz(0.27851963) q[1];
sx q[1];
rz(-1.7428935) q[1];
sx q[1];
rz(-2.1727402) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1135546) q[0];
sx q[0];
rz(-2.4062706) q[0];
sx q[0];
rz(-0.57082973) q[0];
rz(-pi) q[1];
rz(-2.3511841) q[2];
sx q[2];
rz(-1.9876754) q[2];
sx q[2];
rz(2.0415155) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.232695) q[1];
sx q[1];
rz(-1.2746342) q[1];
sx q[1];
rz(-0.044065899) q[1];
rz(-2.4163867) q[3];
sx q[3];
rz(-1.7659597) q[3];
sx q[3];
rz(-0.31828398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1442673) q[2];
sx q[2];
rz(-2.6582025) q[2];
sx q[2];
rz(-0.17511314) q[2];
rz(-1.8363606) q[3];
sx q[3];
rz(-0.94279083) q[3];
sx q[3];
rz(-1.2828627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99821943) q[0];
sx q[0];
rz(-2.0046736) q[0];
sx q[0];
rz(-1.0293707) q[0];
rz(0.81946212) q[1];
sx q[1];
rz(-1.9792604) q[1];
sx q[1];
rz(-0.79684657) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5982847) q[0];
sx q[0];
rz(-0.48439041) q[0];
sx q[0];
rz(-2.7192857) q[0];
x q[1];
rz(-1.5810851) q[2];
sx q[2];
rz(-2.6136189) q[2];
sx q[2];
rz(-0.99425232) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0174344) q[1];
sx q[1];
rz(-0.34338152) q[1];
sx q[1];
rz(0.92059532) q[1];
rz(2.4791777) q[3];
sx q[3];
rz(-1.2235421) q[3];
sx q[3];
rz(-1.2353237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0113819) q[2];
sx q[2];
rz(-0.036616651) q[2];
sx q[2];
rz(1.2288564) q[2];
rz(-2.7025488) q[3];
sx q[3];
rz(-1.5758347) q[3];
sx q[3];
rz(-1.2307833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48290408) q[0];
sx q[0];
rz(-1.6084325) q[0];
sx q[0];
rz(-0.80832344) q[0];
rz(-1.3574379) q[1];
sx q[1];
rz(-2.3517377) q[1];
sx q[1];
rz(2.435991) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9997415) q[0];
sx q[0];
rz(-1.6054285) q[0];
sx q[0];
rz(2.8409728) q[0];
x q[1];
rz(-2.2735177) q[2];
sx q[2];
rz(-1.7484312) q[2];
sx q[2];
rz(-2.9736819) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8562634) q[1];
sx q[1];
rz(-1.3176022) q[1];
sx q[1];
rz(0.002753792) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93813809) q[3];
sx q[3];
rz(-2.2626503) q[3];
sx q[3];
rz(0.16975887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6610873) q[2];
sx q[2];
rz(-1.8578119) q[2];
sx q[2];
rz(2.535848) q[2];
rz(3.0211966) q[3];
sx q[3];
rz(-1.2984637) q[3];
sx q[3];
rz(-0.81407434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.8814988) q[0];
sx q[0];
rz(-2.3277178) q[0];
sx q[0];
rz(-0.0023728097) q[0];
rz(-2.782605) q[1];
sx q[1];
rz(-1.9406043) q[1];
sx q[1];
rz(2.733309) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9914382) q[0];
sx q[0];
rz(-1.6065734) q[0];
sx q[0];
rz(-1.6973901) q[0];
rz(-pi) q[1];
rz(-0.26770182) q[2];
sx q[2];
rz(-2.1609339) q[2];
sx q[2];
rz(-0.4877643) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3107016) q[1];
sx q[1];
rz(-2.9017555) q[1];
sx q[1];
rz(0.28530052) q[1];
rz(-pi) q[2];
rz(-0.95127912) q[3];
sx q[3];
rz(-0.86182098) q[3];
sx q[3];
rz(-3.0326642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4957054) q[2];
sx q[2];
rz(-0.79621035) q[2];
sx q[2];
rz(3.0118946) q[2];
rz(2.7221223) q[3];
sx q[3];
rz(-1.2129815) q[3];
sx q[3];
rz(2.3278842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7523338) q[0];
sx q[0];
rz(-0.77521721) q[0];
sx q[0];
rz(2.4580521) q[0];
rz(-0.93841249) q[1];
sx q[1];
rz(-1.3886195) q[1];
sx q[1];
rz(1.2260812) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9447362) q[0];
sx q[0];
rz(-2.4312907) q[0];
sx q[0];
rz(2.5596749) q[0];
rz(-pi) q[1];
x q[1];
rz(0.60003321) q[2];
sx q[2];
rz(-1.4156315) q[2];
sx q[2];
rz(-2.4324291) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6437373) q[1];
sx q[1];
rz(-1.2853936) q[1];
sx q[1];
rz(-1.1021815) q[1];
rz(2.5874373) q[3];
sx q[3];
rz(-0.99923493) q[3];
sx q[3];
rz(-3.1031648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.73416) q[2];
sx q[2];
rz(-1.7075044) q[2];
sx q[2];
rz(-2.4343991) q[2];
rz(1.6678984) q[3];
sx q[3];
rz(-0.67631045) q[3];
sx q[3];
rz(-1.428712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68607512) q[0];
sx q[0];
rz(-2.7954743) q[0];
sx q[0];
rz(-0.92380512) q[0];
rz(-1.0651945) q[1];
sx q[1];
rz(-1.1895475) q[1];
sx q[1];
rz(-1.1346029) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89399983) q[0];
sx q[0];
rz(-1.1957279) q[0];
sx q[0];
rz(-3.0984237) q[0];
rz(2.258111) q[2];
sx q[2];
rz(-1.0687912) q[2];
sx q[2];
rz(-2.4740681) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.75443134) q[1];
sx q[1];
rz(-2.608077) q[1];
sx q[1];
rz(-1.6759765) q[1];
x q[2];
rz(-1.9530961) q[3];
sx q[3];
rz(-2.144882) q[3];
sx q[3];
rz(-0.36575438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.137546) q[2];
sx q[2];
rz(-1.4038439) q[2];
sx q[2];
rz(2.498632) q[2];
rz(-2.9709587) q[3];
sx q[3];
rz(-2.4891487) q[3];
sx q[3];
rz(1.0954866) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74649015) q[0];
sx q[0];
rz(-3.009142) q[0];
sx q[0];
rz(-2.3353031) q[0];
rz(3.131648) q[1];
sx q[1];
rz(-2.5237623) q[1];
sx q[1];
rz(0.26201216) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1073955) q[0];
sx q[0];
rz(-0.55343628) q[0];
sx q[0];
rz(-1.2518) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7315032) q[2];
sx q[2];
rz(-2.2101058) q[2];
sx q[2];
rz(-0.29603816) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.16740522) q[1];
sx q[1];
rz(-1.7767118) q[1];
sx q[1];
rz(2.0108825) q[1];
rz(2.0724107) q[3];
sx q[3];
rz(-1.2805689) q[3];
sx q[3];
rz(-2.888916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.749873) q[2];
sx q[2];
rz(-0.55730692) q[2];
sx q[2];
rz(2.9987175) q[2];
rz(2.2376132) q[3];
sx q[3];
rz(-1.4986135) q[3];
sx q[3];
rz(-0.8814632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7664465) q[0];
sx q[0];
rz(-1.9117993) q[0];
sx q[0];
rz(2.4850856) q[0];
rz(-0.87908602) q[1];
sx q[1];
rz(-1.9706985) q[1];
sx q[1];
rz(-0.62969977) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59081599) q[0];
sx q[0];
rz(-1.970298) q[0];
sx q[0];
rz(2.4062064) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26816396) q[2];
sx q[2];
rz(-1.0536453) q[2];
sx q[2];
rz(-1.6059396) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.15805298) q[1];
sx q[1];
rz(-3.088542) q[1];
sx q[1];
rz(1.0178465) q[1];
rz(2.8318226) q[3];
sx q[3];
rz(-2.0783721) q[3];
sx q[3];
rz(-0.47981706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5188344) q[2];
sx q[2];
rz(-2.3164985) q[2];
sx q[2];
rz(2.4614914) q[2];
rz(0.71637362) q[3];
sx q[3];
rz(-0.10321897) q[3];
sx q[3];
rz(1.1187925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5685365) q[0];
sx q[0];
rz(-1.9087044) q[0];
sx q[0];
rz(-2.9317324) q[0];
rz(-0.14280351) q[1];
sx q[1];
rz(-2.0695984) q[1];
sx q[1];
rz(-2.4048068) q[1];
rz(-0.057351107) q[2];
sx q[2];
rz(-1.3290559) q[2];
sx q[2];
rz(0.17964687) q[2];
rz(3.1013641) q[3];
sx q[3];
rz(-2.1798402) q[3];
sx q[3];
rz(-1.3932651) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
