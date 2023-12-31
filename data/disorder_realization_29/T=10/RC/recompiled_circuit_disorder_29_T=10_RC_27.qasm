OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1155137) q[0];
sx q[0];
rz(-1.4839412) q[0];
sx q[0];
rz(-0.32615647) q[0];
rz(1.9510608) q[1];
sx q[1];
rz(-1.7915373) q[1];
sx q[1];
rz(-1.5426558) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21539772) q[0];
sx q[0];
rz(-0.30456671) q[0];
sx q[0];
rz(0.79415168) q[0];
rz(-pi) q[1];
rz(-1.8250263) q[2];
sx q[2];
rz(-0.98346522) q[2];
sx q[2];
rz(-1.4884782) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1118288) q[1];
sx q[1];
rz(-1.3297237) q[1];
sx q[1];
rz(1.8087216) q[1];
rz(-0.013734038) q[3];
sx q[3];
rz(-2.3547958) q[3];
sx q[3];
rz(2.9247583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.8618384) q[2];
sx q[2];
rz(-2.162231) q[2];
sx q[2];
rz(-0.88511434) q[2];
rz(2.4195813) q[3];
sx q[3];
rz(-1.6885898) q[3];
sx q[3];
rz(3.1341781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2136114) q[0];
sx q[0];
rz(-2.1827224) q[0];
sx q[0];
rz(-2.0425178) q[0];
rz(-2.4765769) q[1];
sx q[1];
rz(-1.7275093) q[1];
sx q[1];
rz(0.87759334) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35560247) q[0];
sx q[0];
rz(-0.88760932) q[0];
sx q[0];
rz(1.8049294) q[0];
rz(2.8009731) q[2];
sx q[2];
rz(-0.38481958) q[2];
sx q[2];
rz(-2.91586) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7922389) q[1];
sx q[1];
rz(-1.5374447) q[1];
sx q[1];
rz(-1.6778212) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1869933) q[3];
sx q[3];
rz(-1.882453) q[3];
sx q[3];
rz(-0.66031885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7426804) q[2];
sx q[2];
rz(-1.8360527) q[2];
sx q[2];
rz(2.0111283) q[2];
rz(1.8418664) q[3];
sx q[3];
rz(-1.9176509) q[3];
sx q[3];
rz(1.6931504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-2.995342) q[0];
sx q[0];
rz(-1.8595707) q[0];
sx q[0];
rz(-0.28999844) q[0];
rz(0.6668123) q[1];
sx q[1];
rz(-2.1077483) q[1];
sx q[1];
rz(-0.07382948) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2436115) q[0];
sx q[0];
rz(-1.5092761) q[0];
sx q[0];
rz(-3.0520526) q[0];
rz(-pi) q[1];
rz(-0.87186558) q[2];
sx q[2];
rz(-2.7060894) q[2];
sx q[2];
rz(-2.0603902) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2245582) q[1];
sx q[1];
rz(-2.9552166) q[1];
sx q[1];
rz(-1.264155) q[1];
rz(-pi) q[2];
rz(0.35145268) q[3];
sx q[3];
rz(-1.6822527) q[3];
sx q[3];
rz(-0.66672882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4905711) q[2];
sx q[2];
rz(-1.9475503) q[2];
sx q[2];
rz(-2.4948965) q[2];
rz(2.0329287) q[3];
sx q[3];
rz(-0.78287786) q[3];
sx q[3];
rz(2.0126608) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17764238) q[0];
sx q[0];
rz(-2.9681866) q[0];
sx q[0];
rz(1.9529163) q[0];
rz(2.1229318) q[1];
sx q[1];
rz(-2.1689292) q[1];
sx q[1];
rz(1.4368988) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1543717) q[0];
sx q[0];
rz(-0.73723307) q[0];
sx q[0];
rz(0.1371951) q[0];
rz(-pi) q[1];
rz(1.2595348) q[2];
sx q[2];
rz(-1.8831823) q[2];
sx q[2];
rz(-2.600008) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1632299) q[1];
sx q[1];
rz(-2.0320738) q[1];
sx q[1];
rz(-2.7510838) q[1];
x q[2];
rz(3.0152263) q[3];
sx q[3];
rz(-1.7003945) q[3];
sx q[3];
rz(1.453293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.58549762) q[2];
sx q[2];
rz(-1.7079587) q[2];
sx q[2];
rz(-2.0193224) q[2];
rz(2.1155817) q[3];
sx q[3];
rz(-0.75338537) q[3];
sx q[3];
rz(2.1508353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4595903) q[0];
sx q[0];
rz(-2.2408709) q[0];
sx q[0];
rz(-2.7686152) q[0];
rz(-2.9176118) q[1];
sx q[1];
rz(-1.1898899) q[1];
sx q[1];
rz(-1.8251098) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3217357) q[0];
sx q[0];
rz(-0.61367354) q[0];
sx q[0];
rz(0.86921285) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26222783) q[2];
sx q[2];
rz(-1.3272459) q[2];
sx q[2];
rz(-1.3893931) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.74354467) q[1];
sx q[1];
rz(-3.0882356) q[1];
sx q[1];
rz(1.3392901) q[1];
rz(-pi) q[2];
rz(2.3718194) q[3];
sx q[3];
rz(-0.7036182) q[3];
sx q[3];
rz(0.15042703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.10107723) q[2];
sx q[2];
rz(-0.83156362) q[2];
sx q[2];
rz(-2.3357847) q[2];
rz(-0.53330437) q[3];
sx q[3];
rz(-2.0093982) q[3];
sx q[3];
rz(2.1300952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39078113) q[0];
sx q[0];
rz(-1.823714) q[0];
sx q[0];
rz(0.090963013) q[0];
rz(2.2816351) q[1];
sx q[1];
rz(-1.1227612) q[1];
sx q[1];
rz(1.3202753) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0549714) q[0];
sx q[0];
rz(-1.939659) q[0];
sx q[0];
rz(-0.47199179) q[0];
rz(-pi) q[1];
rz(0.27310246) q[2];
sx q[2];
rz(-2.5540076) q[2];
sx q[2];
rz(-2.7762129) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.49530416) q[1];
sx q[1];
rz(-2.1332281) q[1];
sx q[1];
rz(2.5763987) q[1];
x q[2];
rz(2.8517013) q[3];
sx q[3];
rz(-0.93512669) q[3];
sx q[3];
rz(1.0155201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0423353) q[2];
sx q[2];
rz(-0.96747413) q[2];
sx q[2];
rz(0.60097224) q[2];
rz(0.48505923) q[3];
sx q[3];
rz(-2.9197013) q[3];
sx q[3];
rz(1.6962359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27959529) q[0];
sx q[0];
rz(-1.1789362) q[0];
sx q[0];
rz(2.5860508) q[0];
rz(-0.034596054) q[1];
sx q[1];
rz(-0.75841537) q[1];
sx q[1];
rz(1.3909891) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44156528) q[0];
sx q[0];
rz(-2.1040943) q[0];
sx q[0];
rz(-1.9576859) q[0];
rz(-1.820307) q[2];
sx q[2];
rz(-1.2148641) q[2];
sx q[2];
rz(-2.3552259) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.52303752) q[1];
sx q[1];
rz(-2.4589775) q[1];
sx q[1];
rz(-1.7726266) q[1];
rz(-2.3659412) q[3];
sx q[3];
rz(-1.6240048) q[3];
sx q[3];
rz(-0.42078161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.81327072) q[2];
sx q[2];
rz(-2.6997824) q[2];
sx q[2];
rz(0.39548809) q[2];
rz(1.288712) q[3];
sx q[3];
rz(-1.6059395) q[3];
sx q[3];
rz(0.66974631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.678858) q[0];
sx q[0];
rz(-0.33927074) q[0];
sx q[0];
rz(-1.6495552) q[0];
rz(0.95343268) q[1];
sx q[1];
rz(-2.0326734) q[1];
sx q[1];
rz(-1.7038201) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9490818) q[0];
sx q[0];
rz(-1.7532945) q[0];
sx q[0];
rz(1.048868) q[0];
rz(-pi) q[1];
rz(1.7069125) q[2];
sx q[2];
rz(-1.7394749) q[2];
sx q[2];
rz(-2.8240734) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9287712) q[1];
sx q[1];
rz(-1.0911687) q[1];
sx q[1];
rz(-2.4368083) q[1];
rz(-0.020047232) q[3];
sx q[3];
rz(-1.123395) q[3];
sx q[3];
rz(-2.945154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.64547223) q[2];
sx q[2];
rz(-0.43473736) q[2];
sx q[2];
rz(-0.17871857) q[2];
rz(-0.86137613) q[3];
sx q[3];
rz(-1.9390315) q[3];
sx q[3];
rz(-2.7698959) q[3];
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
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20539595) q[0];
sx q[0];
rz(-1.2048756) q[0];
sx q[0];
rz(1.0937011) q[0];
rz(2.4049092) q[1];
sx q[1];
rz(-1.8700347) q[1];
sx q[1];
rz(2.0827983) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5537542) q[0];
sx q[0];
rz(-0.61843473) q[0];
sx q[0];
rz(-1.0879602) q[0];
rz(-pi) q[1];
rz(0.40760298) q[2];
sx q[2];
rz(-2.0098445) q[2];
sx q[2];
rz(1.1328732) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.72653786) q[1];
sx q[1];
rz(-1.6931567) q[1];
sx q[1];
rz(1.1916257) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9731673) q[3];
sx q[3];
rz(-2.0109004) q[3];
sx q[3];
rz(-2.6228867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.27292353) q[2];
sx q[2];
rz(-0.69512853) q[2];
sx q[2];
rz(1.6607364) q[2];
rz(0.41040928) q[3];
sx q[3];
rz(-1.7216262) q[3];
sx q[3];
rz(0.15795344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2720298) q[0];
sx q[0];
rz(-2.8700888) q[0];
sx q[0];
rz(-0.29125443) q[0];
rz(2.5323396) q[1];
sx q[1];
rz(-1.6758502) q[1];
sx q[1];
rz(1.4321009) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8490484) q[0];
sx q[0];
rz(-1.8869072) q[0];
sx q[0];
rz(0.51878099) q[0];
rz(2.6860793) q[2];
sx q[2];
rz(-1.3580772) q[2];
sx q[2];
rz(-2.6754975) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.013195136) q[1];
sx q[1];
rz(-2.5880475) q[1];
sx q[1];
rz(-0.08298666) q[1];
rz(-2.0026783) q[3];
sx q[3];
rz(-1.6500041) q[3];
sx q[3];
rz(1.6457886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.46135205) q[2];
sx q[2];
rz(-1.1169008) q[2];
sx q[2];
rz(-1.1414026) q[2];
rz(-1.6067778) q[3];
sx q[3];
rz(-1.1780058) q[3];
sx q[3];
rz(2.6721568) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5948982) q[0];
sx q[0];
rz(-1.6225157) q[0];
sx q[0];
rz(-1.7058104) q[0];
rz(0.36874157) q[1];
sx q[1];
rz(-1.8992966) q[1];
sx q[1];
rz(-0.13175838) q[1];
rz(-1.776406) q[2];
sx q[2];
rz(-1.793135) q[2];
sx q[2];
rz(-1.8829913) q[2];
rz(-2.7183919) q[3];
sx q[3];
rz(-1.7508218) q[3];
sx q[3];
rz(-1.4765061) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
