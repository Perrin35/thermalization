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
rz(-1.5440829) q[0];
sx q[0];
rz(-1.7680327) q[0];
sx q[0];
rz(-1.0184259) q[0];
rz(2.6234558) q[1];
sx q[1];
rz(-0.75704804) q[1];
sx q[1];
rz(2.5098324) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9874532) q[0];
sx q[0];
rz(-1.69603) q[0];
sx q[0];
rz(1.3628598) q[0];
rz(-pi) q[1];
rz(1.2680797) q[2];
sx q[2];
rz(-1.5470501) q[2];
sx q[2];
rz(2.077092) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85447299) q[1];
sx q[1];
rz(-1.2911671) q[1];
sx q[1];
rz(-0.92052144) q[1];
rz(-pi) q[2];
rz(1.3034091) q[3];
sx q[3];
rz(-0.33743706) q[3];
sx q[3];
rz(0.22663675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6112001) q[2];
sx q[2];
rz(-2.7860614) q[2];
sx q[2];
rz(1.6543039) q[2];
rz(-2.2330331) q[3];
sx q[3];
rz(-2.9049951) q[3];
sx q[3];
rz(-2.1366185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0205883) q[0];
sx q[0];
rz(-1.7089184) q[0];
sx q[0];
rz(-2.9793136) q[0];
rz(0.24457112) q[1];
sx q[1];
rz(-1.2030315) q[1];
sx q[1];
rz(2.8499106) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8732018) q[0];
sx q[0];
rz(-1.2884198) q[0];
sx q[0];
rz(-1.7712405) q[0];
x q[1];
rz(2.9180718) q[2];
sx q[2];
rz(-1.2940027) q[2];
sx q[2];
rz(1.8794488) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4710992) q[1];
sx q[1];
rz(-2.3355977) q[1];
sx q[1];
rz(-2.8189895) q[1];
rz(2.347885) q[3];
sx q[3];
rz(-2.5849403) q[3];
sx q[3];
rz(1.8451549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4633999) q[2];
sx q[2];
rz(-0.86898154) q[2];
sx q[2];
rz(2.0657067) q[2];
rz(1.2470657) q[3];
sx q[3];
rz(-1.5445292) q[3];
sx q[3];
rz(1.4472848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4397044) q[0];
sx q[0];
rz(-2.0151558) q[0];
sx q[0];
rz(-2.9578748) q[0];
rz(1.5048997) q[1];
sx q[1];
rz(-0.74438649) q[1];
sx q[1];
rz(0.16477975) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8863109) q[0];
sx q[0];
rz(-1.5802529) q[0];
sx q[0];
rz(2.1566118) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7585331) q[2];
sx q[2];
rz(-0.7128517) q[2];
sx q[2];
rz(-2.0682316) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2307869) q[1];
sx q[1];
rz(-1.1890829) q[1];
sx q[1];
rz(2.4934105) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0298877) q[3];
sx q[3];
rz(-2.5383484) q[3];
sx q[3];
rz(1.3689976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.034721) q[2];
sx q[2];
rz(-1.160459) q[2];
sx q[2];
rz(1.1064233) q[2];
rz(2.4404081) q[3];
sx q[3];
rz(-1.8132352) q[3];
sx q[3];
rz(-2.6760694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3727386) q[0];
sx q[0];
rz(-1.1356069) q[0];
sx q[0];
rz(0.94451529) q[0];
rz(-1.6230029) q[1];
sx q[1];
rz(-2.1606162) q[1];
sx q[1];
rz(1.6015582) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042451579) q[0];
sx q[0];
rz(-0.97223982) q[0];
sx q[0];
rz(0.70229806) q[0];
rz(-2.3643199) q[2];
sx q[2];
rz(-2.9748355) q[2];
sx q[2];
rz(-0.47698944) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2112) q[1];
sx q[1];
rz(-0.89387776) q[1];
sx q[1];
rz(-0.48447123) q[1];
rz(-pi) q[2];
x q[2];
rz(0.047961162) q[3];
sx q[3];
rz(-0.75100079) q[3];
sx q[3];
rz(2.6329071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.010178415) q[2];
sx q[2];
rz(-2.5203036) q[2];
sx q[2];
rz(0.16769257) q[2];
rz(3.0883279) q[3];
sx q[3];
rz(-2.0361418) q[3];
sx q[3];
rz(1.3767327) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57529706) q[0];
sx q[0];
rz(-3.0360041) q[0];
sx q[0];
rz(-2.8422624) q[0];
rz(0.37295595) q[1];
sx q[1];
rz(-1.8920218) q[1];
sx q[1];
rz(1.6711055) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9540967) q[0];
sx q[0];
rz(-2.4135655) q[0];
sx q[0];
rz(2.4690829) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46385455) q[2];
sx q[2];
rz(-1.4632483) q[2];
sx q[2];
rz(-2.7908418) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8229482) q[1];
sx q[1];
rz(-0.80893436) q[1];
sx q[1];
rz(-0.72973324) q[1];
rz(-pi) q[2];
rz(0.17584189) q[3];
sx q[3];
rz(-2.1851563) q[3];
sx q[3];
rz(1.9906438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2028929) q[2];
sx q[2];
rz(-2.4788269) q[2];
sx q[2];
rz(-2.0311484) q[2];
rz(-0.86197305) q[3];
sx q[3];
rz(-1.5098666) q[3];
sx q[3];
rz(-1.8814258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92457572) q[0];
sx q[0];
rz(-0.68704263) q[0];
sx q[0];
rz(0.4050912) q[0];
rz(-0.336126) q[1];
sx q[1];
rz(-1.7205709) q[1];
sx q[1];
rz(-2.1902671) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38430957) q[0];
sx q[0];
rz(-1.3157789) q[0];
sx q[0];
rz(-2.7268305) q[0];
rz(-pi) q[1];
rz(0.73619618) q[2];
sx q[2];
rz(-1.412782) q[2];
sx q[2];
rz(-0.73537961) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0266708) q[1];
sx q[1];
rz(-1.7975866) q[1];
sx q[1];
rz(0.58057745) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5441465) q[3];
sx q[3];
rz(-2.8077112) q[3];
sx q[3];
rz(-0.88874879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.090791) q[2];
sx q[2];
rz(-1.7369221) q[2];
sx q[2];
rz(-0.23078272) q[2];
rz(-2.9595621) q[3];
sx q[3];
rz(-2.5735276) q[3];
sx q[3];
rz(0.52869421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7886605) q[0];
sx q[0];
rz(-0.039529888) q[0];
sx q[0];
rz(0.93210644) q[0];
rz(2.508029) q[1];
sx q[1];
rz(-2.0220058) q[1];
sx q[1];
rz(1.8036141) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9553878) q[0];
sx q[0];
rz(-3.0078631) q[0];
sx q[0];
rz(-0.9132847) q[0];
rz(2.7510178) q[2];
sx q[2];
rz(-2.8468067) q[2];
sx q[2];
rz(2.8975671) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2436314) q[1];
sx q[1];
rz(-0.96968953) q[1];
sx q[1];
rz(-0.48378418) q[1];
x q[2];
rz(-0.49832817) q[3];
sx q[3];
rz(-1.5428203) q[3];
sx q[3];
rz(1.8701815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6771217) q[2];
sx q[2];
rz(-1.2268343) q[2];
sx q[2];
rz(1.202549) q[2];
rz(-2.7770212) q[3];
sx q[3];
rz(-2.4489844) q[3];
sx q[3];
rz(-0.83474368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6637591) q[0];
sx q[0];
rz(-2.5680225) q[0];
sx q[0];
rz(-2.9649576) q[0];
rz(2.7015576) q[1];
sx q[1];
rz(-1.1853848) q[1];
sx q[1];
rz(-2.1766591) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4266708) q[0];
sx q[0];
rz(-1.3765125) q[0];
sx q[0];
rz(1.4358836) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60243269) q[2];
sx q[2];
rz(-1.1083853) q[2];
sx q[2];
rz(-1.7830069) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7495887) q[1];
sx q[1];
rz(-2.4676305) q[1];
sx q[1];
rz(2.1753009) q[1];
rz(-pi) q[2];
rz(0.075841622) q[3];
sx q[3];
rz(-1.8245398) q[3];
sx q[3];
rz(0.50204078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.21260103) q[2];
sx q[2];
rz(-1.6182199) q[2];
sx q[2];
rz(-3.1316481) q[2];
rz(3.0250004) q[3];
sx q[3];
rz(-2.8023585) q[3];
sx q[3];
rz(-2.5998083) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5803439) q[0];
sx q[0];
rz(-1.9191701) q[0];
sx q[0];
rz(2.5196581) q[0];
rz(-1.4382582) q[1];
sx q[1];
rz(-2.56918) q[1];
sx q[1];
rz(-0.66351801) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2292762) q[0];
sx q[0];
rz(-1.352542) q[0];
sx q[0];
rz(0.7383607) q[0];
rz(1.8542669) q[2];
sx q[2];
rz(-1.6890235) q[2];
sx q[2];
rz(2.9663939) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.92163699) q[1];
sx q[1];
rz(-1.2531279) q[1];
sx q[1];
rz(2.2457473) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8597707) q[3];
sx q[3];
rz(-1.4803866) q[3];
sx q[3];
rz(-2.2304931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5180987) q[2];
sx q[2];
rz(-1.3891209) q[2];
sx q[2];
rz(2.7790879) q[2];
rz(-2.6643961) q[3];
sx q[3];
rz(-2.0878744) q[3];
sx q[3];
rz(2.0188735) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0839194) q[0];
sx q[0];
rz(-0.73129439) q[0];
sx q[0];
rz(-0.90173632) q[0];
rz(-2.4665191) q[1];
sx q[1];
rz(-2.156064) q[1];
sx q[1];
rz(2.9973082) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3659199) q[0];
sx q[0];
rz(-0.36078851) q[0];
sx q[0];
rz(1.3350639) q[0];
x q[1];
rz(-0.4298019) q[2];
sx q[2];
rz(-0.24082213) q[2];
sx q[2];
rz(2.4483829) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.044202494) q[1];
sx q[1];
rz(-0.69141885) q[1];
sx q[1];
rz(-0.39908646) q[1];
rz(-pi) q[2];
rz(0.38548174) q[3];
sx q[3];
rz(-2.6691438) q[3];
sx q[3];
rz(-2.2466898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5054063) q[2];
sx q[2];
rz(-1.9289086) q[2];
sx q[2];
rz(-2.9409161) q[2];
rz(-1.1096654) q[3];
sx q[3];
rz(-2.6789013) q[3];
sx q[3];
rz(-2.5698575) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5889482) q[0];
sx q[0];
rz(-2.1727967) q[0];
sx q[0];
rz(-1.1217242) q[0];
rz(-0.48925346) q[1];
sx q[1];
rz(-1.6901292) q[1];
sx q[1];
rz(2.1314175) q[1];
rz(1.1560925) q[2];
sx q[2];
rz(-1.4857875) q[2];
sx q[2];
rz(-0.93405741) q[2];
rz(-0.060754178) q[3];
sx q[3];
rz(-2.698632) q[3];
sx q[3];
rz(-2.5198577) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
