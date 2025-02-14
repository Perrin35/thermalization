OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3040721) q[0];
sx q[0];
rz(-2.4763595) q[0];
sx q[0];
rz(-1.2256149) q[0];
rz(-0.6839112) q[1];
sx q[1];
rz(3.5332503) q[1];
sx q[1];
rz(11.864301) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4944899) q[0];
sx q[0];
rz(-2.4328874) q[0];
sx q[0];
rz(-2.9604218) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6984387) q[2];
sx q[2];
rz(-2.9753471) q[2];
sx q[2];
rz(0.12285168) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7154123) q[1];
sx q[1];
rz(-0.75580929) q[1];
sx q[1];
rz(2.4017357) q[1];
rz(1.3340215) q[3];
sx q[3];
rz(-0.70565685) q[3];
sx q[3];
rz(1.2232085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2012653) q[2];
sx q[2];
rz(-1.6028812) q[2];
sx q[2];
rz(-0.24844696) q[2];
rz(1.1343608) q[3];
sx q[3];
rz(-3.0076707) q[3];
sx q[3];
rz(1.1321446) q[3];
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
rz(-3.0487173) q[0];
sx q[0];
rz(-2.5889914) q[0];
sx q[0];
rz(-1.5930814) q[0];
rz(2.6379207) q[1];
sx q[1];
rz(-1.2232989) q[1];
sx q[1];
rz(0.37685397) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5232613) q[0];
sx q[0];
rz(-2.0428162) q[0];
sx q[0];
rz(1.5923772) q[0];
rz(-0.84555618) q[2];
sx q[2];
rz(-2.2396002) q[2];
sx q[2];
rz(-1.9921274) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.2959006) q[1];
sx q[1];
rz(-1.0949425) q[1];
sx q[1];
rz(-1.2888539) q[1];
x q[2];
rz(-1.4941856) q[3];
sx q[3];
rz(-1.239068) q[3];
sx q[3];
rz(-2.2463617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.58439955) q[2];
sx q[2];
rz(-0.69584766) q[2];
sx q[2];
rz(-0.86563555) q[2];
rz(-1.4731167) q[3];
sx q[3];
rz(-2.1560463) q[3];
sx q[3];
rz(2.1540811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5675548) q[0];
sx q[0];
rz(-1.2873298) q[0];
sx q[0];
rz(2.0903184) q[0];
rz(1.6846664) q[1];
sx q[1];
rz(-1.8345865) q[1];
sx q[1];
rz(1.5164703) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3932289) q[0];
sx q[0];
rz(-2.3277475) q[0];
sx q[0];
rz(2.7500344) q[0];
rz(2.462492) q[2];
sx q[2];
rz(-1.0179986) q[2];
sx q[2];
rz(-0.9169842) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1976002) q[1];
sx q[1];
rz(-1.9514582) q[1];
sx q[1];
rz(-1.0684816) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6246689) q[3];
sx q[3];
rz(-1.6317211) q[3];
sx q[3];
rz(-2.1046154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7463344) q[2];
sx q[2];
rz(-2.0863775) q[2];
sx q[2];
rz(-2.9193817) q[2];
rz(0.5018417) q[3];
sx q[3];
rz(-1.4072489) q[3];
sx q[3];
rz(2.3327904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.89268452) q[0];
sx q[0];
rz(-2.6539256) q[0];
sx q[0];
rz(-1.1647613) q[0];
rz(-0.89272967) q[1];
sx q[1];
rz(-1.7612532) q[1];
sx q[1];
rz(0.28183118) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.062317693) q[0];
sx q[0];
rz(-1.6910416) q[0];
sx q[0];
rz(-2.6577302) q[0];
rz(-pi) q[1];
rz(-0.71433432) q[2];
sx q[2];
rz(-1.1118982) q[2];
sx q[2];
rz(-2.4752576) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1077256) q[1];
sx q[1];
rz(-1.4078137) q[1];
sx q[1];
rz(3.0708583) q[1];
rz(-pi) q[2];
rz(1.124566) q[3];
sx q[3];
rz(-1.5939624) q[3];
sx q[3];
rz(-1.5678101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0622327) q[2];
sx q[2];
rz(-2.5793109) q[2];
sx q[2];
rz(-2.6083561) q[2];
rz(2.0594635) q[3];
sx q[3];
rz(-0.60492587) q[3];
sx q[3];
rz(0.81348872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.767652) q[0];
sx q[0];
rz(-0.39893183) q[0];
sx q[0];
rz(1.3237413) q[0];
rz(0.81870493) q[1];
sx q[1];
rz(-2.3981514) q[1];
sx q[1];
rz(-2.0735819) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0518347) q[0];
sx q[0];
rz(-1.4115507) q[0];
sx q[0];
rz(1.6510886) q[0];
x q[1];
rz(0.83663656) q[2];
sx q[2];
rz(-1.0435487) q[2];
sx q[2];
rz(0.0076779445) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9449759) q[1];
sx q[1];
rz(-2.2240337) q[1];
sx q[1];
rz(1.7004299) q[1];
x q[2];
rz(2.976494) q[3];
sx q[3];
rz(-1.1037276) q[3];
sx q[3];
rz(1.9593173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9615122) q[2];
sx q[2];
rz(-2.027812) q[2];
sx q[2];
rz(-0.86165825) q[2];
rz(-1.0263475) q[3];
sx q[3];
rz(-1.9828826) q[3];
sx q[3];
rz(-1.1177184) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79477972) q[0];
sx q[0];
rz(-0.81951278) q[0];
sx q[0];
rz(0.89282194) q[0];
rz(-0.18094856) q[1];
sx q[1];
rz(-2.4632958) q[1];
sx q[1];
rz(-0.13042626) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2184046) q[0];
sx q[0];
rz(-1.4124065) q[0];
sx q[0];
rz(0.10938258) q[0];
rz(-pi) q[1];
rz(2.4685235) q[2];
sx q[2];
rz(-2.4166738) q[2];
sx q[2];
rz(2.6672305) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.45589009) q[1];
sx q[1];
rz(-1.8079385) q[1];
sx q[1];
rz(0.030045998) q[1];
rz(-pi) q[2];
rz(2.942286) q[3];
sx q[3];
rz(-0.79108566) q[3];
sx q[3];
rz(-1.6419322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.22455198) q[2];
sx q[2];
rz(-1.2265393) q[2];
sx q[2];
rz(-0.87635931) q[2];
rz(1.2419491) q[3];
sx q[3];
rz(-2.2010937) q[3];
sx q[3];
rz(2.131264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7414311) q[0];
sx q[0];
rz(-0.09859666) q[0];
sx q[0];
rz(0.36369351) q[0];
rz(-0.74186507) q[1];
sx q[1];
rz(-1.0262841) q[1];
sx q[1];
rz(2.738764) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3198217) q[0];
sx q[0];
rz(-1.1603702) q[0];
sx q[0];
rz(-0.26217802) q[0];
rz(-pi) q[1];
rz(2.3260457) q[2];
sx q[2];
rz(-2.6167653) q[2];
sx q[2];
rz(1.8458402) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0707222) q[1];
sx q[1];
rz(-1.6477343) q[1];
sx q[1];
rz(1.6926195) q[1];
rz(-pi) q[2];
rz(2.7807523) q[3];
sx q[3];
rz(-0.66698241) q[3];
sx q[3];
rz(2.9426394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7824629) q[2];
sx q[2];
rz(-0.368258) q[2];
sx q[2];
rz(-1.5472319) q[2];
rz(2.3063229) q[3];
sx q[3];
rz(-2.1850977) q[3];
sx q[3];
rz(-0.48867759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7816724) q[0];
sx q[0];
rz(-1.8052354) q[0];
sx q[0];
rz(-0.51505995) q[0];
rz(-1.4393648) q[1];
sx q[1];
rz(-1.2934338) q[1];
sx q[1];
rz(-0.39915592) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3418509) q[0];
sx q[0];
rz(-3.0106595) q[0];
sx q[0];
rz(-3.127742) q[0];
rz(-pi) q[1];
rz(2.299526) q[2];
sx q[2];
rz(-1.7532323) q[2];
sx q[2];
rz(1.6853756) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6281575) q[1];
sx q[1];
rz(-1.6353459) q[1];
sx q[1];
rz(-1.8316395) q[1];
rz(-pi) q[2];
rz(3.0759381) q[3];
sx q[3];
rz(-1.5695509) q[3];
sx q[3];
rz(-0.37855442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9814375) q[2];
sx q[2];
rz(-1.258054) q[2];
sx q[2];
rz(-2.0951648) q[2];
rz(-2.4753172) q[3];
sx q[3];
rz(-0.25686887) q[3];
sx q[3];
rz(2.090914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(2.1995354) q[0];
sx q[0];
rz(-3.0362447) q[0];
sx q[0];
rz(-2.5355205) q[0];
rz(-0.82707682) q[1];
sx q[1];
rz(-1.8111633) q[1];
sx q[1];
rz(-1.3333295) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24902835) q[0];
sx q[0];
rz(-1.4672534) q[0];
sx q[0];
rz(2.7825481) q[0];
rz(-pi) q[1];
rz(2.3982348) q[2];
sx q[2];
rz(-1.9089235) q[2];
sx q[2];
rz(-1.0184792) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.575653) q[1];
sx q[1];
rz(-2.2350532) q[1];
sx q[1];
rz(-0.14845005) q[1];
x q[2];
rz(-0.17154947) q[3];
sx q[3];
rz(-2.6772873) q[3];
sx q[3];
rz(2.7363079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0102319) q[2];
sx q[2];
rz(-2.2183245) q[2];
sx q[2];
rz(0.38491797) q[2];
rz(-2.5108003) q[3];
sx q[3];
rz(-1.2814949) q[3];
sx q[3];
rz(1.7990641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.887562) q[0];
sx q[0];
rz(-1.7429202) q[0];
sx q[0];
rz(1.7015464) q[0];
rz(0.74132672) q[1];
sx q[1];
rz(-1.7404375) q[1];
sx q[1];
rz(2.8828566) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3869303) q[0];
sx q[0];
rz(-1.5748236) q[0];
sx q[0];
rz(1.8167102) q[0];
rz(-0.3157987) q[2];
sx q[2];
rz(-2.0673896) q[2];
sx q[2];
rz(1.1441355) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8347296) q[1];
sx q[1];
rz(-2.266699) q[1];
sx q[1];
rz(-0.45180288) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1734937) q[3];
sx q[3];
rz(-2.3081995) q[3];
sx q[3];
rz(1.6034338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4197293) q[2];
sx q[2];
rz(-1.1121007) q[2];
sx q[2];
rz(1.8224243) q[2];
rz(-0.81926528) q[3];
sx q[3];
rz(-2.0115439) q[3];
sx q[3];
rz(0.54660249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7405613) q[0];
sx q[0];
rz(-2.2673829) q[0];
sx q[0];
rz(-0.54947214) q[0];
rz(1.7026547) q[1];
sx q[1];
rz(-2.4733652) q[1];
sx q[1];
rz(0.53818902) q[1];
rz(2.2808711) q[2];
sx q[2];
rz(-1.5624983) q[2];
sx q[2];
rz(2.7539243) q[2];
rz(1.0092173) q[3];
sx q[3];
rz(-1.781096) q[3];
sx q[3];
rz(-2.5566035) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
