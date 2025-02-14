OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.34484997) q[0];
sx q[0];
rz(-0.27422187) q[0];
sx q[0];
rz(0.56871498) q[0];
rz(-1.93058) q[1];
sx q[1];
rz(-0.99744263) q[1];
sx q[1];
rz(-2.8740191) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5481603) q[0];
sx q[0];
rz(-0.58565564) q[0];
sx q[0];
rz(-0.36958739) q[0];
rz(-3.028439) q[2];
sx q[2];
rz(-1.238816) q[2];
sx q[2];
rz(-1.5205713) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.97155826) q[1];
sx q[1];
rz(-1.5624678) q[1];
sx q[1];
rz(0.002408601) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4300284) q[3];
sx q[3];
rz(-1.5091617) q[3];
sx q[3];
rz(0.34256645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.25207511) q[2];
sx q[2];
rz(-2.0740261) q[2];
sx q[2];
rz(-2.0102823) q[2];
rz(0.45025292) q[3];
sx q[3];
rz(-2.4501652) q[3];
sx q[3];
rz(-2.1972307) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64304072) q[0];
sx q[0];
rz(-2.2096071) q[0];
sx q[0];
rz(1.407628) q[0];
rz(-2.6990926) q[1];
sx q[1];
rz(-1.4235556) q[1];
sx q[1];
rz(-2.5462467) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2059518) q[0];
sx q[0];
rz(-1.3229587) q[0];
sx q[0];
rz(1.3252844) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6518455) q[2];
sx q[2];
rz(-2.1648266) q[2];
sx q[2];
rz(2.5777532) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1658948) q[1];
sx q[1];
rz(-1.5656316) q[1];
sx q[1];
rz(-2.7126606) q[1];
rz(-pi) q[2];
rz(0.54259681) q[3];
sx q[3];
rz(-0.68477453) q[3];
sx q[3];
rz(-2.7921576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2596316) q[2];
sx q[2];
rz(-0.61820784) q[2];
sx q[2];
rz(-2.372443) q[2];
rz(1.361557) q[3];
sx q[3];
rz(-0.96746126) q[3];
sx q[3];
rz(0.59534016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24460569) q[0];
sx q[0];
rz(-0.91490442) q[0];
sx q[0];
rz(-2.8047674) q[0];
rz(-1.7103851) q[1];
sx q[1];
rz(-2.2957048) q[1];
sx q[1];
rz(-3.0444042) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021193209) q[0];
sx q[0];
rz(-1.0836856) q[0];
sx q[0];
rz(0.22986408) q[0];
rz(-pi) q[1];
rz(1.4732185) q[2];
sx q[2];
rz(-1.9505902) q[2];
sx q[2];
rz(-1.1792091) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1657451) q[1];
sx q[1];
rz(-2.7677508) q[1];
sx q[1];
rz(-2.7553619) q[1];
x q[2];
rz(2.187927) q[3];
sx q[3];
rz(-2.3535471) q[3];
sx q[3];
rz(-2.5833481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.97239697) q[2];
sx q[2];
rz(-1.7652067) q[2];
sx q[2];
rz(-2.3248559) q[2];
rz(-2.2190602) q[3];
sx q[3];
rz(-2.7042992) q[3];
sx q[3];
rz(1.9492662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25794491) q[0];
sx q[0];
rz(-1.3379931) q[0];
sx q[0];
rz(-2.2747967) q[0];
rz(1.2358933) q[1];
sx q[1];
rz(-1.1150603) q[1];
sx q[1];
rz(2.8996276) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4075027) q[0];
sx q[0];
rz(-0.76341141) q[0];
sx q[0];
rz(1.1795189) q[0];
x q[1];
rz(-3.0680806) q[2];
sx q[2];
rz(-1.4097555) q[2];
sx q[2];
rz(2.0025776) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0556524) q[1];
sx q[1];
rz(-1.4341518) q[1];
sx q[1];
rz(0.40811347) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6446043) q[3];
sx q[3];
rz(-1.6145633) q[3];
sx q[3];
rz(0.62326335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.26291004) q[2];
sx q[2];
rz(-1.5107369) q[2];
sx q[2];
rz(2.6061457) q[2];
rz(0.49324909) q[3];
sx q[3];
rz(-0.89087629) q[3];
sx q[3];
rz(2.6935327) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.066147476) q[0];
sx q[0];
rz(-2.9904521) q[0];
sx q[0];
rz(-2.2976663) q[0];
rz(-1.4752202) q[1];
sx q[1];
rz(-1.6553144) q[1];
sx q[1];
rz(0.39316887) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61680865) q[0];
sx q[0];
rz(-2.3032585) q[0];
sx q[0];
rz(-0.34466593) q[0];
rz(-2.5049823) q[2];
sx q[2];
rz(-1.746897) q[2];
sx q[2];
rz(-2.2409093) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0991126) q[1];
sx q[1];
rz(-2.5951324) q[1];
sx q[1];
rz(0.55621712) q[1];
rz(-1.0814905) q[3];
sx q[3];
rz(-0.84546465) q[3];
sx q[3];
rz(0.368834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7273442) q[2];
sx q[2];
rz(-2.9372637) q[2];
sx q[2];
rz(0.3981398) q[2];
rz(-2.5967755) q[3];
sx q[3];
rz(-2.3530493) q[3];
sx q[3];
rz(1.8383693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9141465) q[0];
sx q[0];
rz(-1.0108203) q[0];
sx q[0];
rz(0.8557125) q[0];
rz(-2.4489467) q[1];
sx q[1];
rz(-0.99450642) q[1];
sx q[1];
rz(-2.8725502) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1616283) q[0];
sx q[0];
rz(-1.3792896) q[0];
sx q[0];
rz(3.1167517) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5085717) q[2];
sx q[2];
rz(-0.50743689) q[2];
sx q[2];
rz(-0.98558805) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0768041) q[1];
sx q[1];
rz(-2.2858983) q[1];
sx q[1];
rz(0.24311693) q[1];
rz(1.3790491) q[3];
sx q[3];
rz(-1.4095528) q[3];
sx q[3];
rz(1.92056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1143703) q[2];
sx q[2];
rz(-1.3193069) q[2];
sx q[2];
rz(-3.031292) q[2];
rz(2.2731764) q[3];
sx q[3];
rz(-1.9649558) q[3];
sx q[3];
rz(-1.7051914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7480943) q[0];
sx q[0];
rz(-2.6415934) q[0];
sx q[0];
rz(-0.86135832) q[0];
rz(1.512108) q[1];
sx q[1];
rz(-1.6810828) q[1];
sx q[1];
rz(-2.3177564) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55144009) q[0];
sx q[0];
rz(-0.9853029) q[0];
sx q[0];
rz(-3.0072938) q[0];
rz(-pi) q[1];
rz(1.9806978) q[2];
sx q[2];
rz(-0.86386743) q[2];
sx q[2];
rz(0.17690578) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.87005918) q[1];
sx q[1];
rz(-2.4932269) q[1];
sx q[1];
rz(-0.35761498) q[1];
rz(-1.8855612) q[3];
sx q[3];
rz(-0.94292313) q[3];
sx q[3];
rz(1.0726269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6350101) q[2];
sx q[2];
rz(-2.625605) q[2];
sx q[2];
rz(0.79279509) q[2];
rz(0.005216287) q[3];
sx q[3];
rz(-2.3514533) q[3];
sx q[3];
rz(-1.691157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1714627) q[0];
sx q[0];
rz(-1.1928394) q[0];
sx q[0];
rz(-0.18950732) q[0];
rz(0.7729404) q[1];
sx q[1];
rz(-2.645292) q[1];
sx q[1];
rz(2.5453087) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8533289) q[0];
sx q[0];
rz(-1.402308) q[0];
sx q[0];
rz(-0.69126076) q[0];
rz(0.79769602) q[2];
sx q[2];
rz(-2.0828649) q[2];
sx q[2];
rz(1.860581) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.3459732) q[1];
sx q[1];
rz(-1.4782149) q[1];
sx q[1];
rz(3.0365192) q[1];
rz(-pi) q[2];
rz(2.2828388) q[3];
sx q[3];
rz(-1.158445) q[3];
sx q[3];
rz(-1.2996197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.72558713) q[2];
sx q[2];
rz(-0.013898762) q[2];
sx q[2];
rz(1.928398) q[2];
rz(0.98617918) q[3];
sx q[3];
rz(-1.4158019) q[3];
sx q[3];
rz(1.8825611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.1213433) q[0];
sx q[0];
rz(-1.5859402) q[0];
sx q[0];
rz(-2.0315309) q[0];
rz(0.32866651) q[1];
sx q[1];
rz(-1.5549436) q[1];
sx q[1];
rz(1.2967671) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4121108) q[0];
sx q[0];
rz(-2.4769432) q[0];
sx q[0];
rz(-1.7230711) q[0];
rz(-pi) q[1];
x q[1];
rz(0.062053238) q[2];
sx q[2];
rz(-1.1750571) q[2];
sx q[2];
rz(1.6153276) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4420524) q[1];
sx q[1];
rz(-1.6332383) q[1];
sx q[1];
rz(-2.6477835) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66547439) q[3];
sx q[3];
rz(-1.7413119) q[3];
sx q[3];
rz(-1.0118893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0729596) q[2];
sx q[2];
rz(-2.1559842) q[2];
sx q[2];
rz(3.0375321) q[2];
rz(2.0067298) q[3];
sx q[3];
rz(-1.3645423) q[3];
sx q[3];
rz(0.22741905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65570152) q[0];
sx q[0];
rz(-0.98485297) q[0];
sx q[0];
rz(-0.73053288) q[0];
rz(2.5841374) q[1];
sx q[1];
rz(-1.971222) q[1];
sx q[1];
rz(-2.7117859) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24048478) q[0];
sx q[0];
rz(-0.76178023) q[0];
sx q[0];
rz(3.1408589) q[0];
rz(2.7676959) q[2];
sx q[2];
rz(-1.3329525) q[2];
sx q[2];
rz(1.6872981) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0954206) q[1];
sx q[1];
rz(-0.37994994) q[1];
sx q[1];
rz(-1.7871961) q[1];
rz(-pi) q[2];
rz(-3.0176875) q[3];
sx q[3];
rz(-1.2031816) q[3];
sx q[3];
rz(-2.1218561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.80959117) q[2];
sx q[2];
rz(-2.1198699) q[2];
sx q[2];
rz(-0.29279718) q[2];
rz(3.0012567) q[3];
sx q[3];
rz(-2.1122746) q[3];
sx q[3];
rz(-0.50104195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3705227) q[0];
sx q[0];
rz(-1.4650383) q[0];
sx q[0];
rz(-2.9248206) q[0];
rz(-0.71939214) q[1];
sx q[1];
rz(-1.8665301) q[1];
sx q[1];
rz(-2.9449609) q[1];
rz(-1.2540631) q[2];
sx q[2];
rz(-2.6120196) q[2];
sx q[2];
rz(2.0092464) q[2];
rz(-1.9543129) q[3];
sx q[3];
rz(-0.34579943) q[3];
sx q[3];
rz(-0.75172206) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
