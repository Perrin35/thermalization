OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1239399) q[0];
sx q[0];
rz(-0.91683638) q[0];
sx q[0];
rz(2.7066948) q[0];
rz(1.5974367) q[1];
sx q[1];
rz(3.6606001) q[1];
sx q[1];
rz(11.984348) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8911444) q[0];
sx q[0];
rz(-1.8316557) q[0];
sx q[0];
rz(-2.9820739) q[0];
rz(-0.97442128) q[2];
sx q[2];
rz(-2.4108464) q[2];
sx q[2];
rz(1.3517018) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4705321) q[1];
sx q[1];
rz(-1.3249505) q[1];
sx q[1];
rz(1.125294) q[1];
rz(-pi) q[2];
x q[2];
rz(1.852096) q[3];
sx q[3];
rz(-1.7419551) q[3];
sx q[3];
rz(0.25118235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.55328289) q[2];
sx q[2];
rz(-1.8841691) q[2];
sx q[2];
rz(2.8498939) q[2];
rz(-0.45082539) q[3];
sx q[3];
rz(-2.8901633) q[3];
sx q[3];
rz(0.60744557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70158231) q[0];
sx q[0];
rz(-0.99047438) q[0];
sx q[0];
rz(1.095358) q[0];
rz(-1.9502684) q[1];
sx q[1];
rz(-2.2215863) q[1];
sx q[1];
rz(-2.2390168) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94287017) q[0];
sx q[0];
rz(-2.1784221) q[0];
sx q[0];
rz(-2.3234419) q[0];
rz(-pi) q[1];
rz(2.2392989) q[2];
sx q[2];
rz(-1.8040787) q[2];
sx q[2];
rz(1.7238613) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.892977) q[1];
sx q[1];
rz(-1.4505523) q[1];
sx q[1];
rz(-0.29386947) q[1];
rz(-pi) q[2];
rz(-0.91530494) q[3];
sx q[3];
rz(-2.4065131) q[3];
sx q[3];
rz(-0.32441586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9537182) q[2];
sx q[2];
rz(-2.1672921) q[2];
sx q[2];
rz(0.11889674) q[2];
rz(-2.4925354) q[3];
sx q[3];
rz(-2.9552112) q[3];
sx q[3];
rz(1.6868748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4897937) q[0];
sx q[0];
rz(-1.2540023) q[0];
sx q[0];
rz(0.5157665) q[0];
rz(-1.0924529) q[1];
sx q[1];
rz(-2.7383883) q[1];
sx q[1];
rz(-1.8663503) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78071252) q[0];
sx q[0];
rz(-3.096252) q[0];
sx q[0];
rz(-0.41098292) q[0];
x q[1];
rz(2.8280806) q[2];
sx q[2];
rz(-1.3517153) q[2];
sx q[2];
rz(0.30922019) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0787492) q[1];
sx q[1];
rz(-0.8991407) q[1];
sx q[1];
rz(0.25154227) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3475092) q[3];
sx q[3];
rz(-1.6424254) q[3];
sx q[3];
rz(2.7290727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.065993) q[2];
sx q[2];
rz(-1.3174572) q[2];
sx q[2];
rz(1.9817748) q[2];
rz(-0.48948151) q[3];
sx q[3];
rz(-1.5533841) q[3];
sx q[3];
rz(-0.71973962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8664261) q[0];
sx q[0];
rz(-0.32298276) q[0];
sx q[0];
rz(0.48165709) q[0];
rz(0.9306759) q[1];
sx q[1];
rz(-1.8447256) q[1];
sx q[1];
rz(-3.1226588) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0136495) q[0];
sx q[0];
rz(-3.0873723) q[0];
sx q[0];
rz(-1.9677866) q[0];
rz(1.5140947) q[2];
sx q[2];
rz(-2.4783274) q[2];
sx q[2];
rz(2.5692289) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0684114) q[1];
sx q[1];
rz(-0.8340237) q[1];
sx q[1];
rz(0.72552105) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3778524) q[3];
sx q[3];
rz(-2.1472048) q[3];
sx q[3];
rz(-1.9595722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.21861741) q[2];
sx q[2];
rz(-2.7851892) q[2];
sx q[2];
rz(-2.3962928) q[2];
rz(-0.37799147) q[3];
sx q[3];
rz(-1.9972921) q[3];
sx q[3];
rz(2.2530344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22076631) q[0];
sx q[0];
rz(-0.31196088) q[0];
sx q[0];
rz(-1.1908603) q[0];
rz(0.4885172) q[1];
sx q[1];
rz(-0.91594511) q[1];
sx q[1];
rz(-0.8078422) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2388326) q[0];
sx q[0];
rz(-1.897398) q[0];
sx q[0];
rz(3.045911) q[0];
rz(-pi) q[1];
rz(-1.135181) q[2];
sx q[2];
rz(-2.6382338) q[2];
sx q[2];
rz(1.2762918) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.96236189) q[1];
sx q[1];
rz(-0.46474248) q[1];
sx q[1];
rz(2.9443113) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2000974) q[3];
sx q[3];
rz(-0.95733023) q[3];
sx q[3];
rz(-2.4870949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0025582) q[2];
sx q[2];
rz(-0.99271861) q[2];
sx q[2];
rz(-2.9774418) q[2];
rz(-0.74603355) q[3];
sx q[3];
rz(-1.7697216) q[3];
sx q[3];
rz(0.073808864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1604851) q[0];
sx q[0];
rz(-1.2657413) q[0];
sx q[0];
rz(-2.4142081) q[0];
rz(2.6630867) q[1];
sx q[1];
rz(-1.452927) q[1];
sx q[1];
rz(0.7368288) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9449483) q[0];
sx q[0];
rz(-1.5153312) q[0];
sx q[0];
rz(1.4329628) q[0];
x q[1];
rz(-1.9987891) q[2];
sx q[2];
rz(-2.8341132) q[2];
sx q[2];
rz(0.82915074) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.068253156) q[1];
sx q[1];
rz(-0.96594772) q[1];
sx q[1];
rz(0.44954919) q[1];
rz(-pi) q[2];
rz(-1.6049006) q[3];
sx q[3];
rz(-0.38785245) q[3];
sx q[3];
rz(2.8451827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9937146) q[2];
sx q[2];
rz(-2.1174049) q[2];
sx q[2];
rz(-0.68515879) q[2];
rz(-2.1327175) q[3];
sx q[3];
rz(-2.4515371) q[3];
sx q[3];
rz(-2.8907997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.93609) q[0];
sx q[0];
rz(-0.093955366) q[0];
sx q[0];
rz(1.093338) q[0];
rz(-3.1178442) q[1];
sx q[1];
rz(-0.7370342) q[1];
sx q[1];
rz(-2.645983) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28195295) q[0];
sx q[0];
rz(-0.067159979) q[0];
sx q[0];
rz(0.18224506) q[0];
rz(-pi) q[1];
rz(0.43735403) q[2];
sx q[2];
rz(-1.7582571) q[2];
sx q[2];
rz(-2.1700493) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.62140761) q[1];
sx q[1];
rz(-1.9193135) q[1];
sx q[1];
rz(-1.3775311) q[1];
rz(-pi) q[2];
rz(0.43066671) q[3];
sx q[3];
rz(-1.1209295) q[3];
sx q[3];
rz(-1.3632185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2451943) q[2];
sx q[2];
rz(-2.3485025) q[2];
sx q[2];
rz(2.8450361) q[2];
rz(-2.631393) q[3];
sx q[3];
rz(-1.5254131) q[3];
sx q[3];
rz(0.27142522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1475913) q[0];
sx q[0];
rz(-2.1864102) q[0];
sx q[0];
rz(1.1543132) q[0];
rz(1.9860024) q[1];
sx q[1];
rz(-2.4842333) q[1];
sx q[1];
rz(-1.4591699) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87932779) q[0];
sx q[0];
rz(-1.5478857) q[0];
sx q[0];
rz(2.4266116) q[0];
x q[1];
rz(-1.5695249) q[2];
sx q[2];
rz(-2.4172815) q[2];
sx q[2];
rz(1.6189761) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7717114) q[1];
sx q[1];
rz(-2.071131) q[1];
sx q[1];
rz(2.9369563) q[1];
rz(-pi) q[2];
rz(-2.1064214) q[3];
sx q[3];
rz(-1.0058306) q[3];
sx q[3];
rz(1.21711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3079754) q[2];
sx q[2];
rz(-0.54215446) q[2];
sx q[2];
rz(1.3758434) q[2];
rz(-0.086325072) q[3];
sx q[3];
rz(-2.3089246) q[3];
sx q[3];
rz(1.8858006) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9432705) q[0];
sx q[0];
rz(-2.4477796) q[0];
sx q[0];
rz(-0.86945239) q[0];
rz(-2.4122639) q[1];
sx q[1];
rz(-2.2980233) q[1];
sx q[1];
rz(-2.9076911) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9651523) q[0];
sx q[0];
rz(-1.5230706) q[0];
sx q[0];
rz(2.7244669) q[0];
rz(-1.4748739) q[2];
sx q[2];
rz(-1.7800343) q[2];
sx q[2];
rz(-2.8411381) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.053596951) q[1];
sx q[1];
rz(-1.1302179) q[1];
sx q[1];
rz(1.7338899) q[1];
rz(-0.85373803) q[3];
sx q[3];
rz(-1.3265564) q[3];
sx q[3];
rz(-1.8239886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.087661155) q[2];
sx q[2];
rz(-1.4940741) q[2];
sx q[2];
rz(1.1762478) q[2];
rz(-2.4704399) q[3];
sx q[3];
rz(-1.2495557) q[3];
sx q[3];
rz(-1.5161071) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7072356) q[0];
sx q[0];
rz(-2.2800627) q[0];
sx q[0];
rz(2.0980515) q[0];
rz(0.097298233) q[1];
sx q[1];
rz(-1.0568591) q[1];
sx q[1];
rz(1.1204488) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9981619) q[0];
sx q[0];
rz(-2.1658346) q[0];
sx q[0];
rz(3.0686994) q[0];
x q[1];
rz(-0.78586833) q[2];
sx q[2];
rz(-2.181727) q[2];
sx q[2];
rz(-0.20139209) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0617025) q[1];
sx q[1];
rz(-1.3933059) q[1];
sx q[1];
rz(1.6203141) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0529641) q[3];
sx q[3];
rz(-1.7381765) q[3];
sx q[3];
rz(2.8958304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6843159) q[2];
sx q[2];
rz(-0.61351675) q[2];
sx q[2];
rz(1.3555917) q[2];
rz(1.5761624) q[3];
sx q[3];
rz(-2.0578945) q[3];
sx q[3];
rz(-2.1255597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9044357) q[0];
sx q[0];
rz(-0.56023993) q[0];
sx q[0];
rz(2.3518363) q[0];
rz(0.65658983) q[1];
sx q[1];
rz(-1.1142535) q[1];
sx q[1];
rz(2.5744892) q[1];
rz(0.91083679) q[2];
sx q[2];
rz(-1.177236) q[2];
sx q[2];
rz(0.81753035) q[2];
rz(0.58224596) q[3];
sx q[3];
rz(-0.34366321) q[3];
sx q[3];
rz(-0.2413078) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
