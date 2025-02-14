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
rz(-2.0105536) q[0];
sx q[0];
rz(3.717489) q[0];
sx q[0];
rz(5.1562638) q[0];
rz(2.1908886) q[1];
sx q[1];
rz(3.7428441) q[1];
sx q[1];
rz(9.7584702) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93984825) q[0];
sx q[0];
rz(-0.44838312) q[0];
sx q[0];
rz(-1.9223619) q[0];
rz(2.0804685) q[2];
sx q[2];
rz(-2.6303419) q[2];
sx q[2];
rz(-2.0927657) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5383) q[1];
sx q[1];
rz(-1.913402) q[1];
sx q[1];
rz(-1.3450772) q[1];
rz(1.8435616) q[3];
sx q[3];
rz(-0.78744167) q[3];
sx q[3];
rz(1.8026601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2600962) q[2];
sx q[2];
rz(-1.0755971) q[2];
sx q[2];
rz(1.0502226) q[2];
rz(-2.8460734) q[3];
sx q[3];
rz(-2.3727356) q[3];
sx q[3];
rz(-1.7723005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1693901) q[0];
sx q[0];
rz(-1.0193595) q[0];
sx q[0];
rz(-0.66725677) q[0];
rz(-2.9908906) q[1];
sx q[1];
rz(-1.0548016) q[1];
sx q[1];
rz(2.8299455) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3903605) q[0];
sx q[0];
rz(-1.5003107) q[0];
sx q[0];
rz(-2.0930392) q[0];
rz(-0.30782757) q[2];
sx q[2];
rz(-1.1500949) q[2];
sx q[2];
rz(1.6604648) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.63285349) q[1];
sx q[1];
rz(-1.3876186) q[1];
sx q[1];
rz(1.4289209) q[1];
x q[2];
rz(-0.39763173) q[3];
sx q[3];
rz(-1.6379698) q[3];
sx q[3];
rz(0.45715082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.47505891) q[2];
sx q[2];
rz(-0.89343137) q[2];
sx q[2];
rz(0.7676355) q[2];
rz(2.660699) q[3];
sx q[3];
rz(-2.3700263) q[3];
sx q[3];
rz(-1.5546999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2956706) q[0];
sx q[0];
rz(-1.6039811) q[0];
sx q[0];
rz(0.77958244) q[0];
rz(0.37173158) q[1];
sx q[1];
rz(-2.1340243) q[1];
sx q[1];
rz(-2.380611) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62255732) q[0];
sx q[0];
rz(-1.2815968) q[0];
sx q[0];
rz(1.9054806) q[0];
rz(2.2566608) q[2];
sx q[2];
rz(-1.2239309) q[2];
sx q[2];
rz(2.318813) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0798988) q[1];
sx q[1];
rz(-1.4053646) q[1];
sx q[1];
rz(2.295619) q[1];
rz(-pi) q[2];
rz(-1.5651302) q[3];
sx q[3];
rz(-0.71810371) q[3];
sx q[3];
rz(2.6942962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.80428213) q[2];
sx q[2];
rz(-0.90039841) q[2];
sx q[2];
rz(-2.6864181) q[2];
rz(-2.4496487) q[3];
sx q[3];
rz(-0.71440905) q[3];
sx q[3];
rz(-0.80879319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-2.1728425) q[0];
sx q[0];
rz(-1.7946365) q[0];
sx q[0];
rz(2.8853048) q[0];
rz(1.434727) q[1];
sx q[1];
rz(-1.4204357) q[1];
sx q[1];
rz(-3.0228379) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1349604) q[0];
sx q[0];
rz(-1.7259571) q[0];
sx q[0];
rz(-2.6895903) q[0];
rz(3.0135113) q[2];
sx q[2];
rz(-1.2799147) q[2];
sx q[2];
rz(-1.1156991) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2251216) q[1];
sx q[1];
rz(-0.82092228) q[1];
sx q[1];
rz(-0.31912843) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3865755) q[3];
sx q[3];
rz(-1.8510071) q[3];
sx q[3];
rz(1.8108944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8541096) q[2];
sx q[2];
rz(-0.27366769) q[2];
sx q[2];
rz(1.0556833) q[2];
rz(1.6915197) q[3];
sx q[3];
rz(-1.4472716) q[3];
sx q[3];
rz(-0.84901866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56181041) q[0];
sx q[0];
rz(-0.17794839) q[0];
sx q[0];
rz(-1.6135038) q[0];
rz(-1.1771419) q[1];
sx q[1];
rz(-1.8714995) q[1];
sx q[1];
rz(1.5632163) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6230683) q[0];
sx q[0];
rz(-2.1954064) q[0];
sx q[0];
rz(1.4853857) q[0];
rz(-2.4776978) q[2];
sx q[2];
rz(-0.50506578) q[2];
sx q[2];
rz(0.31873825) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8635378) q[1];
sx q[1];
rz(-1.7124933) q[1];
sx q[1];
rz(0.63512953) q[1];
rz(-pi) q[2];
rz(1.568119) q[3];
sx q[3];
rz(-1.9802046) q[3];
sx q[3];
rz(0.97312991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0057808) q[2];
sx q[2];
rz(-1.9848738) q[2];
sx q[2];
rz(1.8394252) q[2];
rz(2.8016413) q[3];
sx q[3];
rz(-2.3348742) q[3];
sx q[3];
rz(0.66238856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0078916773) q[0];
sx q[0];
rz(-1.5444724) q[0];
sx q[0];
rz(1.9418465) q[0];
rz(-1.7802736) q[1];
sx q[1];
rz(-2.168455) q[1];
sx q[1];
rz(-0.57788411) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21450522) q[0];
sx q[0];
rz(-1.5772235) q[0];
sx q[0];
rz(-1.5570197) q[0];
rz(-pi) q[1];
rz(-0.33249929) q[2];
sx q[2];
rz(-2.6869171) q[2];
sx q[2];
rz(-2.5408059) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.51398811) q[1];
sx q[1];
rz(-1.1179233) q[1];
sx q[1];
rz(1.7220366) q[1];
x q[2];
rz(-1.6541566) q[3];
sx q[3];
rz(-2.8417086) q[3];
sx q[3];
rz(-0.56902992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.63152385) q[2];
sx q[2];
rz(-2.5670467) q[2];
sx q[2];
rz(-2.6681382) q[2];
rz(1.8691285) q[3];
sx q[3];
rz(-1.748964) q[3];
sx q[3];
rz(-0.22245358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1943787) q[0];
sx q[0];
rz(-0.65776238) q[0];
sx q[0];
rz(2.8084602) q[0];
rz(-2.9130452) q[1];
sx q[1];
rz(-1.3401778) q[1];
sx q[1];
rz(2.5659335) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7299907) q[0];
sx q[0];
rz(-1.4648998) q[0];
sx q[0];
rz(-2.3292755) q[0];
rz(-pi) q[1];
rz(2.6645053) q[2];
sx q[2];
rz(-1.2735575) q[2];
sx q[2];
rz(-1.5971668) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.043283431) q[1];
sx q[1];
rz(-2.7554338) q[1];
sx q[1];
rz(2.1257867) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6778464) q[3];
sx q[3];
rz(-0.52846842) q[3];
sx q[3];
rz(-1.0077882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8754742) q[2];
sx q[2];
rz(-0.55142752) q[2];
sx q[2];
rz(1.4761338) q[2];
rz(2.0406145) q[3];
sx q[3];
rz(-2.0783547) q[3];
sx q[3];
rz(0.5180009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.5440893) q[0];
sx q[0];
rz(-0.88633716) q[0];
sx q[0];
rz(-2.8344717) q[0];
rz(2.8111474) q[1];
sx q[1];
rz(-1.5978866) q[1];
sx q[1];
rz(-1.7764567) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48956051) q[0];
sx q[0];
rz(-0.083261641) q[0];
sx q[0];
rz(0.21647446) q[0];
rz(-pi) q[1];
rz(-0.48319419) q[2];
sx q[2];
rz(-2.7527134) q[2];
sx q[2];
rz(2.4318397) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4933984) q[1];
sx q[1];
rz(-1.9983564) q[1];
sx q[1];
rz(2.8570064) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4151814) q[3];
sx q[3];
rz(-0.74374108) q[3];
sx q[3];
rz(1.2228325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8660628) q[2];
sx q[2];
rz(-0.58774647) q[2];
sx q[2];
rz(2.8928939) q[2];
rz(1.9281049) q[3];
sx q[3];
rz(-1.6971734) q[3];
sx q[3];
rz(0.061633751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1520749) q[0];
sx q[0];
rz(-0.9062506) q[0];
sx q[0];
rz(-0.686598) q[0];
rz(-2.0286512) q[1];
sx q[1];
rz(-0.80250347) q[1];
sx q[1];
rz(0.31563219) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0811809) q[0];
sx q[0];
rz(-0.32627772) q[0];
sx q[0];
rz(1.7148036) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8191843) q[2];
sx q[2];
rz(-0.82652107) q[2];
sx q[2];
rz(-2.7374817) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3179847) q[1];
sx q[1];
rz(-0.66901842) q[1];
sx q[1];
rz(-1.4186335) q[1];
rz(-pi) q[2];
x q[2];
rz(0.078048869) q[3];
sx q[3];
rz(-2.4025318) q[3];
sx q[3];
rz(-0.8821677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.687872) q[2];
sx q[2];
rz(-2.5747955) q[2];
sx q[2];
rz(2.4570214) q[2];
rz(1.941393) q[3];
sx q[3];
rz(-1.129351) q[3];
sx q[3];
rz(2.9620192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0612563) q[0];
sx q[0];
rz(-0.48114023) q[0];
sx q[0];
rz(0.0048333724) q[0];
rz(1.5176516) q[1];
sx q[1];
rz(-2.3976517) q[1];
sx q[1];
rz(0.25371107) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1954945) q[0];
sx q[0];
rz(-1.4854447) q[0];
sx q[0];
rz(-1.781989) q[0];
rz(-pi) q[1];
rz(0.51400916) q[2];
sx q[2];
rz(-2.3367662) q[2];
sx q[2];
rz(0.45586205) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.0066792329) q[1];
sx q[1];
rz(-0.90613922) q[1];
sx q[1];
rz(0.24941872) q[1];
rz(-pi) q[2];
x q[2];
rz(0.55048914) q[3];
sx q[3];
rz(-1.7409678) q[3];
sx q[3];
rz(-2.6941534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.97810513) q[2];
sx q[2];
rz(-1.9065964) q[2];
sx q[2];
rz(1.9592436) q[2];
rz(2.4649418) q[3];
sx q[3];
rz(-0.38656056) q[3];
sx q[3];
rz(-2.1277908) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61426281) q[0];
sx q[0];
rz(-2.3557721) q[0];
sx q[0];
rz(-2.8489805) q[0];
rz(-1.0489427) q[1];
sx q[1];
rz(-1.7221778) q[1];
sx q[1];
rz(-2.4377951) q[1];
rz(-0.47164698) q[2];
sx q[2];
rz(-1.8744962) q[2];
sx q[2];
rz(2.9203109) q[2];
rz(-1.9378035) q[3];
sx q[3];
rz(-1.4784151) q[3];
sx q[3];
rz(-2.9611931) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
