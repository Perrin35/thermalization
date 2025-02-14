OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.51337564) q[0];
sx q[0];
rz(-1.6802508) q[0];
sx q[0];
rz(-0.914855) q[0];
rz(2.272361) q[1];
sx q[1];
rz(2.4183122) q[1];
sx q[1];
rz(8.0719168) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.992743) q[0];
sx q[0];
rz(-2.2701716) q[0];
sx q[0];
rz(-2.4464038) q[0];
x q[1];
rz(1.9454141) q[2];
sx q[2];
rz(-2.2785419) q[2];
sx q[2];
rz(1.5951235) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.1654403) q[1];
sx q[1];
rz(-2.4399477) q[1];
sx q[1];
rz(-1.2734423) q[1];
rz(-pi) q[2];
rz(2.1800024) q[3];
sx q[3];
rz(-1.7477027) q[3];
sx q[3];
rz(-0.78395236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4065518) q[2];
sx q[2];
rz(-2.328379) q[2];
sx q[2];
rz(-1.0068033) q[2];
rz(-1.6793647) q[3];
sx q[3];
rz(-1.6961325) q[3];
sx q[3];
rz(1.6026976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8524858) q[0];
sx q[0];
rz(-0.91133457) q[0];
sx q[0];
rz(0.12595969) q[0];
rz(2.0055298) q[1];
sx q[1];
rz(-1.2964396) q[1];
sx q[1];
rz(1.0596641) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9618767) q[0];
sx q[0];
rz(-1.5715412) q[0];
sx q[0];
rz(-0.0019093328) q[0];
rz(-pi) q[1];
rz(0.95086348) q[2];
sx q[2];
rz(-2.1427832) q[2];
sx q[2];
rz(2.7754155) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.48908614) q[1];
sx q[1];
rz(-0.87372696) q[1];
sx q[1];
rz(2.1537495) q[1];
rz(-pi) q[2];
rz(2.7823388) q[3];
sx q[3];
rz(-1.5240961) q[3];
sx q[3];
rz(0.31622313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.58553592) q[2];
sx q[2];
rz(-1.6126596) q[2];
sx q[2];
rz(2.5118828) q[2];
rz(1.9017879) q[3];
sx q[3];
rz(-2.3972378) q[3];
sx q[3];
rz(-2.0201717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.020551) q[0];
sx q[0];
rz(-2.910055) q[0];
sx q[0];
rz(1.9112021) q[0];
rz(-2.4379099) q[1];
sx q[1];
rz(-0.92862248) q[1];
sx q[1];
rz(-1.5603265) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11008292) q[0];
sx q[0];
rz(-0.39387273) q[0];
sx q[0];
rz(-1.5316233) q[0];
x q[1];
rz(2.4034924) q[2];
sx q[2];
rz(-1.4261275) q[2];
sx q[2];
rz(1.5371145) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3250503) q[1];
sx q[1];
rz(-0.61707622) q[1];
sx q[1];
rz(-2.3338334) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37902351) q[3];
sx q[3];
rz(-0.82250094) q[3];
sx q[3];
rz(-1.8387972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8065717) q[2];
sx q[2];
rz(-2.0609914) q[2];
sx q[2];
rz(-0.054277167) q[2];
rz(2.055376) q[3];
sx q[3];
rz(-0.79038668) q[3];
sx q[3];
rz(-2.6495892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7739968) q[0];
sx q[0];
rz(-3.1258899) q[0];
sx q[0];
rz(-2.1570461) q[0];
rz(-1.7355851) q[1];
sx q[1];
rz(-1.6008629) q[1];
sx q[1];
rz(-0.23449177) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2024883) q[0];
sx q[0];
rz(-1.2959769) q[0];
sx q[0];
rz(-2.2780134) q[0];
x q[1];
rz(-2.2735944) q[2];
sx q[2];
rz(-0.36938399) q[2];
sx q[2];
rz(-0.59472769) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.76012661) q[1];
sx q[1];
rz(-1.0262118) q[1];
sx q[1];
rz(-1.8167956) q[1];
rz(-pi) q[2];
rz(2.1470039) q[3];
sx q[3];
rz(-1.5458917) q[3];
sx q[3];
rz(1.2664317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8977114) q[2];
sx q[2];
rz(-2.477024) q[2];
sx q[2];
rz(0.74771869) q[2];
rz(-2.054935) q[3];
sx q[3];
rz(-2.1472411) q[3];
sx q[3];
rz(-2.7896816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20027593) q[0];
sx q[0];
rz(-0.87012297) q[0];
sx q[0];
rz(1.548832) q[0];
rz(3.0039655) q[1];
sx q[1];
rz(-2.177114) q[1];
sx q[1];
rz(-2.4639938) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0726822) q[0];
sx q[0];
rz(-3.125801) q[0];
sx q[0];
rz(0.48245211) q[0];
rz(-pi) q[1];
rz(0.65998544) q[2];
sx q[2];
rz(-2.8980245) q[2];
sx q[2];
rz(3.0483831) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2210834) q[1];
sx q[1];
rz(-1.5766605) q[1];
sx q[1];
rz(1.7076034) q[1];
x q[2];
rz(-1.3278392) q[3];
sx q[3];
rz(-2.2928709) q[3];
sx q[3];
rz(2.9606282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4981093) q[2];
sx q[2];
rz(-2.6124239) q[2];
sx q[2];
rz(-0.37437487) q[2];
rz(-2.3937285) q[3];
sx q[3];
rz(-1.6467843) q[3];
sx q[3];
rz(-1.098746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13136524) q[0];
sx q[0];
rz(-1.7759198) q[0];
sx q[0];
rz(2.9170872) q[0];
rz(2.7911216) q[1];
sx q[1];
rz(-1.1843362) q[1];
sx q[1];
rz(-1.1856461) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6939589) q[0];
sx q[0];
rz(-0.15800755) q[0];
sx q[0];
rz(-1.9290646) q[0];
rz(-pi) q[1];
rz(-2.9509928) q[2];
sx q[2];
rz(-2.6707044) q[2];
sx q[2];
rz(-0.48656551) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.396186) q[1];
sx q[1];
rz(-0.90382677) q[1];
sx q[1];
rz(1.0140258) q[1];
rz(-2.3697183) q[3];
sx q[3];
rz(-2.3079685) q[3];
sx q[3];
rz(-1.1590568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.4993569) q[2];
sx q[2];
rz(-2.2580052) q[2];
sx q[2];
rz(2.9767766) q[2];
rz(2.3757101) q[3];
sx q[3];
rz(-0.82847786) q[3];
sx q[3];
rz(-0.59144056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0931452) q[0];
sx q[0];
rz(-3.1146545) q[0];
sx q[0];
rz(-0.41937605) q[0];
rz(-1.5834437) q[1];
sx q[1];
rz(-2.7132468) q[1];
sx q[1];
rz(1.3126119) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87154138) q[0];
sx q[0];
rz(-1.4082552) q[0];
sx q[0];
rz(0.97515653) q[0];
x q[1];
rz(0.68676853) q[2];
sx q[2];
rz(-2.1109952) q[2];
sx q[2];
rz(1.6322002) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2129107) q[1];
sx q[1];
rz(-2.8281459) q[1];
sx q[1];
rz(-2.2261966) q[1];
rz(-pi) q[2];
rz(0.44371554) q[3];
sx q[3];
rz(-2.8464918) q[3];
sx q[3];
rz(2.091311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7340362) q[2];
sx q[2];
rz(-2.4702256) q[2];
sx q[2];
rz(2.7835795) q[2];
rz(0.19896209) q[3];
sx q[3];
rz(-0.81148654) q[3];
sx q[3];
rz(3.0437886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1426549) q[0];
sx q[0];
rz(-2.1079347) q[0];
sx q[0];
rz(2.348483) q[0];
rz(-2.2456887) q[1];
sx q[1];
rz(-1.221012) q[1];
sx q[1];
rz(2.6633247) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9288328) q[0];
sx q[0];
rz(-0.46657714) q[0];
sx q[0];
rz(-0.27884941) q[0];
x q[1];
rz(-2.7061126) q[2];
sx q[2];
rz(-1.9009095) q[2];
sx q[2];
rz(-0.94089905) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7859808) q[1];
sx q[1];
rz(-2.0352239) q[1];
sx q[1];
rz(2.8465438) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0317814) q[3];
sx q[3];
rz(-2.0834907) q[3];
sx q[3];
rz(-1.5231665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6200977) q[2];
sx q[2];
rz(-1.564881) q[2];
sx q[2];
rz(3.133797) q[2];
rz(-1.9589641) q[3];
sx q[3];
rz(-1.3820442) q[3];
sx q[3];
rz(1.8462605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.1496534) q[0];
sx q[0];
rz(-2.6872334) q[0];
sx q[0];
rz(0.40865189) q[0];
rz(-2.0463792) q[1];
sx q[1];
rz(-1.4827671) q[1];
sx q[1];
rz(2.2832787) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3811144) q[0];
sx q[0];
rz(-1.4234626) q[0];
sx q[0];
rz(3.0744746) q[0];
x q[1];
rz(0.13225358) q[2];
sx q[2];
rz(-2.0200854) q[2];
sx q[2];
rz(1.9455036) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.62095642) q[1];
sx q[1];
rz(-0.22527105) q[1];
sx q[1];
rz(-1.0566637) q[1];
x q[2];
rz(-0.59767234) q[3];
sx q[3];
rz(-1.7558756) q[3];
sx q[3];
rz(-1.5633893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3506763) q[2];
sx q[2];
rz(-2.3886267) q[2];
sx q[2];
rz(0.6616627) q[2];
rz(2.3412797) q[3];
sx q[3];
rz(-1.1626264) q[3];
sx q[3];
rz(-0.77320981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3043587) q[0];
sx q[0];
rz(-0.14166129) q[0];
sx q[0];
rz(-0.71453553) q[0];
rz(-0.47572687) q[1];
sx q[1];
rz(-2.5560502) q[1];
sx q[1];
rz(0.48019662) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59360945) q[0];
sx q[0];
rz(-2.3385323) q[0];
sx q[0];
rz(-0.70269967) q[0];
x q[1];
rz(1.9942029) q[2];
sx q[2];
rz(-2.4443279) q[2];
sx q[2];
rz(0.84691511) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3298397) q[1];
sx q[1];
rz(-1.8564714) q[1];
sx q[1];
rz(-2.4655254) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1350182) q[3];
sx q[3];
rz(-1.8284855) q[3];
sx q[3];
rz(2.1439752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6314038) q[2];
sx q[2];
rz(-1.1978585) q[2];
sx q[2];
rz(-0.27604827) q[2];
rz(-2.7035642) q[3];
sx q[3];
rz(-1.6249526) q[3];
sx q[3];
rz(-0.62694222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.940687) q[0];
sx q[0];
rz(-1.154366) q[0];
sx q[0];
rz(0.05703297) q[0];
rz(3.1055462) q[1];
sx q[1];
rz(-1.6135975) q[1];
sx q[1];
rz(1.5560908) q[1];
rz(2.8677058) q[2];
sx q[2];
rz(-1.0928434) q[2];
sx q[2];
rz(2.6437987) q[2];
rz(1.7071758) q[3];
sx q[3];
rz(-1.3188667) q[3];
sx q[3];
rz(1.8018166) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
