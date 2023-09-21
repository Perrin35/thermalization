OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.723974347114563) q[0];
sx q[0];
rz(1.48998776276643) q[0];
sx q[0];
rz(8.49433538912936) q[0];
rz(0.62970495223999) q[1];
sx q[1];
rz(4.27607444127137) q[1];
sx q[1];
rz(8.31743357180759) q[1];
cx q[1],q[0];
rz(-0.212702766060829) q[0];
sx q[0];
rz(5.03758564789826) q[0];
sx q[0];
rz(13.4297241926114) q[0];
rz(0.532221138477325) q[2];
sx q[2];
rz(3.93630436261232) q[2];
sx q[2];
rz(9.88671711682483) q[2];
cx q[2],q[1];
rz(2.10792303085327) q[1];
sx q[1];
rz(2.84069970448548) q[1];
sx q[1];
rz(6.66987965106174) q[1];
rz(-0.635251522064209) q[3];
sx q[3];
rz(4.9971218426996) q[3];
sx q[3];
rz(10.3487162351529) q[3];
cx q[3],q[2];
rz(0.0636231079697609) q[2];
sx q[2];
rz(3.87063232262666) q[2];
sx q[2];
rz(8.09675059317752) q[2];
rz(-0.320874571800232) q[3];
sx q[3];
rz(4.12754801114137) q[3];
sx q[3];
rz(12.4343950509946) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.48224738240242) q[0];
sx q[0];
rz(6.17079916794831) q[0];
sx q[0];
rz(8.54415712355777) q[0];
rz(4.43567132949829) q[1];
sx q[1];
rz(6.70113650162751) q[1];
sx q[1];
rz(10.2420418023984) q[1];
cx q[1],q[0];
rz(-1.19992315769196) q[0];
sx q[0];
rz(3.51620000799234) q[0];
sx q[0];
rz(9.03075311183139) q[0];
rz(-4.93059968948364) q[2];
sx q[2];
rz(5.193015011149) q[2];
sx q[2];
rz(13.0461513757627) q[2];
cx q[2],q[1];
rz(-0.686899900436401) q[1];
sx q[1];
rz(5.341776521998) q[1];
sx q[1];
rz(15.3901409864347) q[1];
rz(3.49861025810242) q[3];
sx q[3];
rz(3.64178005059297) q[3];
sx q[3];
rz(6.86069295405551) q[3];
cx q[3],q[2];
rz(4.05943489074707) q[2];
sx q[2];
rz(5.6023327430063) q[2];
sx q[2];
rz(9.78880534171268) q[2];
rz(-4.12797260284424) q[3];
sx q[3];
rz(4.86634448369081) q[3];
sx q[3];
rz(7.74782905577823) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.77469754219055) q[0];
sx q[0];
rz(5.45083990891511) q[0];
sx q[0];
rz(11.6000554323117) q[0];
rz(-0.192933827638626) q[1];
sx q[1];
rz(5.19455185730989) q[1];
sx q[1];
rz(14.0134415388028) q[1];
cx q[1],q[0];
rz(-0.236117854714394) q[0];
sx q[0];
rz(3.35731543798978) q[0];
sx q[0];
rz(11.0777879714887) q[0];
rz(2.85084939002991) q[2];
sx q[2];
rz(4.63313332398469) q[2];
sx q[2];
rz(8.79850170611545) q[2];
cx q[2],q[1];
rz(3.3822979927063) q[1];
sx q[1];
rz(5.08656147320802) q[1];
sx q[1];
rz(11.64527962207) q[1];
rz(-0.514794766902924) q[3];
sx q[3];
rz(0.47009530861909) q[3];
sx q[3];
rz(11.0445415735166) q[3];
cx q[3],q[2];
rz(0.438592791557312) q[2];
sx q[2];
rz(2.6575227697664) q[2];
sx q[2];
rz(10.5299995899121) q[2];
rz(2.39537191390991) q[3];
sx q[3];
rz(4.63096288044984) q[3];
sx q[3];
rz(8.44905022381946) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.00639092922211) q[0];
sx q[0];
rz(0.524403007822581) q[0];
sx q[0];
rz(10.890723681442) q[0];
rz(-0.284945964813232) q[1];
sx q[1];
rz(5.21286383469636) q[1];
sx q[1];
rz(9.03579143284961) q[1];
cx q[1],q[0];
rz(1.32724142074585) q[0];
sx q[0];
rz(4.4301359971338) q[0];
sx q[0];
rz(9.18591727911636) q[0];
rz(0.35819274187088) q[2];
sx q[2];
rz(3.65573296149308) q[2];
sx q[2];
rz(8.91124764680072) q[2];
cx q[2],q[1];
rz(3.04918384552002) q[1];
sx q[1];
rz(5.25748697121675) q[1];
sx q[1];
rz(8.92428359984561) q[1];
rz(1.34974706172943) q[3];
sx q[3];
rz(3.90182951291139) q[3];
sx q[3];
rz(7.86190674304172) q[3];
cx q[3],q[2];
rz(-1.41855621337891) q[2];
sx q[2];
rz(4.30523875554139) q[2];
sx q[2];
rz(13.0230767488401) q[2];
rz(-1.51519536972046) q[3];
sx q[3];
rz(4.09066608746583) q[3];
sx q[3];
rz(12.2840277910154) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.919637680053711) q[0];
sx q[0];
rz(4.2813448031717) q[0];
sx q[0];
rz(9.78480253218814) q[0];
rz(2.49417638778687) q[1];
sx q[1];
rz(4.67081657250459) q[1];
sx q[1];
rz(8.93027037977382) q[1];
cx q[1],q[0];
rz(1.76171028614044) q[0];
sx q[0];
rz(3.44870618184144) q[0];
sx q[0];
rz(10.5733743667524) q[0];
rz(1.61022293567657) q[2];
sx q[2];
rz(4.84774294694001) q[2];
sx q[2];
rz(8.89212754964038) q[2];
cx q[2],q[1];
rz(-4.42118358612061) q[1];
sx q[1];
rz(1.51085666020448) q[1];
sx q[1];
rz(10.3651407718579) q[1];
rz(3.11616063117981) q[3];
sx q[3];
rz(4.26841464837129) q[3];
sx q[3];
rz(8.62448481320545) q[3];
cx q[3],q[2];
rz(-2.43137240409851) q[2];
sx q[2];
rz(3.7974757870012) q[2];
sx q[2];
rz(5.98467824458286) q[2];
rz(-0.31202894449234) q[3];
sx q[3];
rz(8.09407153924043) q[3];
sx q[3];
rz(8.7862792968671) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.145377799868584) q[0];
sx q[0];
rz(5.56010428269441) q[0];
sx q[0];
rz(9.22887082993194) q[0];
rz(3.16267514228821) q[1];
sx q[1];
rz(7.68177285988862) q[1];
sx q[1];
rz(10.6599769353788) q[1];
cx q[1],q[0];
rz(0.933419525623322) q[0];
sx q[0];
rz(6.74467030365998) q[0];
sx q[0];
rz(13.3911149263303) q[0];
rz(5.94450187683105) q[2];
sx q[2];
rz(2.24656960566575) q[2];
sx q[2];
rz(8.37287793158695) q[2];
cx q[2],q[1];
rz(-2.13933086395264) q[1];
sx q[1];
rz(1.98632231553132) q[1];
sx q[1];
rz(11.6321182012479) q[1];
rz(-1.33448684215546) q[3];
sx q[3];
rz(1.65225389798219) q[3];
sx q[3];
rz(12.8443357705991) q[3];
cx q[3],q[2];
rz(3.63491773605347) q[2];
sx q[2];
rz(5.35705605347688) q[2];
sx q[2];
rz(8.99308524130985) q[2];
rz(4.87068700790405) q[3];
sx q[3];
rz(3.8639641722017) q[3];
sx q[3];
rz(9.28878154455825) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.391111850738525) q[0];
sx q[0];
rz(4.31040480931336) q[0];
sx q[0];
rz(9.53689434229537) q[0];
rz(-0.215133666992188) q[1];
sx q[1];
rz(4.72261038621003) q[1];
sx q[1];
rz(8.31137189864322) q[1];
cx q[1],q[0];
rz(2.45569014549255) q[0];
sx q[0];
rz(3.61262795527513) q[0];
sx q[0];
rz(7.25334546565219) q[0];
rz(1.25540113449097) q[2];
sx q[2];
rz(1.8176259120279) q[2];
sx q[2];
rz(8.46722528933688) q[2];
cx q[2],q[1];
rz(-1.09431624412537) q[1];
sx q[1];
rz(5.39183321793611) q[1];
sx q[1];
rz(15.3132824659269) q[1];
rz(1.35735356807709) q[3];
sx q[3];
rz(5.15078988869721) q[3];
sx q[3];
rz(11.4060896396558) q[3];
cx q[3],q[2];
rz(0.00117428798694164) q[2];
sx q[2];
rz(5.55053821404512) q[2];
sx q[2];
rz(9.29874911009475) q[2];
rz(1.04729390144348) q[3];
sx q[3];
rz(4.96244874795014) q[3];
sx q[3];
rz(11.8581690549771) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.660131812095642) q[0];
sx q[0];
rz(5.52449837525422) q[0];
sx q[0];
rz(7.96460196971103) q[0];
rz(-1.24496448040009) q[1];
sx q[1];
rz(4.23591283162171) q[1];
sx q[1];
rz(10.633680677406) q[1];
cx q[1],q[0];
rz(-0.932232439517975) q[0];
sx q[0];
rz(2.86738422711427) q[0];
sx q[0];
rz(9.12822116016551) q[0];
rz(2.40215468406677) q[2];
sx q[2];
rz(5.23242560227449) q[2];
sx q[2];
rz(8.06405827998325) q[2];
cx q[2],q[1];
rz(-4.47871828079224) q[1];
sx q[1];
rz(3.99941650231416) q[1];
sx q[1];
rz(8.98253399728938) q[1];
rz(1.79024791717529) q[3];
sx q[3];
rz(2.19005742867524) q[3];
sx q[3];
rz(7.37962172030612) q[3];
cx q[3],q[2];
rz(-5.9047417640686) q[2];
sx q[2];
rz(1.90387895901734) q[2];
sx q[2];
rz(14.3883633375089) q[2];
rz(0.59213262796402) q[3];
sx q[3];
rz(1.72546461422975) q[3];
sx q[3];
rz(9.38963625802799) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.232934400439262) q[0];
sx q[0];
rz(2.61805513699586) q[0];
sx q[0];
rz(7.64428780078098) q[0];
rz(-1.21104419231415) q[1];
sx q[1];
rz(5.37855306466157) q[1];
sx q[1];
rz(9.81646517514392) q[1];
cx q[1],q[0];
rz(2.21013331413269) q[0];
sx q[0];
rz(5.14691582520539) q[0];
sx q[0];
rz(6.27673313616916) q[0];
rz(-4.26007080078125) q[2];
sx q[2];
rz(1.75121334393556) q[2];
sx q[2];
rz(9.84663850664302) q[2];
cx q[2],q[1];
rz(1.0810250043869) q[1];
sx q[1];
rz(4.00681927998597) q[1];
sx q[1];
rz(7.41482493876621) q[1];
rz(3.2146909236908) q[3];
sx q[3];
rz(4.18785575230653) q[3];
sx q[3];
rz(11.6526989698331) q[3];
cx q[3],q[2];
rz(0.520355105400085) q[2];
sx q[2];
rz(1.34277990658815) q[2];
sx q[2];
rz(11.2295924186628) q[2];
rz(0.761986017227173) q[3];
sx q[3];
rz(2.82189941604669) q[3];
sx q[3];
rz(8.62440154551669) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(15/(14*pi)) q[0];
sx q[0];
rz(2.83881920774514) q[0];
sx q[0];
rz(8.85387853383228) q[0];
rz(-1.42924976348877) q[1];
sx q[1];
rz(5.21928301652009) q[1];
sx q[1];
rz(9.58671874403163) q[1];
cx q[1],q[0];
rz(0.256259441375732) q[0];
sx q[0];
rz(2.46180394490296) q[0];
sx q[0];
rz(9.89321819543048) q[0];
rz(-0.0477242320775986) q[2];
sx q[2];
rz(2.95279103715951) q[2];
sx q[2];
rz(10.827086544029) q[2];
cx q[2],q[1];
rz(0.754412770271301) q[1];
sx q[1];
rz(4.54091611702973) q[1];
sx q[1];
rz(10.865037059776) q[1];
rz(-0.0333719290792942) q[3];
sx q[3];
rz(4.18951609929139) q[3];
sx q[3];
rz(12.4996402025144) q[3];
cx q[3],q[2];
rz(-4.0990138053894) q[2];
sx q[2];
rz(4.24676838715608) q[2];
sx q[2];
rz(11.1381457805555) q[2];
rz(4.34102201461792) q[3];
sx q[3];
rz(4.2349410374933) q[3];
sx q[3];
rz(10.795412516586) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.93903076648712) q[0];
sx q[0];
rz(1.99555972416932) q[0];
sx q[0];
rz(13.9166788816373) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(0.908441722393036) q[1];
sx q[1];
rz(3.97938135464723) q[1];
sx q[1];
rz(8.5176692366521) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(2.25605392456055) q[2];
sx q[2];
rz(2.93209403951699) q[2];
sx q[2];
rz(10.202290570728) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-0.343266844749451) q[3];
sx q[3];
rz(1.97231367428834) q[3];
sx q[3];
rz(9.20666669904395) q[3];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];