OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3192531) q[0];
sx q[0];
rz(-2.1847794) q[0];
sx q[0];
rz(-1.4332888) q[0];
rz(2.9046471) q[1];
sx q[1];
rz(-1.2304996) q[1];
sx q[1];
rz(1.1448316) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.618453) q[0];
sx q[0];
rz(-1.8123527) q[0];
sx q[0];
rz(-0.068058204) q[0];
rz(-pi) q[1];
rz(-0.28540622) q[2];
sx q[2];
rz(-1.1661582) q[2];
sx q[2];
rz(1.0813431) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.634234) q[1];
sx q[1];
rz(-0.45295742) q[1];
sx q[1];
rz(-1.5276315) q[1];
x q[2];
rz(-2.5979795) q[3];
sx q[3];
rz(-1.8279148) q[3];
sx q[3];
rz(1.4527814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.11674374) q[2];
sx q[2];
rz(-1.3329093) q[2];
sx q[2];
rz(-1.9072745) q[2];
rz(2.0862789) q[3];
sx q[3];
rz(-1.0588667) q[3];
sx q[3];
rz(1.9887259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.9556483) q[0];
sx q[0];
rz(-1.9444822) q[0];
sx q[0];
rz(0.82988513) q[0];
rz(2.8886967) q[1];
sx q[1];
rz(-0.77650944) q[1];
sx q[1];
rz(0.84709644) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0711489) q[0];
sx q[0];
rz(-1.1304454) q[0];
sx q[0];
rz(-1.1308934) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8663835) q[2];
sx q[2];
rz(-2.327773) q[2];
sx q[2];
rz(0.57397599) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9287195) q[1];
sx q[1];
rz(-1.6391014) q[1];
sx q[1];
rz(2.5961155) q[1];
rz(-0.3470207) q[3];
sx q[3];
rz(-1.1565398) q[3];
sx q[3];
rz(1.8300717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0236686) q[2];
sx q[2];
rz(-2.7682722) q[2];
sx q[2];
rz(0.71933293) q[2];
rz(1.2891278) q[3];
sx q[3];
rz(-0.23580655) q[3];
sx q[3];
rz(-1.4891362) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2704724) q[0];
sx q[0];
rz(-1.8914762) q[0];
sx q[0];
rz(-2.5701994) q[0];
rz(2.1748623) q[1];
sx q[1];
rz(-1.2582018) q[1];
sx q[1];
rz(-2.6142696) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8910599) q[0];
sx q[0];
rz(-2.8322304) q[0];
sx q[0];
rz(3.1191349) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8965365) q[2];
sx q[2];
rz(-2.2176761) q[2];
sx q[2];
rz(1.5145472) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.61422435) q[1];
sx q[1];
rz(-2.3389611) q[1];
sx q[1];
rz(2.7781092) q[1];
rz(-pi) q[2];
rz(2.4089291) q[3];
sx q[3];
rz(-1.3244197) q[3];
sx q[3];
rz(0.25783595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.30806914) q[2];
sx q[2];
rz(-1.6922733) q[2];
sx q[2];
rz(-2.4148338) q[2];
rz(-2.1440992) q[3];
sx q[3];
rz(-0.62524978) q[3];
sx q[3];
rz(-1.6231026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24546394) q[0];
sx q[0];
rz(-1.6587057) q[0];
sx q[0];
rz(2.1380651) q[0];
rz(0.040680496) q[1];
sx q[1];
rz(-2.0122806) q[1];
sx q[1];
rz(-0.8262659) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59506455) q[0];
sx q[0];
rz(-1.0142769) q[0];
sx q[0];
rz(0.80842774) q[0];
rz(-pi) q[1];
rz(-3.0777061) q[2];
sx q[2];
rz(-0.70110828) q[2];
sx q[2];
rz(-2.137616) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.77376765) q[1];
sx q[1];
rz(-2.7959679) q[1];
sx q[1];
rz(1.7961545) q[1];
rz(-0.63185933) q[3];
sx q[3];
rz(-1.251289) q[3];
sx q[3];
rz(-2.7416122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.56759175) q[2];
sx q[2];
rz(-1.2446128) q[2];
sx q[2];
rz(-2.2272002) q[2];
rz(-0.10522035) q[3];
sx q[3];
rz(-1.6243694) q[3];
sx q[3];
rz(2.5523394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2545664) q[0];
sx q[0];
rz(-1.6396602) q[0];
sx q[0];
rz(2.0078833) q[0];
rz(-3.1359613) q[1];
sx q[1];
rz(-1.7040323) q[1];
sx q[1];
rz(0.61757913) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.782398) q[0];
sx q[0];
rz(-0.054464666) q[0];
sx q[0];
rz(-2.1468303) q[0];
rz(-0.56474833) q[2];
sx q[2];
rz(-2.4167633) q[2];
sx q[2];
rz(0.85541475) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2684106) q[1];
sx q[1];
rz(-2.1408484) q[1];
sx q[1];
rz(2.8084055) q[1];
rz(-pi) q[2];
rz(-0.17554749) q[3];
sx q[3];
rz(-1.5581589) q[3];
sx q[3];
rz(-0.73568425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3474943) q[2];
sx q[2];
rz(-1.2604875) q[2];
sx q[2];
rz(2.9079672) q[2];
rz(0.99308333) q[3];
sx q[3];
rz(-2.1916094) q[3];
sx q[3];
rz(-0.63794199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1722906) q[0];
sx q[0];
rz(-3.0811221) q[0];
sx q[0];
rz(-3.0517975) q[0];
rz(0.88109294) q[1];
sx q[1];
rz(-0.52905622) q[1];
sx q[1];
rz(-2.9615013) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0952547) q[0];
sx q[0];
rz(-1.006608) q[0];
sx q[0];
rz(2.6020223) q[0];
rz(-1.300235) q[2];
sx q[2];
rz(-0.62108835) q[2];
sx q[2];
rz(-1.5121216) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5345726) q[1];
sx q[1];
rz(-1.8587451) q[1];
sx q[1];
rz(1.401591) q[1];
x q[2];
rz(-3.1230934) q[3];
sx q[3];
rz(-2.0634994) q[3];
sx q[3];
rz(0.25445709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.2580516) q[2];
sx q[2];
rz(-1.5641314) q[2];
sx q[2];
rz(1.1759261) q[2];
rz(-0.073143395) q[3];
sx q[3];
rz(-0.72435838) q[3];
sx q[3];
rz(0.49889645) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95773762) q[0];
sx q[0];
rz(-0.65559214) q[0];
sx q[0];
rz(-1.4181597) q[0];
rz(1.9765967) q[1];
sx q[1];
rz(-0.55924758) q[1];
sx q[1];
rz(3.1013536) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9366074) q[0];
sx q[0];
rz(-1.4892704) q[0];
sx q[0];
rz(1.8100912) q[0];
rz(-1.4637429) q[2];
sx q[2];
rz(-2.7333626) q[2];
sx q[2];
rz(-1.6183311) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.096958728) q[1];
sx q[1];
rz(-2.5068336) q[1];
sx q[1];
rz(1.8083014) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2651029) q[3];
sx q[3];
rz(-1.2672658) q[3];
sx q[3];
rz(-0.85206735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0030901) q[2];
sx q[2];
rz(-0.78531992) q[2];
sx q[2];
rz(-2.9677532) q[2];
rz(1.396817) q[3];
sx q[3];
rz(-1.4649748) q[3];
sx q[3];
rz(2.2824536) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1180856) q[0];
sx q[0];
rz(-0.21723391) q[0];
sx q[0];
rz(1.7370976) q[0];
rz(-1.2069758) q[1];
sx q[1];
rz(-1.6563621) q[1];
sx q[1];
rz(1.5054024) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0818427) q[0];
sx q[0];
rz(-1.2568226) q[0];
sx q[0];
rz(-0.44072515) q[0];
rz(-2.5152399) q[2];
sx q[2];
rz(-0.92712958) q[2];
sx q[2];
rz(-2.6477637) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9291032) q[1];
sx q[1];
rz(-0.87506064) q[1];
sx q[1];
rz(1.6131669) q[1];
x q[2];
rz(-1.1106311) q[3];
sx q[3];
rz(-0.65350973) q[3];
sx q[3];
rz(1.489952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76190844) q[2];
sx q[2];
rz(-2.6613993) q[2];
sx q[2];
rz(-0.15052477) q[2];
rz(1.5395509) q[3];
sx q[3];
rz(-1.9463041) q[3];
sx q[3];
rz(0.62301821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6983011) q[0];
sx q[0];
rz(-2.0069831) q[0];
sx q[0];
rz(-0.64569965) q[0];
rz(0.49939108) q[1];
sx q[1];
rz(-1.4285409) q[1];
sx q[1];
rz(1.2649149) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9547573) q[0];
sx q[0];
rz(-1.1823913) q[0];
sx q[0];
rz(-2.3507621) q[0];
rz(-pi) q[1];
rz(0.55391295) q[2];
sx q[2];
rz(-1.1914807) q[2];
sx q[2];
rz(1.5253138) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.39667323) q[1];
sx q[1];
rz(-1.6066237) q[1];
sx q[1];
rz(0.14573914) q[1];
rz(1.3690788) q[3];
sx q[3];
rz(-1.9864051) q[3];
sx q[3];
rz(-2.6430074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.50679961) q[2];
sx q[2];
rz(-2.1719077) q[2];
sx q[2];
rz(0.52465049) q[2];
rz(0.55650416) q[3];
sx q[3];
rz(-1.6289214) q[3];
sx q[3];
rz(2.7867253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5234579) q[0];
sx q[0];
rz(-2.3045461) q[0];
sx q[0];
rz(0.24542228) q[0];
rz(-0.55631176) q[1];
sx q[1];
rz(-1.6106771) q[1];
sx q[1];
rz(1.6773178) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3905555) q[0];
sx q[0];
rz(-1.7883736) q[0];
sx q[0];
rz(-1.6839954) q[0];
x q[1];
rz(-1.9015058) q[2];
sx q[2];
rz(-0.61243528) q[2];
sx q[2];
rz(1.2182747) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.97800868) q[1];
sx q[1];
rz(-2.2469236) q[1];
sx q[1];
rz(2.7917557) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.085300307) q[3];
sx q[3];
rz(-1.7461667) q[3];
sx q[3];
rz(2.8773957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2422553) q[2];
sx q[2];
rz(-1.238404) q[2];
sx q[2];
rz(2.396092) q[2];
rz(-1.7177104) q[3];
sx q[3];
rz(-0.29885492) q[3];
sx q[3];
rz(-1.1283114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2765008) q[0];
sx q[0];
rz(-1.4958953) q[0];
sx q[0];
rz(-2.4687742) q[0];
rz(-1.8854234) q[1];
sx q[1];
rz(-0.80614631) q[1];
sx q[1];
rz(2.0731906) q[1];
rz(-2.4297498) q[2];
sx q[2];
rz(-1.9528452) q[2];
sx q[2];
rz(-0.14609329) q[2];
rz(0.57337702) q[3];
sx q[3];
rz(-1.8918512) q[3];
sx q[3];
rz(0.20791114) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
