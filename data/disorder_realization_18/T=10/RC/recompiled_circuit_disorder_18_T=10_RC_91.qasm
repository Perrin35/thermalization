OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.82233959) q[0];
sx q[0];
rz(5.3263721) q[0];
sx q[0];
rz(10.858067) q[0];
rz(2.9046471) q[1];
sx q[1];
rz(-1.2304996) q[1];
sx q[1];
rz(-1.9967611) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0776318) q[0];
sx q[0];
rz(-1.6368757) q[0];
sx q[0];
rz(-1.8128916) q[0];
x q[1];
rz(-0.28540622) q[2];
sx q[2];
rz(-1.9754344) q[2];
sx q[2];
rz(2.0602496) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1022542) q[1];
sx q[1];
rz(-1.551911) q[1];
sx q[1];
rz(2.0233872) q[1];
rz(1.8688698) q[3];
sx q[3];
rz(-1.0469336) q[3];
sx q[3];
rz(0.03447547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0248489) q[2];
sx q[2];
rz(-1.8086834) q[2];
sx q[2];
rz(-1.2343181) q[2];
rz(-1.0553137) q[3];
sx q[3];
rz(-2.0827259) q[3];
sx q[3];
rz(1.1528667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9556483) q[0];
sx q[0];
rz(-1.1971104) q[0];
sx q[0];
rz(-2.3117075) q[0];
rz(0.25289598) q[1];
sx q[1];
rz(-0.77650944) q[1];
sx q[1];
rz(-0.84709644) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0711489) q[0];
sx q[0];
rz(-2.0111472) q[0];
sx q[0];
rz(2.0106993) q[0];
rz(-pi) q[1];
rz(2.3624624) q[2];
sx q[2];
rz(-1.7841633) q[2];
sx q[2];
rz(-1.9386292) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.21287316) q[1];
sx q[1];
rz(-1.5024912) q[1];
sx q[1];
rz(2.5961155) q[1];
rz(-2.0081677) q[3];
sx q[3];
rz(-1.2542033) q[3];
sx q[3];
rz(-0.40382995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0236686) q[2];
sx q[2];
rz(-0.37332049) q[2];
sx q[2];
rz(2.4222597) q[2];
rz(-1.2891278) q[3];
sx q[3];
rz(-0.23580655) q[3];
sx q[3];
rz(1.4891362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8711202) q[0];
sx q[0];
rz(-1.8914762) q[0];
sx q[0];
rz(0.57139325) q[0];
rz(-2.1748623) q[1];
sx q[1];
rz(-1.2582018) q[1];
sx q[1];
rz(2.6142696) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2505328) q[0];
sx q[0];
rz(-0.30936229) q[0];
sx q[0];
rz(-0.022457794) q[0];
rz(-pi) q[1];
rz(-2.4685523) q[2];
sx q[2];
rz(-1.8290142) q[2];
sx q[2];
rz(2.8845127) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5273683) q[1];
sx q[1];
rz(-0.80263153) q[1];
sx q[1];
rz(-2.7781092) q[1];
x q[2];
rz(2.4089291) q[3];
sx q[3];
rz(-1.3244197) q[3];
sx q[3];
rz(-2.8837567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.30806914) q[2];
sx q[2];
rz(-1.4493194) q[2];
sx q[2];
rz(-2.4148338) q[2];
rz(0.99749342) q[3];
sx q[3];
rz(-0.62524978) q[3];
sx q[3];
rz(1.5184901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24546394) q[0];
sx q[0];
rz(-1.4828869) q[0];
sx q[0];
rz(-2.1380651) q[0];
rz(-3.1009122) q[1];
sx q[1];
rz(-1.129312) q[1];
sx q[1];
rz(0.8262659) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6711133) q[0];
sx q[0];
rz(-0.90958909) q[0];
sx q[0];
rz(0.83755042) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0777061) q[2];
sx q[2];
rz(-2.4404844) q[2];
sx q[2];
rz(2.137616) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1288209) q[1];
sx q[1];
rz(-1.2342617) q[1];
sx q[1];
rz(3.0613042) q[1];
rz(-pi) q[2];
rz(-0.63185933) q[3];
sx q[3];
rz(-1.8903036) q[3];
sx q[3];
rz(-0.39998049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.56759175) q[2];
sx q[2];
rz(-1.8969798) q[2];
sx q[2];
rz(-0.91439247) q[2];
rz(0.10522035) q[3];
sx q[3];
rz(-1.5172232) q[3];
sx q[3];
rz(-0.58925327) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88702622) q[0];
sx q[0];
rz(-1.6396602) q[0];
sx q[0];
rz(1.1337093) q[0];
rz(-0.0056313593) q[1];
sx q[1];
rz(-1.7040323) q[1];
sx q[1];
rz(2.5240135) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3546346) q[0];
sx q[0];
rz(-1.5411396) q[0];
sx q[0];
rz(1.5251072) q[0];
x q[1];
rz(-2.0134301) q[2];
sx q[2];
rz(-0.97634041) q[2];
sx q[2];
rz(0.15304676) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2684106) q[1];
sx q[1];
rz(-2.1408484) q[1];
sx q[1];
rz(2.8084055) q[1];
rz(-pi) q[2];
rz(-1.5836309) q[3];
sx q[3];
rz(-1.395263) q[3];
sx q[3];
rz(-0.83735355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.79409838) q[2];
sx q[2];
rz(-1.2604875) q[2];
sx q[2];
rz(0.23362544) q[2];
rz(2.1485093) q[3];
sx q[3];
rz(-0.94998327) q[3];
sx q[3];
rz(2.5036507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1722906) q[0];
sx q[0];
rz(-3.0811221) q[0];
sx q[0];
rz(-3.0517975) q[0];
rz(-2.2604997) q[1];
sx q[1];
rz(-0.52905622) q[1];
sx q[1];
rz(0.18009137) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20443944) q[0];
sx q[0];
rz(-0.75980543) q[0];
sx q[0];
rz(2.252749) q[0];
rz(-pi) q[1];
rz(-2.9526268) q[2];
sx q[2];
rz(-2.1660888) q[2];
sx q[2];
rz(-1.8408066) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1485968) q[1];
sx q[1];
rz(-2.8088048) q[1];
sx q[1];
rz(-2.6246043) q[1];
rz(-pi) q[2];
rz(-1.0780219) q[3];
sx q[3];
rz(-1.587095) q[3];
sx q[3];
rz(1.3075882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.883541) q[2];
sx q[2];
rz(-1.5641314) q[2];
sx q[2];
rz(-1.9656666) q[2];
rz(0.073143395) q[3];
sx q[3];
rz(-0.72435838) q[3];
sx q[3];
rz(-0.49889645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.183855) q[0];
sx q[0];
rz(-2.4860005) q[0];
sx q[0];
rz(-1.7234329) q[0];
rz(1.9765967) q[1];
sx q[1];
rz(-2.5823451) q[1];
sx q[1];
rz(-3.1013536) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0979472) q[0];
sx q[0];
rz(-0.25254927) q[0];
sx q[0];
rz(-1.2388242) q[0];
rz(-pi) q[1];
rz(-1.6778498) q[2];
sx q[2];
rz(-0.40823001) q[2];
sx q[2];
rz(-1.6183311) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.28133) q[1];
sx q[1];
rz(-1.7107692) q[1];
sx q[1];
rz(-2.1920188) q[1];
rz(-pi) q[2];
rz(-0.76544806) q[3];
sx q[3];
rz(-0.42740373) q[3];
sx q[3];
rz(-3.1020853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13850257) q[2];
sx q[2];
rz(-0.78531992) q[2];
sx q[2];
rz(0.17383943) q[2];
rz(1.7447757) q[3];
sx q[3];
rz(-1.4649748) q[3];
sx q[3];
rz(-2.2824536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.023507) q[0];
sx q[0];
rz(-2.9243587) q[0];
sx q[0];
rz(-1.404495) q[0];
rz(1.9346168) q[1];
sx q[1];
rz(-1.6563621) q[1];
sx q[1];
rz(1.5054024) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0818427) q[0];
sx q[0];
rz(-1.2568226) q[0];
sx q[0];
rz(-0.44072515) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3178188) q[2];
sx q[2];
rz(-2.0588377) q[2];
sx q[2];
rz(1.6549695) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7561188) q[1];
sx q[1];
rz(-1.5382775) q[1];
sx q[1];
rz(-2.4454152) q[1];
x q[2];
rz(2.1721341) q[3];
sx q[3];
rz(-1.8441895) q[3];
sx q[3];
rz(-2.6858342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3796842) q[2];
sx q[2];
rz(-2.6613993) q[2];
sx q[2];
rz(-2.9910679) q[2];
rz(1.6020417) q[3];
sx q[3];
rz(-1.1952885) q[3];
sx q[3];
rz(0.62301821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6983011) q[0];
sx q[0];
rz(-1.1346096) q[0];
sx q[0];
rz(-2.495893) q[0];
rz(2.6422016) q[1];
sx q[1];
rz(-1.4285409) q[1];
sx q[1];
rz(1.8766778) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0255233) q[0];
sx q[0];
rz(-2.2795296) q[0];
sx q[0];
rz(-0.52225964) q[0];
x q[1];
rz(1.132498) q[2];
sx q[2];
rz(-2.0813137) q[2];
sx q[2];
rz(0.27062182) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9727271) q[1];
sx q[1];
rz(-1.4251514) q[1];
sx q[1];
rz(-1.6070073) q[1];
x q[2];
rz(-0.42322741) q[3];
sx q[3];
rz(-1.3864577) q[3];
sx q[3];
rz(1.9870027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.634793) q[2];
sx q[2];
rz(-0.96968499) q[2];
sx q[2];
rz(-0.52465049) q[2];
rz(0.55650416) q[3];
sx q[3];
rz(-1.5126712) q[3];
sx q[3];
rz(0.35486737) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5234579) q[0];
sx q[0];
rz(-0.83704656) q[0];
sx q[0];
rz(2.8961704) q[0];
rz(-2.5852809) q[1];
sx q[1];
rz(-1.6106771) q[1];
sx q[1];
rz(-1.6773178) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2357764) q[0];
sx q[0];
rz(-0.24484867) q[0];
sx q[0];
rz(2.66923) q[0];
rz(-2.9173031) q[2];
sx q[2];
rz(-2.1456492) q[2];
sx q[2];
rz(1.526051) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.36832419) q[1];
sx q[1];
rz(-1.8414294) q[1];
sx q[1];
rz(0.86398217) q[1];
x q[2];
rz(-1.1226095) q[3];
sx q[3];
rz(-2.9467694) q[3];
sx q[3];
rz(-0.71988718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8993373) q[2];
sx q[2];
rz(-1.9031886) q[2];
sx q[2];
rz(-0.74550068) q[2];
rz(-1.7177104) q[3];
sx q[3];
rz(-2.8427377) q[3];
sx q[3];
rz(-2.0132813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2765008) q[0];
sx q[0];
rz(-1.6456974) q[0];
sx q[0];
rz(0.67281848) q[0];
rz(1.2561692) q[1];
sx q[1];
rz(-0.80614631) q[1];
sx q[1];
rz(2.0731906) q[1];
rz(-0.55143572) q[2];
sx q[2];
rz(-0.79179344) q[2];
sx q[2];
rz(1.0168016) q[2];
rz(-1.1938548) q[3];
sx q[3];
rz(-1.0300763) q[3];
sx q[3];
rz(-1.1618617) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
