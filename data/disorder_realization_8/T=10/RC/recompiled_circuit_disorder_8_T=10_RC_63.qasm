OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8060057) q[0];
sx q[0];
rz(-0.94526362) q[0];
sx q[0];
rz(2.6160016) q[0];
rz(-2.8984012) q[1];
sx q[1];
rz(-1.2326198) q[1];
sx q[1];
rz(-0.90484172) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5072767) q[0];
sx q[0];
rz(-2.0024558) q[0];
sx q[0];
rz(-1.0213724) q[0];
rz(0.71360795) q[2];
sx q[2];
rz(-2.8461694) q[2];
sx q[2];
rz(-1.7507391) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.875784) q[1];
sx q[1];
rz(-1.335653) q[1];
sx q[1];
rz(-1.9782515) q[1];
rz(-pi) q[2];
rz(1.8815133) q[3];
sx q[3];
rz(-1.7508535) q[3];
sx q[3];
rz(0.77375274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9378172) q[2];
sx q[2];
rz(-1.7353461) q[2];
sx q[2];
rz(-3.0453483) q[2];
rz(-1.0359267) q[3];
sx q[3];
rz(-2.7544498) q[3];
sx q[3];
rz(2.9878785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7006943) q[0];
sx q[0];
rz(-0.39114025) q[0];
sx q[0];
rz(0.76517117) q[0];
rz(-1.8493429) q[1];
sx q[1];
rz(-2.6563829) q[1];
sx q[1];
rz(-0.66295019) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5040341) q[0];
sx q[0];
rz(-0.090888977) q[0];
sx q[0];
rz(0.81764098) q[0];
x q[1];
rz(2.2234369) q[2];
sx q[2];
rz(-1.895972) q[2];
sx q[2];
rz(-2.6174389) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4274802) q[1];
sx q[1];
rz(-2.6926827) q[1];
sx q[1];
rz(2.950579) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29942056) q[3];
sx q[3];
rz(-0.50588183) q[3];
sx q[3];
rz(-2.8477856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.092676) q[2];
sx q[2];
rz(-1.0485704) q[2];
sx q[2];
rz(0.32901397) q[2];
rz(-0.66550231) q[3];
sx q[3];
rz(-0.21829675) q[3];
sx q[3];
rz(-1.3177419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5927785) q[0];
sx q[0];
rz(-2.9187262) q[0];
sx q[0];
rz(2.9192525) q[0];
rz(1.0173343) q[1];
sx q[1];
rz(-0.72128123) q[1];
sx q[1];
rz(2.6229048) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1919353) q[0];
sx q[0];
rz(-1.2149095) q[0];
sx q[0];
rz(0.8215254) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2661744) q[2];
sx q[2];
rz(-1.4783579) q[2];
sx q[2];
rz(0.12649525) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.31836244) q[1];
sx q[1];
rz(-1.743411) q[1];
sx q[1];
rz(2.361972) q[1];
rz(-pi) q[2];
rz(-1.5699925) q[3];
sx q[3];
rz(-1.5161848) q[3];
sx q[3];
rz(-2.2664203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2805933) q[2];
sx q[2];
rz(-2.7797647) q[2];
sx q[2];
rz(-3.0730491) q[2];
rz(-2.5391501) q[3];
sx q[3];
rz(-0.76255637) q[3];
sx q[3];
rz(3.0025735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.65072) q[0];
sx q[0];
rz(-0.77712178) q[0];
sx q[0];
rz(-0.17424507) q[0];
rz(-0.53025591) q[1];
sx q[1];
rz(-1.4825876) q[1];
sx q[1];
rz(0.51309103) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1888694) q[0];
sx q[0];
rz(-2.1892622) q[0];
sx q[0];
rz(-1.6678715) q[0];
rz(-pi) q[1];
rz(-2.1999173) q[2];
sx q[2];
rz(-0.86949124) q[2];
sx q[2];
rz(2.6710682) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9070223) q[1];
sx q[1];
rz(-1.4325732) q[1];
sx q[1];
rz(2.6837818) q[1];
rz(-pi) q[2];
rz(1.3644049) q[3];
sx q[3];
rz(-1.6009814) q[3];
sx q[3];
rz(-2.9914732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.47485581) q[2];
sx q[2];
rz(-2.7754144) q[2];
sx q[2];
rz(0.22988698) q[2];
rz(2.722548) q[3];
sx q[3];
rz(-1.34904) q[3];
sx q[3];
rz(-0.45927799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6722365) q[0];
sx q[0];
rz(-2.3903963) q[0];
sx q[0];
rz(2.4647734) q[0];
rz(-2.6485486) q[1];
sx q[1];
rz(-2.1926011) q[1];
sx q[1];
rz(0.61606032) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6989243) q[0];
sx q[0];
rz(-1.6232423) q[0];
sx q[0];
rz(-3.1057182) q[0];
rz(-2.8370503) q[2];
sx q[2];
rz(-0.38123044) q[2];
sx q[2];
rz(2.9074557) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7459813) q[1];
sx q[1];
rz(-2.1484054) q[1];
sx q[1];
rz(-1.5773768) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89415278) q[3];
sx q[3];
rz(-1.812462) q[3];
sx q[3];
rz(1.5622996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.73413509) q[2];
sx q[2];
rz(-0.55441684) q[2];
sx q[2];
rz(-2.8862254) q[2];
rz(-1.5363961) q[3];
sx q[3];
rz(-2.1891749) q[3];
sx q[3];
rz(-2.3674964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42246321) q[0];
sx q[0];
rz(-0.96619773) q[0];
sx q[0];
rz(0.47250026) q[0];
rz(2.6155112) q[1];
sx q[1];
rz(-0.20985797) q[1];
sx q[1];
rz(2.2568259) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.103325) q[0];
sx q[0];
rz(-1.9511576) q[0];
sx q[0];
rz(3.0618219) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9315345) q[2];
sx q[2];
rz(-1.5569599) q[2];
sx q[2];
rz(-0.42523281) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6073508) q[1];
sx q[1];
rz(-2.631175) q[1];
sx q[1];
rz(0.95153248) q[1];
x q[2];
rz(0.24598908) q[3];
sx q[3];
rz(-1.0503328) q[3];
sx q[3];
rz(-0.99508475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3970967) q[2];
sx q[2];
rz(-0.8049736) q[2];
sx q[2];
rz(2.8302144) q[2];
rz(-1.7729676) q[3];
sx q[3];
rz(-0.45752782) q[3];
sx q[3];
rz(0.5293203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41806528) q[0];
sx q[0];
rz(-0.27856809) q[0];
sx q[0];
rz(0.061070651) q[0];
rz(-0.04018499) q[1];
sx q[1];
rz(-1.9804852) q[1];
sx q[1];
rz(-0.73289245) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7948579) q[0];
sx q[0];
rz(-1.96694) q[0];
sx q[0];
rz(-2.5651155) q[0];
x q[1];
rz(0.47809017) q[2];
sx q[2];
rz(-2.2149137) q[2];
sx q[2];
rz(-0.23020506) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1454754) q[1];
sx q[1];
rz(-0.62578326) q[1];
sx q[1];
rz(-1.5596418) q[1];
rz(-pi) q[2];
rz(1.93768) q[3];
sx q[3];
rz(-3.0266579) q[3];
sx q[3];
rz(0.68655187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0682893) q[2];
sx q[2];
rz(-1.1969593) q[2];
sx q[2];
rz(-0.51458365) q[2];
rz(1.2375281) q[3];
sx q[3];
rz(-0.5830183) q[3];
sx q[3];
rz(-2.5966743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.056203689) q[0];
sx q[0];
rz(-1.4852925) q[0];
sx q[0];
rz(-0.7094267) q[0];
rz(-1.5052694) q[1];
sx q[1];
rz(-2.067833) q[1];
sx q[1];
rz(2.8628796) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80517171) q[0];
sx q[0];
rz(-1.780691) q[0];
sx q[0];
rz(1.7091442) q[0];
x q[1];
rz(1.8144572) q[2];
sx q[2];
rz(-2.116034) q[2];
sx q[2];
rz(-1.3512163) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0540788) q[1];
sx q[1];
rz(-1.8808108) q[1];
sx q[1];
rz(2.9127496) q[1];
rz(-1.7940984) q[3];
sx q[3];
rz(-1.4947065) q[3];
sx q[3];
rz(3.1267816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8401106) q[2];
sx q[2];
rz(-1.8436517) q[2];
sx q[2];
rz(-0.056079496) q[2];
rz(0.85514832) q[3];
sx q[3];
rz(-2.6931098) q[3];
sx q[3];
rz(-2.7364031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1466325) q[0];
sx q[0];
rz(-2.9652847) q[0];
sx q[0];
rz(-2.1726998) q[0];
rz(-0.47337198) q[1];
sx q[1];
rz(-0.7946161) q[1];
sx q[1];
rz(1.999058) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8965217) q[0];
sx q[0];
rz(-2.8327496) q[0];
sx q[0];
rz(1.2678498) q[0];
rz(-0.70334401) q[2];
sx q[2];
rz(-1.5726657) q[2];
sx q[2];
rz(1.6714931) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7361703) q[1];
sx q[1];
rz(-2.1076964) q[1];
sx q[1];
rz(-1.7169397) q[1];
rz(2.9386018) q[3];
sx q[3];
rz(-0.91068017) q[3];
sx q[3];
rz(2.4329894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.896495) q[2];
sx q[2];
rz(-1.2197887) q[2];
sx q[2];
rz(-1.0207821) q[2];
rz(2.8178689) q[3];
sx q[3];
rz(-2.3886069) q[3];
sx q[3];
rz(-0.011172115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97994119) q[0];
sx q[0];
rz(-0.027898235) q[0];
sx q[0];
rz(-2.4401869) q[0];
rz(0.91570634) q[1];
sx q[1];
rz(-2.1332824) q[1];
sx q[1];
rz(-1.9030301) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0093875) q[0];
sx q[0];
rz(-1.3552249) q[0];
sx q[0];
rz(-1.042004) q[0];
rz(-1.5194703) q[2];
sx q[2];
rz(-2.3231299) q[2];
sx q[2];
rz(2.3527956) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.30675948) q[1];
sx q[1];
rz(-0.66654897) q[1];
sx q[1];
rz(1.761761) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11390399) q[3];
sx q[3];
rz(-2.1602727) q[3];
sx q[3];
rz(-1.9437499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.68676585) q[2];
sx q[2];
rz(-1.0964311) q[2];
sx q[2];
rz(-3.0977541) q[2];
rz(1.2010126) q[3];
sx q[3];
rz(-2.4062556) q[3];
sx q[3];
rz(-2.1380077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4713521) q[0];
sx q[0];
rz(-2.4155407) q[0];
sx q[0];
rz(1.7759905) q[0];
rz(3.1162221) q[1];
sx q[1];
rz(-1.3062968) q[1];
sx q[1];
rz(1.2702373) q[1];
rz(-0.30517898) q[2];
sx q[2];
rz(-1.1883493) q[2];
sx q[2];
rz(-2.3408008) q[2];
rz(-1.745789) q[3];
sx q[3];
rz(-2.3621515) q[3];
sx q[3];
rz(3.0099517) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];