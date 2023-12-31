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
rz(0.2431915) q[1];
sx q[1];
rz(4.3742124) q[1];
sx q[1];
rz(10.32962) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63431595) q[0];
sx q[0];
rz(-2.0024558) q[0];
sx q[0];
rz(-2.1202203) q[0];
rz(-pi) q[1];
rz(-0.71360795) q[2];
sx q[2];
rz(-0.2954233) q[2];
sx q[2];
rz(-1.7507391) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.79996586) q[1];
sx q[1];
rz(-0.46712616) q[1];
sx q[1];
rz(-1.0270236) q[1];
rz(-pi) q[2];
rz(-2.1078029) q[3];
sx q[3];
rz(-2.783943) q[3];
sx q[3];
rz(1.3057749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.20377542) q[2];
sx q[2];
rz(-1.4062466) q[2];
sx q[2];
rz(-3.0453483) q[2];
rz(-2.105666) q[3];
sx q[3];
rz(-2.7544498) q[3];
sx q[3];
rz(-2.9878785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44089833) q[0];
sx q[0];
rz(-0.39114025) q[0];
sx q[0];
rz(2.3764215) q[0];
rz(-1.8493429) q[1];
sx q[1];
rz(-0.48520979) q[1];
sx q[1];
rz(0.66295019) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63755858) q[0];
sx q[0];
rz(-3.0507037) q[0];
sx q[0];
rz(-0.81764098) q[0];
x q[1];
rz(2.2234369) q[2];
sx q[2];
rz(-1.895972) q[2];
sx q[2];
rz(-2.6174389) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.029164974) q[1];
sx q[1];
rz(-1.6532835) q[1];
sx q[1];
rz(-2.6998181) q[1];
rz(2.8421721) q[3];
sx q[3];
rz(-2.6357108) q[3];
sx q[3];
rz(0.29380709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.092676) q[2];
sx q[2];
rz(-1.0485704) q[2];
sx q[2];
rz(-0.32901397) q[2];
rz(2.4760903) q[3];
sx q[3];
rz(-0.21829675) q[3];
sx q[3];
rz(1.8238508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5488141) q[0];
sx q[0];
rz(-0.22286649) q[0];
sx q[0];
rz(-0.22234017) q[0];
rz(1.0173343) q[1];
sx q[1];
rz(-2.4203114) q[1];
sx q[1];
rz(0.51868784) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1919353) q[0];
sx q[0];
rz(-1.9266832) q[0];
sx q[0];
rz(0.8215254) q[0];
x q[1];
rz(-1.475004) q[2];
sx q[2];
rz(-1.8358069) q[2];
sx q[2];
rz(-1.4694627) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0842523) q[1];
sx q[1];
rz(-2.3358314) q[1];
sx q[1];
rz(-1.8112103) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1268901) q[3];
sx q[3];
rz(-0.054617453) q[3];
sx q[3];
rz(0.86044776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8609994) q[2];
sx q[2];
rz(-0.36182797) q[2];
sx q[2];
rz(-3.0730491) q[2];
rz(0.60244256) q[3];
sx q[3];
rz(-0.76255637) q[3];
sx q[3];
rz(3.0025735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4908726) q[0];
sx q[0];
rz(-2.3644709) q[0];
sx q[0];
rz(-0.17424507) q[0];
rz(2.6113367) q[1];
sx q[1];
rz(-1.4825876) q[1];
sx q[1];
rz(0.51309103) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4671191) q[0];
sx q[0];
rz(-1.4917443) q[0];
sx q[0];
rz(2.5208958) q[0];
rz(-pi) q[1];
rz(2.1999173) q[2];
sx q[2];
rz(-0.86949124) q[2];
sx q[2];
rz(-2.6710682) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7375813) q[1];
sx q[1];
rz(-1.1176795) q[1];
sx q[1];
rz(1.416942) q[1];
x q[2];
rz(3.1107535) q[3];
sx q[3];
rz(-1.3645002) q[3];
sx q[3];
rz(1.7145969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6667368) q[2];
sx q[2];
rz(-2.7754144) q[2];
sx q[2];
rz(0.22988698) q[2];
rz(2.722548) q[3];
sx q[3];
rz(-1.34904) q[3];
sx q[3];
rz(2.6823147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4693562) q[0];
sx q[0];
rz(-0.75119632) q[0];
sx q[0];
rz(2.4647734) q[0];
rz(-2.6485486) q[1];
sx q[1];
rz(-2.1926011) q[1];
sx q[1];
rz(0.61606032) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6989243) q[0];
sx q[0];
rz(-1.5183503) q[0];
sx q[0];
rz(3.1057182) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4511756) q[2];
sx q[2];
rz(-1.9336485) q[2];
sx q[2];
rz(-2.5809443) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9628145) q[1];
sx q[1];
rz(-1.5763092) q[1];
sx q[1];
rz(0.57761901) q[1];
rz(-pi) q[2];
rz(1.945799) q[3];
sx q[3];
rz(-2.4295394) q[3];
sx q[3];
rz(-0.28111162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.73413509) q[2];
sx q[2];
rz(-0.55441684) q[2];
sx q[2];
rz(2.8862254) q[2];
rz(-1.6051965) q[3];
sx q[3];
rz(-0.95241773) q[3];
sx q[3];
rz(-2.3674964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7191294) q[0];
sx q[0];
rz(-2.1753949) q[0];
sx q[0];
rz(-0.47250026) q[0];
rz(-2.6155112) q[1];
sx q[1];
rz(-0.20985797) q[1];
sx q[1];
rz(-2.2568259) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5793943) q[0];
sx q[0];
rz(-1.644855) q[0];
sx q[0];
rz(1.9522569) q[0];
rz(0.066263513) q[2];
sx q[2];
rz(-0.21050669) q[2];
sx q[2];
rz(-1.9312242) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2192167) q[1];
sx q[1];
rz(-1.9799385) q[1];
sx q[1];
rz(-2.8273696) q[1];
x q[2];
rz(-2.8956036) q[3];
sx q[3];
rz(-2.0912598) q[3];
sx q[3];
rz(-2.1465079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3970967) q[2];
sx q[2];
rz(-0.8049736) q[2];
sx q[2];
rz(-2.8302144) q[2];
rz(-1.3686251) q[3];
sx q[3];
rz(-0.45752782) q[3];
sx q[3];
rz(2.6122724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7235274) q[0];
sx q[0];
rz(-2.8630246) q[0];
sx q[0];
rz(0.061070651) q[0];
rz(3.1014077) q[1];
sx q[1];
rz(-1.9804852) q[1];
sx q[1];
rz(2.4087002) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7948579) q[0];
sx q[0];
rz(-1.1746527) q[0];
sx q[0];
rz(-0.57647716) q[0];
rz(0.86874666) q[2];
sx q[2];
rz(-1.947543) q[2];
sx q[2];
rz(1.6422611) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1454754) q[1];
sx q[1];
rz(-2.5158094) q[1];
sx q[1];
rz(-1.5596418) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4634499) q[3];
sx q[3];
rz(-1.6119453) q[3];
sx q[3];
rz(-2.6220208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0733033) q[2];
sx q[2];
rz(-1.9446334) q[2];
sx q[2];
rz(-2.627009) q[2];
rz(-1.2375281) q[3];
sx q[3];
rz(-0.5830183) q[3];
sx q[3];
rz(-0.54491836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-1.6363232) q[1];
sx q[1];
rz(-2.067833) q[1];
sx q[1];
rz(0.27871305) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80517171) q[0];
sx q[0];
rz(-1.3609017) q[0];
sx q[0];
rz(-1.4324485) q[0];
x q[1];
rz(-1.8144572) q[2];
sx q[2];
rz(-2.116034) q[2];
sx q[2];
rz(1.3512163) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0875138) q[1];
sx q[1];
rz(-1.2607818) q[1];
sx q[1];
rz(-2.9127496) q[1];
x q[2];
rz(1.9023499) q[3];
sx q[3];
rz(-2.9058876) q[3];
sx q[3];
rz(-1.9086259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8401106) q[2];
sx q[2];
rz(-1.297941) q[2];
sx q[2];
rz(-3.0855132) q[2];
rz(0.85514832) q[3];
sx q[3];
rz(-0.44848281) q[3];
sx q[3];
rz(2.7364031) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99496019) q[0];
sx q[0];
rz(-2.9652847) q[0];
sx q[0];
rz(-2.1726998) q[0];
rz(0.47337198) q[1];
sx q[1];
rz(-2.3469766) q[1];
sx q[1];
rz(1.999058) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8965217) q[0];
sx q[0];
rz(-2.8327496) q[0];
sx q[0];
rz(-1.8737428) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5683453) q[2];
sx q[2];
rz(-2.2741389) q[2];
sx q[2];
rz(0.0991115) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7361703) q[1];
sx q[1];
rz(-1.0338963) q[1];
sx q[1];
rz(-1.7169397) q[1];
x q[2];
rz(-0.20299083) q[3];
sx q[3];
rz(-2.2309125) q[3];
sx q[3];
rz(-2.4329894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.896495) q[2];
sx q[2];
rz(-1.2197887) q[2];
sx q[2];
rz(2.1208105) q[2];
rz(-2.8178689) q[3];
sx q[3];
rz(-0.75298572) q[3];
sx q[3];
rz(-0.011172115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.1616515) q[0];
sx q[0];
rz(-0.027898235) q[0];
sx q[0];
rz(-0.7014057) q[0];
rz(0.91570634) q[1];
sx q[1];
rz(-2.1332824) q[1];
sx q[1];
rz(-1.9030301) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1322051) q[0];
sx q[0];
rz(-1.3552249) q[0];
sx q[0];
rz(-2.0995887) q[0];
rz(-pi) q[1];
rz(1.5194703) q[2];
sx q[2];
rz(-2.3231299) q[2];
sx q[2];
rz(-2.3527956) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.30675948) q[1];
sx q[1];
rz(-0.66654897) q[1];
sx q[1];
rz(-1.761761) q[1];
rz(-0.978312) q[3];
sx q[3];
rz(-1.6654135) q[3];
sx q[3];
rz(-2.832151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.68676585) q[2];
sx q[2];
rz(-1.0964311) q[2];
sx q[2];
rz(-0.043838538) q[2];
rz(1.94058) q[3];
sx q[3];
rz(-0.73533708) q[3];
sx q[3];
rz(-2.1380077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6702406) q[0];
sx q[0];
rz(-0.72605194) q[0];
sx q[0];
rz(-1.3656021) q[0];
rz(3.1162221) q[1];
sx q[1];
rz(-1.3062968) q[1];
sx q[1];
rz(1.2702373) q[1];
rz(0.92924835) q[2];
sx q[2];
rz(-2.6570448) q[2];
sx q[2];
rz(-1.6397283) q[2];
rz(-0.79904859) q[3];
sx q[3];
rz(-1.4481164) q[3];
sx q[3];
rz(-1.8275402) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
