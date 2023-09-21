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
rz(-0.95681325) q[0];
sx q[0];
rz(1.4332888) q[0];
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
rz(1.0639609) q[0];
sx q[0];
rz(-1.6368757) q[0];
sx q[0];
rz(1.328701) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8561864) q[2];
sx q[2];
rz(-1.9754344) q[2];
sx q[2];
rz(2.0602496) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.45935985) q[1];
sx q[1];
rz(-1.1182922) q[1];
sx q[1];
rz(-0.020999055) q[1];
rz(-pi) q[2];
rz(-1.8688698) q[3];
sx q[3];
rz(-1.0469336) q[3];
sx q[3];
rz(3.1071172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0248489) q[2];
sx q[2];
rz(-1.3329093) q[2];
sx q[2];
rz(1.9072745) q[2];
rz(-2.0862789) q[3];
sx q[3];
rz(-2.0827259) q[3];
sx q[3];
rz(1.9887259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18594436) q[0];
sx q[0];
rz(-1.9444822) q[0];
sx q[0];
rz(2.3117075) q[0];
rz(2.8886967) q[1];
sx q[1];
rz(-2.3650832) q[1];
sx q[1];
rz(-0.84709644) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0711489) q[0];
sx q[0];
rz(-1.1304454) q[0];
sx q[0];
rz(-1.1308934) q[0];
rz(2.8424938) q[2];
sx q[2];
rz(-2.3397589) q[2];
sx q[2];
rz(-0.15660827) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6716799) q[1];
sx q[1];
rz(-2.5922853) q[1];
sx q[1];
rz(-3.0104907) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2291525) q[3];
sx q[3];
rz(-2.6077301) q[3];
sx q[3];
rz(-2.562059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0236686) q[2];
sx q[2];
rz(-0.37332049) q[2];
sx q[2];
rz(0.71933293) q[2];
rz(1.8524648) q[3];
sx q[3];
rz(-2.9057861) q[3];
sx q[3];
rz(1.6524564) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2704724) q[0];
sx q[0];
rz(-1.2501165) q[0];
sx q[0];
rz(0.57139325) q[0];
rz(-2.1748623) q[1];
sx q[1];
rz(-1.8833908) q[1];
sx q[1];
rz(0.5273231) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2741094) q[0];
sx q[0];
rz(-1.880078) q[0];
sx q[0];
rz(1.5779737) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2450561) q[2];
sx q[2];
rz(-2.2176761) q[2];
sx q[2];
rz(1.6270454) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.69818245) q[1];
sx q[1];
rz(-1.829362) q[1];
sx q[1];
rz(-2.3727388) q[1];
rz(2.4089291) q[3];
sx q[3];
rz(-1.817173) q[3];
sx q[3];
rz(2.8837567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8335235) q[2];
sx q[2];
rz(-1.6922733) q[2];
sx q[2];
rz(-2.4148338) q[2];
rz(2.1440992) q[3];
sx q[3];
rz(-2.5163429) q[3];
sx q[3];
rz(-1.6231026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8961287) q[0];
sx q[0];
rz(-1.4828869) q[0];
sx q[0];
rz(-1.0035275) q[0];
rz(-0.040680496) q[1];
sx q[1];
rz(-1.129312) q[1];
sx q[1];
rz(2.3153268) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6986615) q[0];
sx q[0];
rz(-2.1974265) q[0];
sx q[0];
rz(-0.71039623) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0777061) q[2];
sx q[2];
rz(-2.4404844) q[2];
sx q[2];
rz(2.137616) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.5845881) q[1];
sx q[1];
rz(-1.4950206) q[1];
sx q[1];
rz(1.2332548) q[1];
x q[2];
rz(-0.63185933) q[3];
sx q[3];
rz(-1.251289) q[3];
sx q[3];
rz(-2.7416122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5740009) q[2];
sx q[2];
rz(-1.8969798) q[2];
sx q[2];
rz(2.2272002) q[2];
rz(3.0363723) q[3];
sx q[3];
rz(-1.5172232) q[3];
sx q[3];
rz(-2.5523394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2545664) q[0];
sx q[0];
rz(-1.5019324) q[0];
sx q[0];
rz(2.0078833) q[0];
rz(0.0056313593) q[1];
sx q[1];
rz(-1.4375604) q[1];
sx q[1];
rz(2.5240135) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3546346) q[0];
sx q[0];
rz(-1.5411396) q[0];
sx q[0];
rz(-1.5251072) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1281625) q[2];
sx q[2];
rz(-0.97634041) q[2];
sx q[2];
rz(-2.9885459) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6981814) q[1];
sx q[1];
rz(-2.4907506) q[1];
sx q[1];
rz(-1.0990259) q[1];
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
rz(-0.99308333) q[3];
sx q[3];
rz(-2.1916094) q[3];
sx q[3];
rz(0.63794199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1722906) q[0];
sx q[0];
rz(-3.0811221) q[0];
sx q[0];
rz(-0.0897952) q[0];
rz(0.88109294) q[1];
sx q[1];
rz(-0.52905622) q[1];
sx q[1];
rz(-2.9615013) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20443944) q[0];
sx q[0];
rz(-2.3817872) q[0];
sx q[0];
rz(2.252749) q[0];
rz(-pi) q[1];
rz(-0.96713709) q[2];
sx q[2];
rz(-1.4146311) q[2];
sx q[2];
rz(-0.16317633) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0568911) q[1];
sx q[1];
rz(-1.7329721) q[1];
sx q[1];
rz(-2.8497036) q[1];
x q[2];
rz(-1.5363541) q[3];
sx q[3];
rz(-0.4930217) q[3];
sx q[3];
rz(0.29355129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.2580516) q[2];
sx q[2];
rz(-1.5641314) q[2];
sx q[2];
rz(-1.9656666) q[2];
rz(-0.073143395) q[3];
sx q[3];
rz(-0.72435838) q[3];
sx q[3];
rz(0.49889645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.183855) q[0];
sx q[0];
rz(-2.4860005) q[0];
sx q[0];
rz(-1.7234329) q[0];
rz(-1.9765967) q[1];
sx q[1];
rz(-0.55924758) q[1];
sx q[1];
rz(0.040239008) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2049853) q[0];
sx q[0];
rz(-1.4892704) q[0];
sx q[0];
rz(-1.3315014) q[0];
x q[1];
rz(3.0954103) q[2];
sx q[2];
rz(-1.1650411) q[2];
sx q[2];
rz(1.501776) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.38899598) q[1];
sx q[1];
rz(-0.95655671) q[1];
sx q[1];
rz(-0.1715626) q[1];
x q[2];
rz(1.8764898) q[3];
sx q[3];
rz(-1.2672658) q[3];
sx q[3];
rz(-2.2895253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.13850257) q[2];
sx q[2];
rz(-0.78531992) q[2];
sx q[2];
rz(-2.9677532) q[2];
rz(1.7447757) q[3];
sx q[3];
rz(-1.6766179) q[3];
sx q[3];
rz(-0.85913908) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1180856) q[0];
sx q[0];
rz(-0.21723391) q[0];
sx q[0];
rz(-1.404495) q[0];
rz(1.2069758) q[1];
sx q[1];
rz(-1.4852306) q[1];
sx q[1];
rz(-1.6361902) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0508142) q[0];
sx q[0];
rz(-2.6065126) q[0];
sx q[0];
rz(-2.4909766) q[0];
rz(-pi) q[1];
rz(0.6263528) q[2];
sx q[2];
rz(-2.2144631) q[2];
sx q[2];
rz(2.6477637) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7561188) q[1];
sx q[1];
rz(-1.6033152) q[1];
sx q[1];
rz(-0.69617747) q[1];
x q[2];
rz(1.1106311) q[3];
sx q[3];
rz(-0.65350973) q[3];
sx q[3];
rz(-1.489952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3796842) q[2];
sx q[2];
rz(-0.48019335) q[2];
sx q[2];
rz(-2.9910679) q[2];
rz(-1.5395509) q[3];
sx q[3];
rz(-1.9463041) q[3];
sx q[3];
rz(-0.62301821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4432916) q[0];
sx q[0];
rz(-1.1346096) q[0];
sx q[0];
rz(-2.495893) q[0];
rz(-0.49939108) q[1];
sx q[1];
rz(-1.4285409) q[1];
sx q[1];
rz(-1.2649149) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3919968) q[0];
sx q[0];
rz(-2.2889334) q[0];
sx q[0];
rz(2.0977661) q[0];
rz(-pi) q[1];
rz(0.64847704) q[2];
sx q[2];
rz(-0.65995526) q[2];
sx q[2];
rz(-2.6476268) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.39667323) q[1];
sx q[1];
rz(-1.5349689) q[1];
sx q[1];
rz(-2.9958535) q[1];
rz(-pi) q[2];
rz(-0.42322741) q[3];
sx q[3];
rz(-1.755135) q[3];
sx q[3];
rz(1.1545899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.634793) q[2];
sx q[2];
rz(-0.96968499) q[2];
sx q[2];
rz(-2.6169422) q[2];
rz(-0.55650416) q[3];
sx q[3];
rz(-1.6289214) q[3];
sx q[3];
rz(-2.7867253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5234579) q[0];
sx q[0];
rz(-2.3045461) q[0];
sx q[0];
rz(0.24542228) q[0];
rz(0.55631176) q[1];
sx q[1];
rz(-1.5309155) q[1];
sx q[1];
rz(1.6773178) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2357764) q[0];
sx q[0];
rz(-2.896744) q[0];
sx q[0];
rz(-2.66923) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1572838) q[2];
sx q[2];
rz(-1.7585635) q[2];
sx q[2];
rz(-0.078660065) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.163584) q[1];
sx q[1];
rz(-0.89466909) q[1];
sx q[1];
rz(-2.7917557) q[1];
x q[2];
rz(-0.085300307) q[3];
sx q[3];
rz(-1.395426) q[3];
sx q[3];
rz(-2.8773957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2422553) q[2];
sx q[2];
rz(-1.238404) q[2];
sx q[2];
rz(0.74550068) q[2];
rz(-1.4238822) q[3];
sx q[3];
rz(-2.8427377) q[3];
sx q[3];
rz(2.0132813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8650919) q[0];
sx q[0];
rz(-1.4958953) q[0];
sx q[0];
rz(-2.4687742) q[0];
rz(-1.2561692) q[1];
sx q[1];
rz(-2.3354463) q[1];
sx q[1];
rz(-1.068402) q[1];
rz(-2.4297498) q[2];
sx q[2];
rz(-1.9528452) q[2];
sx q[2];
rz(-0.14609329) q[2];
rz(-1.9477378) q[3];
sx q[3];
rz(-2.1115163) q[3];
sx q[3];
rz(1.9797309) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
