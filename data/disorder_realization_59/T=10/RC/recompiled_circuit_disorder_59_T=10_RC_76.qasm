OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.74270785) q[0];
sx q[0];
rz(-0.78200114) q[0];
sx q[0];
rz(1.8703823) q[0];
rz(3.4186163) q[1];
sx q[1];
rz(3.613598) q[1];
sx q[1];
rz(9.4233905) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1479552) q[0];
sx q[0];
rz(-1.4076828) q[0];
sx q[0];
rz(1.9440584) q[0];
rz(-pi) q[1];
rz(2.6944461) q[2];
sx q[2];
rz(-1.5329754) q[2];
sx q[2];
rz(0.58085261) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4692987) q[1];
sx q[1];
rz(-2.0672332) q[1];
sx q[1];
rz(1.6507571) q[1];
x q[2];
rz(1.9418342) q[3];
sx q[3];
rz(-1.3317809) q[3];
sx q[3];
rz(2.1109076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9871621) q[2];
sx q[2];
rz(-0.61750948) q[2];
sx q[2];
rz(-0.74938613) q[2];
rz(-1.0162214) q[3];
sx q[3];
rz(-1.9640434) q[3];
sx q[3];
rz(-0.40482503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7063023) q[0];
sx q[0];
rz(-2.3162233) q[0];
sx q[0];
rz(-0.97066561) q[0];
rz(-2.1043815) q[1];
sx q[1];
rz(-1.4379921) q[1];
sx q[1];
rz(0.81545365) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97074189) q[0];
sx q[0];
rz(-2.4903957) q[0];
sx q[0];
rz(3.0424776) q[0];
x q[1];
rz(1.2977807) q[2];
sx q[2];
rz(-0.85999876) q[2];
sx q[2];
rz(3.0836011) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.26317877) q[1];
sx q[1];
rz(-2.8385332) q[1];
sx q[1];
rz(0.11942272) q[1];
x q[2];
rz(-2.4957982) q[3];
sx q[3];
rz(-1.7495973) q[3];
sx q[3];
rz(-1.4327232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4619535) q[2];
sx q[2];
rz(-1.5755499) q[2];
sx q[2];
rz(2.5088076) q[2];
rz(1.9880382) q[3];
sx q[3];
rz(-0.76806918) q[3];
sx q[3];
rz(0.30953428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84045029) q[0];
sx q[0];
rz(-1.9294894) q[0];
sx q[0];
rz(0.87483037) q[0];
rz(-1.8114999) q[1];
sx q[1];
rz(-1.7069838) q[1];
sx q[1];
rz(-0.99951807) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32854983) q[0];
sx q[0];
rz(-1.2721491) q[0];
sx q[0];
rz(-2.9850328) q[0];
rz(-pi) q[1];
rz(-2.3109762) q[2];
sx q[2];
rz(-2.3306371) q[2];
sx q[2];
rz(-0.22731552) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1516583) q[1];
sx q[1];
rz(-0.5608359) q[1];
sx q[1];
rz(2.6480617) q[1];
rz(-1.4286832) q[3];
sx q[3];
rz(-2.6283773) q[3];
sx q[3];
rz(1.3877447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.53753608) q[2];
sx q[2];
rz(-2.2194922) q[2];
sx q[2];
rz(-0.58004722) q[2];
rz(-0.81702685) q[3];
sx q[3];
rz(-1.3823119) q[3];
sx q[3];
rz(-1.1497315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.375305) q[0];
sx q[0];
rz(-1.5505318) q[0];
sx q[0];
rz(-0.91039175) q[0];
rz(0.45122775) q[1];
sx q[1];
rz(-1.5952361) q[1];
sx q[1];
rz(2.8667563) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9592322) q[0];
sx q[0];
rz(-1.4920456) q[0];
sx q[0];
rz(1.6596646) q[0];
x q[1];
rz(1.1103815) q[2];
sx q[2];
rz(-1.5693671) q[2];
sx q[2];
rz(-1.187385) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.53510016) q[1];
sx q[1];
rz(-2.1537158) q[1];
sx q[1];
rz(2.5938354) q[1];
rz(-pi) q[2];
rz(3.0148388) q[3];
sx q[3];
rz(-2.2634014) q[3];
sx q[3];
rz(2.9664489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.794902) q[2];
sx q[2];
rz(-1.1179504) q[2];
sx q[2];
rz(1.5083195) q[2];
rz(-1.9968962) q[3];
sx q[3];
rz(-0.73994023) q[3];
sx q[3];
rz(2.9798853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4836924) q[0];
sx q[0];
rz(-1.2150486) q[0];
sx q[0];
rz(-0.99779469) q[0];
rz(0.18355852) q[1];
sx q[1];
rz(-1.4869556) q[1];
sx q[1];
rz(1.516974) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.957513) q[0];
sx q[0];
rz(-0.74680579) q[0];
sx q[0];
rz(1.0190796) q[0];
rz(0.5337358) q[2];
sx q[2];
rz(-1.9118475) q[2];
sx q[2];
rz(1.595572) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6301873) q[1];
sx q[1];
rz(-1.9747707) q[1];
sx q[1];
rz(-1.9460815) q[1];
rz(2.153271) q[3];
sx q[3];
rz(-1.054793) q[3];
sx q[3];
rz(-2.2367246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8020442) q[2];
sx q[2];
rz(-2.1495154) q[2];
sx q[2];
rz(-0.3240164) q[2];
rz(-1.3230532) q[3];
sx q[3];
rz(-0.75606212) q[3];
sx q[3];
rz(-1.5312622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76535392) q[0];
sx q[0];
rz(-2.1026251) q[0];
sx q[0];
rz(-1.2639686) q[0];
rz(2.2309247) q[1];
sx q[1];
rz(-1.940454) q[1];
sx q[1];
rz(0.34067672) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7621988) q[0];
sx q[0];
rz(-1.5389991) q[0];
sx q[0];
rz(-1.5048774) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3855626) q[2];
sx q[2];
rz(-1.7871734) q[2];
sx q[2];
rz(-2.0331969) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50960474) q[1];
sx q[1];
rz(-1.0972411) q[1];
sx q[1];
rz(1.3118841) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86990279) q[3];
sx q[3];
rz(-2.5330336) q[3];
sx q[3];
rz(0.22939798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.95057758) q[2];
sx q[2];
rz(-0.58379972) q[2];
sx q[2];
rz(-0.77159709) q[2];
rz(2.5937882) q[3];
sx q[3];
rz(-2.1717725) q[3];
sx q[3];
rz(2.5781393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9724378) q[0];
sx q[0];
rz(-1.6945524) q[0];
sx q[0];
rz(-2.7959438) q[0];
rz(-3.0787643) q[1];
sx q[1];
rz(-2.6627916) q[1];
sx q[1];
rz(2.6766434) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5349605) q[0];
sx q[0];
rz(-2.6323942) q[0];
sx q[0];
rz(-1.1688933) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.99636997) q[2];
sx q[2];
rz(-1.1466221) q[2];
sx q[2];
rz(0.043957274) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.17864922) q[1];
sx q[1];
rz(-1.6766251) q[1];
sx q[1];
rz(2.9220198) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.025860272) q[3];
sx q[3];
rz(-1.6169294) q[3];
sx q[3];
rz(-0.24054724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0039625) q[2];
sx q[2];
rz(-1.5045065) q[2];
sx q[2];
rz(2.8239992) q[2];
rz(0.57146227) q[3];
sx q[3];
rz(-1.0390037) q[3];
sx q[3];
rz(-0.28731829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5193609) q[0];
sx q[0];
rz(-1.2943635) q[0];
sx q[0];
rz(0.28433329) q[0];
rz(2.590086) q[1];
sx q[1];
rz(-0.14177828) q[1];
sx q[1];
rz(-0.078358738) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1348304) q[0];
sx q[0];
rz(-1.5015331) q[0];
sx q[0];
rz(-0.80471054) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80348357) q[2];
sx q[2];
rz(-0.33499559) q[2];
sx q[2];
rz(-1.6840881) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8790508) q[1];
sx q[1];
rz(-1.7444376) q[1];
sx q[1];
rz(0.91211984) q[1];
x q[2];
rz(1.7275229) q[3];
sx q[3];
rz(-1.5369475) q[3];
sx q[3];
rz(-2.1050997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7408961) q[2];
sx q[2];
rz(-0.90831465) q[2];
sx q[2];
rz(0.25137869) q[2];
rz(2.5583983) q[3];
sx q[3];
rz(-2.0299032) q[3];
sx q[3];
rz(1.73197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3437929) q[0];
sx q[0];
rz(-3.0594337) q[0];
sx q[0];
rz(-3.0902241) q[0];
rz(0.92357606) q[1];
sx q[1];
rz(-2.4802465) q[1];
sx q[1];
rz(-2.267568) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4955935) q[0];
sx q[0];
rz(-0.80276239) q[0];
sx q[0];
rz(-1.5828703) q[0];
rz(-pi) q[1];
rz(0.85272312) q[2];
sx q[2];
rz(-2.6211779) q[2];
sx q[2];
rz(2.4783217) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.45174949) q[1];
sx q[1];
rz(-0.50536957) q[1];
sx q[1];
rz(1.7187353) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2285352) q[3];
sx q[3];
rz(-2.1517793) q[3];
sx q[3];
rz(2.5653605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41436568) q[2];
sx q[2];
rz(-0.7545158) q[2];
sx q[2];
rz(-2.5218463) q[2];
rz(-1.184458) q[3];
sx q[3];
rz(-1.8871566) q[3];
sx q[3];
rz(1.7782036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1353564) q[0];
sx q[0];
rz(-1.0422491) q[0];
sx q[0];
rz(-2.4172879) q[0];
rz(-0.18877098) q[1];
sx q[1];
rz(-2.962208) q[1];
sx q[1];
rz(-1.9627409) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9261242) q[0];
sx q[0];
rz(-1.7755531) q[0];
sx q[0];
rz(-1.8490851) q[0];
rz(1.4775425) q[2];
sx q[2];
rz(-1.4443195) q[2];
sx q[2];
rz(-1.7699514) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9465543) q[1];
sx q[1];
rz(-0.95609162) q[1];
sx q[1];
rz(2.9123995) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7614818) q[3];
sx q[3];
rz(-2.5862525) q[3];
sx q[3];
rz(1.7172608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0835691) q[2];
sx q[2];
rz(-1.0408137) q[2];
sx q[2];
rz(0.89938346) q[2];
rz(0.87456885) q[3];
sx q[3];
rz(-2.7159297) q[3];
sx q[3];
rz(-1.4609059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62109229) q[0];
sx q[0];
rz(-0.4091456) q[0];
sx q[0];
rz(-2.8950305) q[0];
rz(2.3836366) q[1];
sx q[1];
rz(-1.4823722) q[1];
sx q[1];
rz(-1.4588251) q[1];
rz(-1.4964676) q[2];
sx q[2];
rz(-2.7004514) q[2];
sx q[2];
rz(-2.2218291) q[2];
rz(-1.7189797) q[3];
sx q[3];
rz(-1.2978745) q[3];
sx q[3];
rz(0.10980448) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
