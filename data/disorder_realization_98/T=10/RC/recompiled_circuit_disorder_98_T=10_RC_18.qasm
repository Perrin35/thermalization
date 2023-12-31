OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7712819) q[0];
sx q[0];
rz(-0.9960649) q[0];
sx q[0];
rz(2.2709742) q[0];
rz(-1.0215966) q[1];
sx q[1];
rz(-0.28290132) q[1];
sx q[1];
rz(-0.14970782) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71404845) q[0];
sx q[0];
rz(-0.85435003) q[0];
sx q[0];
rz(-2.2358405) q[0];
rz(-pi) q[1];
rz(-2.6354191) q[2];
sx q[2];
rz(-1.0867599) q[2];
sx q[2];
rz(-1.6137705) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.61923164) q[1];
sx q[1];
rz(-1.751096) q[1];
sx q[1];
rz(0.038832263) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1211795) q[3];
sx q[3];
rz(-1.4074416) q[3];
sx q[3];
rz(2.6640716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.25847882) q[2];
sx q[2];
rz(-3.0443865) q[2];
sx q[2];
rz(-0.70409888) q[2];
rz(-0.95300931) q[3];
sx q[3];
rz(-2.1694031) q[3];
sx q[3];
rz(-1.7378418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54884058) q[0];
sx q[0];
rz(-1.623818) q[0];
sx q[0];
rz(0.63252226) q[0];
rz(-0.44644341) q[1];
sx q[1];
rz(-1.4182785) q[1];
sx q[1];
rz(2.4893563) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3632293) q[0];
sx q[0];
rz(-0.031154545) q[0];
sx q[0];
rz(-2.7345783) q[0];
rz(-pi) q[1];
rz(-1.863443) q[2];
sx q[2];
rz(-1.7990944) q[2];
sx q[2];
rz(1.4676859) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.80174175) q[1];
sx q[1];
rz(-2.7738214) q[1];
sx q[1];
rz(2.4778609) q[1];
x q[2];
rz(0.60450508) q[3];
sx q[3];
rz(-2.1406056) q[3];
sx q[3];
rz(-1.9250211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5474881) q[2];
sx q[2];
rz(-1.0141806) q[2];
sx q[2];
rz(1.9799505) q[2];
rz(1.9836327) q[3];
sx q[3];
rz(-2.0722814) q[3];
sx q[3];
rz(-1.6903711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050425477) q[0];
sx q[0];
rz(-2.7624891) q[0];
sx q[0];
rz(-0.85025775) q[0];
rz(0.49750528) q[1];
sx q[1];
rz(-1.9583227) q[1];
sx q[1];
rz(1.7920378) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2230167) q[0];
sx q[0];
rz(-2.4239459) q[0];
sx q[0];
rz(0.36803228) q[0];
rz(-1.0565874) q[2];
sx q[2];
rz(-3.0658256) q[2];
sx q[2];
rz(0.55759341) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.68223665) q[1];
sx q[1];
rz(-0.56486928) q[1];
sx q[1];
rz(-0.45046803) q[1];
x q[2];
rz(0.0098185929) q[3];
sx q[3];
rz(-1.9968642) q[3];
sx q[3];
rz(2.1992418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1094018) q[2];
sx q[2];
rz(-1.9058062) q[2];
sx q[2];
rz(-2.2375977) q[2];
rz(-2.8404625) q[3];
sx q[3];
rz(-1.7588994) q[3];
sx q[3];
rz(1.7416471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6501453) q[0];
sx q[0];
rz(-0.97147816) q[0];
sx q[0];
rz(-1.7310463) q[0];
rz(2.5097805) q[1];
sx q[1];
rz(-1.8099064) q[1];
sx q[1];
rz(-0.036380336) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.13658) q[0];
sx q[0];
rz(-2.4653325) q[0];
sx q[0];
rz(-1.9146634) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6668947) q[2];
sx q[2];
rz(-2.5050852) q[2];
sx q[2];
rz(2.0091025) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7155647) q[1];
sx q[1];
rz(-1.702311) q[1];
sx q[1];
rz(-0.87042602) q[1];
x q[2];
rz(-2.9070204) q[3];
sx q[3];
rz(-1.6664701) q[3];
sx q[3];
rz(2.5535339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2146384) q[2];
sx q[2];
rz(-1.4321233) q[2];
sx q[2];
rz(1.1882163) q[2];
rz(0.67048091) q[3];
sx q[3];
rz(-1.2256349) q[3];
sx q[3];
rz(-0.59613434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7017512) q[0];
sx q[0];
rz(-1.8816467) q[0];
sx q[0];
rz(2.9751076) q[0];
rz(2.3855551) q[1];
sx q[1];
rz(-0.88587228) q[1];
sx q[1];
rz(-2.9072445) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34412947) q[0];
sx q[0];
rz(-1.2327694) q[0];
sx q[0];
rz(-2.5812838) q[0];
x q[1];
rz(-1.7469823) q[2];
sx q[2];
rz(-0.50054769) q[2];
sx q[2];
rz(-1.2622152) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.105141) q[1];
sx q[1];
rz(-2.3804133) q[1];
sx q[1];
rz(0.21703227) q[1];
rz(-0.078951051) q[3];
sx q[3];
rz(-1.8936994) q[3];
sx q[3];
rz(-0.80432804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4328737) q[2];
sx q[2];
rz(-1.1756228) q[2];
sx q[2];
rz(-0.22053545) q[2];
rz(-2.7045414) q[3];
sx q[3];
rz(-2.1190937) q[3];
sx q[3];
rz(2.3760858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41546145) q[0];
sx q[0];
rz(-2.8227865) q[0];
sx q[0];
rz(0.81714001) q[0];
rz(0.56610402) q[1];
sx q[1];
rz(-1.348446) q[1];
sx q[1];
rz(-1.1436499) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7276579) q[0];
sx q[0];
rz(-1.2000788) q[0];
sx q[0];
rz(1.6137705) q[0];
rz(-pi) q[1];
rz(1.3104865) q[2];
sx q[2];
rz(-2.227265) q[2];
sx q[2];
rz(1.3442163) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1144048) q[1];
sx q[1];
rz(-1.4804375) q[1];
sx q[1];
rz(-2.0143709) q[1];
rz(-pi) q[2];
rz(-2.6753622) q[3];
sx q[3];
rz(-2.5737408) q[3];
sx q[3];
rz(-0.0044435244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3879261) q[2];
sx q[2];
rz(-1.5619229) q[2];
sx q[2];
rz(2.6521519) q[2];
rz(-0.22805452) q[3];
sx q[3];
rz(-1.2585879) q[3];
sx q[3];
rz(0.42603809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.56931) q[0];
sx q[0];
rz(-2.4991878) q[0];
sx q[0];
rz(1.2868767) q[0];
rz(-2.4781748) q[1];
sx q[1];
rz(-1.5692915) q[1];
sx q[1];
rz(-1.2333262) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7466) q[0];
sx q[0];
rz(-1.9548423) q[0];
sx q[0];
rz(-3.1296484) q[0];
rz(-pi) q[1];
rz(0.6255409) q[2];
sx q[2];
rz(-0.95226804) q[2];
sx q[2];
rz(2.2369838) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7066321) q[1];
sx q[1];
rz(-3.0295277) q[1];
sx q[1];
rz(-0.12607615) q[1];
rz(-pi) q[2];
rz(2.535378) q[3];
sx q[3];
rz(-2.5698235) q[3];
sx q[3];
rz(-1.1796463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.037501637) q[2];
sx q[2];
rz(-2.9064894) q[2];
sx q[2];
rz(2.1255778) q[2];
rz(-3.0715023) q[3];
sx q[3];
rz(-1.9349808) q[3];
sx q[3];
rz(2.0751374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12524097) q[0];
sx q[0];
rz(-2.4588983) q[0];
sx q[0];
rz(-1.4455147) q[0];
rz(0.21487543) q[1];
sx q[1];
rz(-2.3863249) q[1];
sx q[1];
rz(-1.8833556) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1930775) q[0];
sx q[0];
rz(-2.6876039) q[0];
sx q[0];
rz(-1.8817188) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3636742) q[2];
sx q[2];
rz(-0.81233835) q[2];
sx q[2];
rz(1.1493491) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3360577) q[1];
sx q[1];
rz(-1.6748669) q[1];
sx q[1];
rz(-0.54688262) q[1];
x q[2];
rz(-2.1636837) q[3];
sx q[3];
rz(-0.80855723) q[3];
sx q[3];
rz(0.18329328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4593279) q[2];
sx q[2];
rz(-2.9569646) q[2];
sx q[2];
rz(-2.5788467) q[2];
rz(2.941926) q[3];
sx q[3];
rz(-0.66185799) q[3];
sx q[3];
rz(1.1220042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0257618) q[0];
sx q[0];
rz(-2.0985726) q[0];
sx q[0];
rz(-0.2510221) q[0];
rz(-2.714278) q[1];
sx q[1];
rz(-1.9117833) q[1];
sx q[1];
rz(3.1138611) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32556191) q[0];
sx q[0];
rz(-2.0566018) q[0];
sx q[0];
rz(-2.3150139) q[0];
rz(-0.26720033) q[2];
sx q[2];
rz(-2.2520817) q[2];
sx q[2];
rz(-0.53182488) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9412649) q[1];
sx q[1];
rz(-0.97201921) q[1];
sx q[1];
rz(1.4501249) q[1];
x q[2];
rz(-1.5042217) q[3];
sx q[3];
rz(-1.5449617) q[3];
sx q[3];
rz(-1.0915826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.015908265) q[2];
sx q[2];
rz(-2.1618844) q[2];
sx q[2];
rz(-1.1431747) q[2];
rz(-2.9987191) q[3];
sx q[3];
rz(-1.5121127) q[3];
sx q[3];
rz(-0.85723248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7509572) q[0];
sx q[0];
rz(-1.9366783) q[0];
sx q[0];
rz(-0.45387682) q[0];
rz(0.67165309) q[1];
sx q[1];
rz(-1.6891054) q[1];
sx q[1];
rz(0.25751105) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47782183) q[0];
sx q[0];
rz(-1.1241233) q[0];
sx q[0];
rz(-1.1429943) q[0];
rz(-pi) q[1];
rz(2.4117878) q[2];
sx q[2];
rz(-0.68241718) q[2];
sx q[2];
rz(2.6953816) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.232302) q[1];
sx q[1];
rz(-1.5564939) q[1];
sx q[1];
rz(-2.4453352) q[1];
rz(1.0912283) q[3];
sx q[3];
rz(-1.7864831) q[3];
sx q[3];
rz(-1.2558162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8132849) q[2];
sx q[2];
rz(-1.4537145) q[2];
sx q[2];
rz(-0.60662398) q[2];
rz(0.47484067) q[3];
sx q[3];
rz(-2.1947221) q[3];
sx q[3];
rz(2.301208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7941147) q[0];
sx q[0];
rz(-1.5119727) q[0];
sx q[0];
rz(2.150362) q[0];
rz(0.22656245) q[1];
sx q[1];
rz(-1.4145874) q[1];
sx q[1];
rz(0.58691595) q[1];
rz(-1.1889585) q[2];
sx q[2];
rz(-2.0094064) q[2];
sx q[2];
rz(-1.1878427) q[2];
rz(-2.0541035) q[3];
sx q[3];
rz(-1.8046422) q[3];
sx q[3];
rz(1.3840152) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
