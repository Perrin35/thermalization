OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.030169686) q[0];
sx q[0];
rz(-1.7813959) q[0];
sx q[0];
rz(2.7786322) q[0];
rz(-1.1765923) q[1];
sx q[1];
rz(-1.7791553) q[1];
sx q[1];
rz(-1.4443614) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1189277) q[0];
sx q[0];
rz(-1.3937409) q[0];
sx q[0];
rz(-1.6657959) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8536804) q[2];
sx q[2];
rz(-1.4511564) q[2];
sx q[2];
rz(1.5515012) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1361268) q[1];
sx q[1];
rz(-2.695752) q[1];
sx q[1];
rz(1.7923097) q[1];
x q[2];
rz(1.6300171) q[3];
sx q[3];
rz(-0.60225669) q[3];
sx q[3];
rz(-1.8730375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1732257) q[2];
sx q[2];
rz(-0.78725615) q[2];
sx q[2];
rz(-0.095189722) q[2];
rz(1.7732874) q[3];
sx q[3];
rz(-0.94822001) q[3];
sx q[3];
rz(2.8424954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66265166) q[0];
sx q[0];
rz(-2.3421685) q[0];
sx q[0];
rz(2.4484402) q[0];
rz(0.33723801) q[1];
sx q[1];
rz(-1.8353029) q[1];
sx q[1];
rz(0.79280028) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1900729) q[0];
sx q[0];
rz(-1.8197104) q[0];
sx q[0];
rz(-2.3547257) q[0];
rz(-pi) q[1];
rz(2.9117965) q[2];
sx q[2];
rz(-1.0655128) q[2];
sx q[2];
rz(1.2838703) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2094272) q[1];
sx q[1];
rz(-1.8575728) q[1];
sx q[1];
rz(2.3580103) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6565142) q[3];
sx q[3];
rz(-2.2612319) q[3];
sx q[3];
rz(2.9196965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1482131) q[2];
sx q[2];
rz(-0.17017636) q[2];
sx q[2];
rz(-1.4519838) q[2];
rz(-2.4827237) q[3];
sx q[3];
rz(-1.4631319) q[3];
sx q[3];
rz(0.39670593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0186998) q[0];
sx q[0];
rz(-1.0188228) q[0];
sx q[0];
rz(0.10957154) q[0];
rz(2.295629) q[1];
sx q[1];
rz(-0.94596091) q[1];
sx q[1];
rz(0.74839655) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1066662) q[0];
sx q[0];
rz(-1.0163808) q[0];
sx q[0];
rz(3.0361452) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0399299) q[2];
sx q[2];
rz(-1.0541087) q[2];
sx q[2];
rz(-2.6831107) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3378395) q[1];
sx q[1];
rz(-3.0782135) q[1];
sx q[1];
rz(-1.4638639) q[1];
rz(-pi) q[2];
rz(-0.99936463) q[3];
sx q[3];
rz(-1.4203705) q[3];
sx q[3];
rz(2.5932464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1154068) q[2];
sx q[2];
rz(-1.7188965) q[2];
sx q[2];
rz(-2.4288948) q[2];
rz(-1.4767492) q[3];
sx q[3];
rz(-1.289117) q[3];
sx q[3];
rz(-3.0570928) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0067714) q[0];
sx q[0];
rz(-1.1399784) q[0];
sx q[0];
rz(0.086932927) q[0];
rz(0.25761071) q[1];
sx q[1];
rz(-2.0517709) q[1];
sx q[1];
rz(1.849252) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9168388) q[0];
sx q[0];
rz(-1.8096905) q[0];
sx q[0];
rz(-0.80015667) q[0];
rz(-pi) q[1];
rz(1.7316225) q[2];
sx q[2];
rz(-1.2431182) q[2];
sx q[2];
rz(0.92152061) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7304157) q[1];
sx q[1];
rz(-2.431793) q[1];
sx q[1];
rz(-1.0214772) q[1];
rz(-1.4216018) q[3];
sx q[3];
rz(-1.215862) q[3];
sx q[3];
rz(-1.9082663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8452235) q[2];
sx q[2];
rz(-1.1901647) q[2];
sx q[2];
rz(2.2947218) q[2];
rz(1.752468) q[3];
sx q[3];
rz(-1.4877157) q[3];
sx q[3];
rz(2.5510767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1883989) q[0];
sx q[0];
rz(-1.1253091) q[0];
sx q[0];
rz(-0.97287792) q[0];
rz(0.94213048) q[1];
sx q[1];
rz(-0.82754389) q[1];
sx q[1];
rz(-3.1408302) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89524901) q[0];
sx q[0];
rz(-2.0813119) q[0];
sx q[0];
rz(0.58175485) q[0];
rz(-1.23998) q[2];
sx q[2];
rz(-2.1789765) q[2];
sx q[2];
rz(1.1229148) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9474802) q[1];
sx q[1];
rz(-1.0092217) q[1];
sx q[1];
rz(0.33163867) q[1];
rz(-pi) q[2];
rz(-1.6490178) q[3];
sx q[3];
rz(-0.97710246) q[3];
sx q[3];
rz(-1.2562615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3804649) q[2];
sx q[2];
rz(-1.7385812) q[2];
sx q[2];
rz(-0.38499704) q[2];
rz(0.72832406) q[3];
sx q[3];
rz(-2.5326122) q[3];
sx q[3];
rz(0.42496067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4094792) q[0];
sx q[0];
rz(-2.3347169) q[0];
sx q[0];
rz(-0.4581067) q[0];
rz(1.4234281) q[1];
sx q[1];
rz(-0.79045311) q[1];
sx q[1];
rz(-0.093553392) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6038743) q[0];
sx q[0];
rz(-0.97386003) q[0];
sx q[0];
rz(0.56054795) q[0];
x q[1];
rz(0.16675253) q[2];
sx q[2];
rz(-1.3534596) q[2];
sx q[2];
rz(-1.8390997) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.72274745) q[1];
sx q[1];
rz(-2.0993352) q[1];
sx q[1];
rz(1.8209807) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5817952) q[3];
sx q[3];
rz(-1.5526532) q[3];
sx q[3];
rz(-1.960235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.88390049) q[2];
sx q[2];
rz(-0.42282405) q[2];
sx q[2];
rz(-0.60919961) q[2];
rz(-0.60112634) q[3];
sx q[3];
rz(-1.5355891) q[3];
sx q[3];
rz(-0.48457423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6236098) q[0];
sx q[0];
rz(-1.6684775) q[0];
sx q[0];
rz(1.8286888) q[0];
rz(-3.0598705) q[1];
sx q[1];
rz(-1.7535968) q[1];
sx q[1];
rz(1.0645701) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3400187) q[0];
sx q[0];
rz(-2.193515) q[0];
sx q[0];
rz(-0.69635551) q[0];
x q[1];
rz(0.95439744) q[2];
sx q[2];
rz(-1.5570882) q[2];
sx q[2];
rz(-0.15132667) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9254522) q[1];
sx q[1];
rz(-1.6419159) q[1];
sx q[1];
rz(0.99477648) q[1];
x q[2];
rz(2.1737868) q[3];
sx q[3];
rz(-1.9120312) q[3];
sx q[3];
rz(1.6399217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1434023) q[2];
sx q[2];
rz(-0.4021796) q[2];
sx q[2];
rz(1.1581988) q[2];
rz(0.68425933) q[3];
sx q[3];
rz(-1.3254157) q[3];
sx q[3];
rz(1.6392596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96712464) q[0];
sx q[0];
rz(-3.0376349) q[0];
sx q[0];
rz(-0.52484584) q[0];
rz(1.9606494) q[1];
sx q[1];
rz(-2.3044105) q[1];
sx q[1];
rz(-2.4139074) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.380881) q[0];
sx q[0];
rz(-0.91479874) q[0];
sx q[0];
rz(-2.8674346) q[0];
x q[1];
rz(1.069017) q[2];
sx q[2];
rz(-2.2282218) q[2];
sx q[2];
rz(1.1595961) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.452949) q[1];
sx q[1];
rz(-0.51876215) q[1];
sx q[1];
rz(-2.5888836) q[1];
rz(0.91638751) q[3];
sx q[3];
rz(-1.4855412) q[3];
sx q[3];
rz(0.55201149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9119447) q[2];
sx q[2];
rz(-1.0934528) q[2];
sx q[2];
rz(1.7952807) q[2];
rz(0.66156578) q[3];
sx q[3];
rz(-1.1450359) q[3];
sx q[3];
rz(-0.0865817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0410556) q[0];
sx q[0];
rz(-2.4424398) q[0];
sx q[0];
rz(-0.93851411) q[0];
rz(-1.2988633) q[1];
sx q[1];
rz(-2.2524736) q[1];
sx q[1];
rz(3.0079957) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4105281) q[0];
sx q[0];
rz(-0.34240568) q[0];
sx q[0];
rz(-0.28241856) q[0];
rz(-pi) q[1];
rz(-0.44914896) q[2];
sx q[2];
rz(-2.5173016) q[2];
sx q[2];
rz(-0.99770944) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0444298) q[1];
sx q[1];
rz(-1.4931152) q[1];
sx q[1];
rz(-2.7809073) q[1];
rz(2.4997098) q[3];
sx q[3];
rz(-1.2716023) q[3];
sx q[3];
rz(3.0282216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9962697) q[2];
sx q[2];
rz(-1.0066373) q[2];
sx q[2];
rz(-0.88627306) q[2];
rz(2.8730872) q[3];
sx q[3];
rz(-1.0849181) q[3];
sx q[3];
rz(-1.6296384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71235424) q[0];
sx q[0];
rz(-2.7463284) q[0];
sx q[0];
rz(0.23671737) q[0];
rz(-2.9780544) q[1];
sx q[1];
rz(-1.6103368) q[1];
sx q[1];
rz(-0.18855655) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9360787) q[0];
sx q[0];
rz(-0.88919836) q[0];
sx q[0];
rz(-0.42968269) q[0];
x q[1];
rz(-2.677146) q[2];
sx q[2];
rz(-2.6147343) q[2];
sx q[2];
rz(-1.3642481) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1592387) q[1];
sx q[1];
rz(-2.6272863) q[1];
sx q[1];
rz(1.4532754) q[1];
x q[2];
rz(-2.7017835) q[3];
sx q[3];
rz(-1.4521003) q[3];
sx q[3];
rz(-2.4921093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4208372) q[2];
sx q[2];
rz(-1.6381702) q[2];
sx q[2];
rz(-2.3756964) q[2];
rz(0.01865538) q[3];
sx q[3];
rz(-1.3466287) q[3];
sx q[3];
rz(-2.6115821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50273773) q[0];
sx q[0];
rz(-1.1549594) q[0];
sx q[0];
rz(1.4422944) q[0];
rz(-2.8522708) q[1];
sx q[1];
rz(-2.0850291) q[1];
sx q[1];
rz(-2.6883968) q[1];
rz(0.53884698) q[2];
sx q[2];
rz(-1.6872023) q[2];
sx q[2];
rz(1.4414773) q[2];
rz(-0.069613386) q[3];
sx q[3];
rz(-1.1021464) q[3];
sx q[3];
rz(-2.7819421) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
