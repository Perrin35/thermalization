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
rz(-1.9967611) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24554907) q[0];
sx q[0];
rz(-2.890812) q[0];
sx q[0];
rz(1.8401237) q[0];
rz(2.8561864) q[2];
sx q[2];
rz(-1.9754344) q[2];
sx q[2];
rz(-1.0813431) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6822328) q[1];
sx q[1];
rz(-1.1182922) q[1];
sx q[1];
rz(-3.1205936) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8688698) q[3];
sx q[3];
rz(-2.094659) q[3];
sx q[3];
rz(-0.03447547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0248489) q[2];
sx q[2];
rz(-1.3329093) q[2];
sx q[2];
rz(-1.9072745) q[2];
rz(-2.0862789) q[3];
sx q[3];
rz(-1.0588667) q[3];
sx q[3];
rz(-1.9887259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9556483) q[0];
sx q[0];
rz(-1.1971104) q[0];
sx q[0];
rz(2.3117075) q[0];
rz(-2.8886967) q[1];
sx q[1];
rz(-0.77650944) q[1];
sx q[1];
rz(-0.84709644) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0704437) q[0];
sx q[0];
rz(-1.1304454) q[0];
sx q[0];
rz(-1.1308934) q[0];
rz(-2.8424938) q[2];
sx q[2];
rz(-0.80183376) q[2];
sx q[2];
rz(-0.15660827) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3165247) q[1];
sx q[1];
rz(-2.1148588) q[1];
sx q[1];
rz(-1.6506509) q[1];
rz(-pi) q[2];
rz(0.91244016) q[3];
sx q[3];
rz(-2.6077301) q[3];
sx q[3];
rz(2.562059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0236686) q[2];
sx q[2];
rz(-2.7682722) q[2];
sx q[2];
rz(-0.71933293) q[2];
rz(-1.8524648) q[3];
sx q[3];
rz(-0.23580655) q[3];
sx q[3];
rz(-1.4891362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
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
rz(2.1748623) q[1];
sx q[1];
rz(-1.2582018) q[1];
sx q[1];
rz(0.5273231) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29887154) q[0];
sx q[0];
rz(-1.5776331) q[0];
sx q[0];
rz(-0.30928916) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8965365) q[2];
sx q[2];
rz(-2.2176761) q[2];
sx q[2];
rz(1.6270454) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4434102) q[1];
sx q[1];
rz(-1.829362) q[1];
sx q[1];
rz(-2.3727388) q[1];
rz(-pi) q[2];
x q[2];
rz(0.73266352) q[3];
sx q[3];
rz(-1.3244197) q[3];
sx q[3];
rz(2.8837567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8335235) q[2];
sx q[2];
rz(-1.4493194) q[2];
sx q[2];
rz(-2.4148338) q[2];
rz(2.1440992) q[3];
sx q[3];
rz(-0.62524978) q[3];
sx q[3];
rz(1.6231026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24546394) q[0];
sx q[0];
rz(-1.6587057) q[0];
sx q[0];
rz(-1.0035275) q[0];
rz(0.040680496) q[1];
sx q[1];
rz(-1.129312) q[1];
sx q[1];
rz(0.8262659) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4429312) q[0];
sx q[0];
rz(-0.94416617) q[0];
sx q[0];
rz(-0.71039623) q[0];
rz(-pi) q[1];
rz(3.0777061) q[2];
sx q[2];
rz(-2.4404844) q[2];
sx q[2];
rz(-2.137616) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1288209) q[1];
sx q[1];
rz(-1.9073309) q[1];
sx q[1];
rz(0.080288447) q[1];
rz(-2.6310001) q[3];
sx q[3];
rz(-2.4435352) q[3];
sx q[3];
rz(-1.576168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5740009) q[2];
sx q[2];
rz(-1.8969798) q[2];
sx q[2];
rz(-0.91439247) q[2];
rz(-0.10522035) q[3];
sx q[3];
rz(-1.6243694) q[3];
sx q[3];
rz(-0.58925327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88702622) q[0];
sx q[0];
rz(-1.6396602) q[0];
sx q[0];
rz(-1.1337093) q[0];
rz(-0.0056313593) q[1];
sx q[1];
rz(-1.7040323) q[1];
sx q[1];
rz(2.5240135) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-2.4992906) q[2];
sx q[2];
rz(-1.2080492) q[2];
sx q[2];
rz(-1.158266) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4434112) q[1];
sx q[1];
rz(-2.4907506) q[1];
sx q[1];
rz(2.0425668) q[1];
rz(-pi) q[2];
x q[2];
rz(0.072237416) q[3];
sx q[3];
rz(-2.9655955) q[3];
sx q[3];
rz(-2.3776059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3474943) q[2];
sx q[2];
rz(-1.2604875) q[2];
sx q[2];
rz(0.23362544) q[2];
rz(-2.1485093) q[3];
sx q[3];
rz(-0.94998327) q[3];
sx q[3];
rz(-2.5036507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1722906) q[0];
sx q[0];
rz(-3.0811221) q[0];
sx q[0];
rz(3.0517975) q[0];
rz(-2.2604997) q[1];
sx q[1];
rz(-0.52905622) q[1];
sx q[1];
rz(0.18009137) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9371532) q[0];
sx q[0];
rz(-2.3817872) q[0];
sx q[0];
rz(2.252749) q[0];
rz(-0.18896582) q[2];
sx q[2];
rz(-2.1660888) q[2];
sx q[2];
rz(1.8408066) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0568911) q[1];
sx q[1];
rz(-1.7329721) q[1];
sx q[1];
rz(2.8497036) q[1];
x q[2];
rz(1.5363541) q[3];
sx q[3];
rz(-0.4930217) q[3];
sx q[3];
rz(2.8480414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.883541) q[2];
sx q[2];
rz(-1.5641314) q[2];
sx q[2];
rz(1.1759261) q[2];
rz(0.073143395) q[3];
sx q[3];
rz(-0.72435838) q[3];
sx q[3];
rz(-0.49889645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.183855) q[0];
sx q[0];
rz(-2.4860005) q[0];
sx q[0];
rz(-1.7234329) q[0];
rz(1.1649959) q[1];
sx q[1];
rz(-2.5823451) q[1];
sx q[1];
rz(3.1013536) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2049853) q[0];
sx q[0];
rz(-1.6523223) q[0];
sx q[0];
rz(-1.3315014) q[0];
rz(1.6778498) q[2];
sx q[2];
rz(-0.40823001) q[2];
sx q[2];
rz(-1.5232616) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.38899598) q[1];
sx q[1];
rz(-2.1850359) q[1];
sx q[1];
rz(-2.9700301) q[1];
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
rz(0.13850257) q[2];
sx q[2];
rz(-2.3562727) q[2];
sx q[2];
rz(0.17383943) q[2];
rz(-1.7447757) q[3];
sx q[3];
rz(-1.4649748) q[3];
sx q[3];
rz(2.2824536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.023507) q[0];
sx q[0];
rz(-2.9243587) q[0];
sx q[0];
rz(1.404495) q[0];
rz(-1.9346168) q[1];
sx q[1];
rz(-1.4852306) q[1];
sx q[1];
rz(1.5054024) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36639402) q[0];
sx q[0];
rz(-1.1530071) q[0];
sx q[0];
rz(1.226107) q[0];
rz(-2.5152399) q[2];
sx q[2];
rz(-2.2144631) q[2];
sx q[2];
rz(2.6477637) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7561188) q[1];
sx q[1];
rz(-1.6033152) q[1];
sx q[1];
rz(2.4454152) q[1];
rz(2.8137915) q[3];
sx q[3];
rz(-0.9947239) q[3];
sx q[3];
rz(0.93186659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.76190844) q[2];
sx q[2];
rz(-2.6613993) q[2];
sx q[2];
rz(2.9910679) q[2];
rz(-1.6020417) q[3];
sx q[3];
rz(-1.9463041) q[3];
sx q[3];
rz(0.62301821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-1.4432916) q[0];
sx q[0];
rz(-2.0069831) q[0];
sx q[0];
rz(0.64569965) q[0];
rz(0.49939108) q[1];
sx q[1];
rz(-1.7130518) q[1];
sx q[1];
rz(1.8766778) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7495959) q[0];
sx q[0];
rz(-0.8526593) q[0];
sx q[0];
rz(1.0438265) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.64847704) q[2];
sx q[2];
rz(-0.65995526) q[2];
sx q[2];
rz(-0.49396587) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9727271) q[1];
sx q[1];
rz(-1.4251514) q[1];
sx q[1];
rz(1.5345854) q[1];
x q[2];
rz(-0.42322741) q[3];
sx q[3];
rz(-1.755135) q[3];
sx q[3];
rz(1.1545899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.50679961) q[2];
sx q[2];
rz(-0.96968499) q[2];
sx q[2];
rz(0.52465049) q[2];
rz(-0.55650416) q[3];
sx q[3];
rz(-1.5126712) q[3];
sx q[3];
rz(-0.35486737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6181347) q[0];
sx q[0];
rz(-2.3045461) q[0];
sx q[0];
rz(0.24542228) q[0];
rz(-0.55631176) q[1];
sx q[1];
rz(-1.5309155) q[1];
sx q[1];
rz(1.4642749) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2357764) q[0];
sx q[0];
rz(-0.24484867) q[0];
sx q[0];
rz(0.4723627) q[0];
rz(-2.9173031) q[2];
sx q[2];
rz(-2.1456492) q[2];
sx q[2];
rz(1.526051) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.163584) q[1];
sx q[1];
rz(-0.89466909) q[1];
sx q[1];
rz(-0.34983695) q[1];
x q[2];
rz(-1.3947992) q[3];
sx q[3];
rz(-1.6547852) q[3];
sx q[3];
rz(-1.849911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2422553) q[2];
sx q[2];
rz(-1.9031886) q[2];
sx q[2];
rz(2.396092) q[2];
rz(1.4238822) q[3];
sx q[3];
rz(-0.29885492) q[3];
sx q[3];
rz(-1.1283114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2765008) q[0];
sx q[0];
rz(-1.6456974) q[0];
sx q[0];
rz(0.67281848) q[0];
rz(1.8854234) q[1];
sx q[1];
rz(-2.3354463) q[1];
sx q[1];
rz(-1.068402) q[1];
rz(-1.0829265) q[2];
sx q[2];
rz(-0.91960533) q[2];
sx q[2];
rz(-1.4056924) q[2];
rz(0.54995723) q[3];
sx q[3];
rz(-2.4933542) q[3];
sx q[3];
rz(-1.8174432) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];