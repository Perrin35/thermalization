OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.174515) q[0];
sx q[0];
rz(-1.6835901) q[0];
sx q[0];
rz(-0.30012497) q[0];
rz(1.1974273) q[1];
sx q[1];
rz(-1.6150885) q[1];
sx q[1];
rz(1.6405029) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7952657) q[0];
sx q[0];
rz(-0.015595868) q[0];
sx q[0];
rz(-2.8898284) q[0];
x q[1];
rz(-1.9674703) q[2];
sx q[2];
rz(-0.80840092) q[2];
sx q[2];
rz(1.8044745) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.056880205) q[1];
sx q[1];
rz(-2.4507629) q[1];
sx q[1];
rz(0.34617306) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5759104) q[3];
sx q[3];
rz(-2.3795914) q[3];
sx q[3];
rz(0.68540547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7370721) q[2];
sx q[2];
rz(-1.1417192) q[2];
sx q[2];
rz(-0.70297757) q[2];
rz(-0.56973714) q[3];
sx q[3];
rz(-0.8786141) q[3];
sx q[3];
rz(1.5841293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0074145929) q[0];
sx q[0];
rz(-0.61892048) q[0];
sx q[0];
rz(-1.9084357) q[0];
rz(0.59457072) q[1];
sx q[1];
rz(-2.1023127) q[1];
sx q[1];
rz(2.0436683) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9482518) q[0];
sx q[0];
rz(-1.3701212) q[0];
sx q[0];
rz(-3.0737123) q[0];
rz(0.64867257) q[2];
sx q[2];
rz(-2.5498516) q[2];
sx q[2];
rz(1.5209188) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.532053) q[1];
sx q[1];
rz(-1.5455478) q[1];
sx q[1];
rz(1.6482501) q[1];
rz(-1.6637832) q[3];
sx q[3];
rz(-0.78130296) q[3];
sx q[3];
rz(-2.6546991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7132831) q[2];
sx q[2];
rz(-1.2998591) q[2];
sx q[2];
rz(-1.5353047) q[2];
rz(-0.71896583) q[3];
sx q[3];
rz(-2.1956367) q[3];
sx q[3];
rz(-2.9276221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3035901) q[0];
sx q[0];
rz(-0.8135697) q[0];
sx q[0];
rz(-0.59233061) q[0];
rz(-1.8968286) q[1];
sx q[1];
rz(-2.0015621) q[1];
sx q[1];
rz(2.9959784) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1942822) q[0];
sx q[0];
rz(-2.6366451) q[0];
sx q[0];
rz(-2.8584252) q[0];
rz(0.12965173) q[2];
sx q[2];
rz(-2.2551422) q[2];
sx q[2];
rz(0.42829542) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9551202) q[1];
sx q[1];
rz(-1.8545517) q[1];
sx q[1];
rz(-2.1452745) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0907902) q[3];
sx q[3];
rz(-2.601805) q[3];
sx q[3];
rz(1.2320428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.738872) q[2];
sx q[2];
rz(-1.0307743) q[2];
sx q[2];
rz(-1.9090451) q[2];
rz(-2.9018719) q[3];
sx q[3];
rz(-1.0545878) q[3];
sx q[3];
rz(0.48503748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1211014) q[0];
sx q[0];
rz(-0.70206577) q[0];
sx q[0];
rz(-3.1109911) q[0];
rz(2.5864511) q[1];
sx q[1];
rz(-1.6629013) q[1];
sx q[1];
rz(1.0928924) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26390938) q[0];
sx q[0];
rz(-1.0315686) q[0];
sx q[0];
rz(-0.79099057) q[0];
rz(0.76538122) q[2];
sx q[2];
rz(-1.5281406) q[2];
sx q[2];
rz(2.9431107) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79697414) q[1];
sx q[1];
rz(-0.94616468) q[1];
sx q[1];
rz(-0.25298869) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3044358) q[3];
sx q[3];
rz(-0.87292307) q[3];
sx q[3];
rz(2.9540889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0838202) q[2];
sx q[2];
rz(-1.3670992) q[2];
sx q[2];
rz(2.8687381) q[2];
rz(1.7504292) q[3];
sx q[3];
rz(-0.83943668) q[3];
sx q[3];
rz(2.1896037) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6735753) q[0];
sx q[0];
rz(-1.3382358) q[0];
sx q[0];
rz(-1.2982298) q[0];
rz(-2.4507554) q[1];
sx q[1];
rz(-1.9419443) q[1];
sx q[1];
rz(1.615049) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0218791) q[0];
sx q[0];
rz(-2.102431) q[0];
sx q[0];
rz(-2.3945532) q[0];
x q[1];
rz(-1.853302) q[2];
sx q[2];
rz(-1.2307248) q[2];
sx q[2];
rz(-1.9486537) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1880875) q[1];
sx q[1];
rz(-1.5966478) q[1];
sx q[1];
rz(0.8698747) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2722716) q[3];
sx q[3];
rz(-1.9179452) q[3];
sx q[3];
rz(2.783028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4096628) q[2];
sx q[2];
rz(-0.80594984) q[2];
sx q[2];
rz(2.9742677) q[2];
rz(-1.7512199) q[3];
sx q[3];
rz(-2.6087587) q[3];
sx q[3];
rz(-0.65756857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.9996027) q[0];
sx q[0];
rz(-0.36407343) q[0];
sx q[0];
rz(1.863119) q[0];
rz(-1.7059884) q[1];
sx q[1];
rz(-1.4878788) q[1];
sx q[1];
rz(2.7659168) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6784793) q[0];
sx q[0];
rz(-2.7857384) q[0];
sx q[0];
rz(-2.2694017) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7406875) q[2];
sx q[2];
rz(-1.4235911) q[2];
sx q[2];
rz(2.3224555) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.76498079) q[1];
sx q[1];
rz(-1.7936447) q[1];
sx q[1];
rz(1.8972023) q[1];
x q[2];
rz(-2.385684) q[3];
sx q[3];
rz(-1.7497853) q[3];
sx q[3];
rz(0.32839113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0756691) q[2];
sx q[2];
rz(-1.064294) q[2];
sx q[2];
rz(-0.86090487) q[2];
rz(1.9979477) q[3];
sx q[3];
rz(-2.836561) q[3];
sx q[3];
rz(2.8633269) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.025573108) q[0];
sx q[0];
rz(-1.8208068) q[0];
sx q[0];
rz(-2.386911) q[0];
rz(0.93983752) q[1];
sx q[1];
rz(-0.87516963) q[1];
sx q[1];
rz(-1.2831203) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1719753) q[0];
sx q[0];
rz(-1.8775058) q[0];
sx q[0];
rz(-3.1344916) q[0];
rz(-2.0286328) q[2];
sx q[2];
rz(-1.2414059) q[2];
sx q[2];
rz(0.068458099) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2312113) q[1];
sx q[1];
rz(-1.2802135) q[1];
sx q[1];
rz(-0.39931764) q[1];
rz(-0.51605172) q[3];
sx q[3];
rz(-2.1911804) q[3];
sx q[3];
rz(-0.32143264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5963886) q[2];
sx q[2];
rz(-1.117492) q[2];
sx q[2];
rz(-2.0580573) q[2];
rz(2.3986473) q[3];
sx q[3];
rz(-1.5321621) q[3];
sx q[3];
rz(1.4325745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1525986) q[0];
sx q[0];
rz(-0.78482634) q[0];
sx q[0];
rz(-0.75491828) q[0];
rz(-0.49916357) q[1];
sx q[1];
rz(-2.1776336) q[1];
sx q[1];
rz(2.4430433) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38235086) q[0];
sx q[0];
rz(-0.97433358) q[0];
sx q[0];
rz(-1.9475219) q[0];
rz(2.8590747) q[2];
sx q[2];
rz(-2.660224) q[2];
sx q[2];
rz(0.11061874) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9883635) q[1];
sx q[1];
rz(-2.6090342) q[1];
sx q[1];
rz(1.531324) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9898765) q[3];
sx q[3];
rz(-1.5718096) q[3];
sx q[3];
rz(0.59194293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2373206) q[2];
sx q[2];
rz(-2.1492683) q[2];
sx q[2];
rz(2.3392056) q[2];
rz(2.5896416) q[3];
sx q[3];
rz(-0.80058432) q[3];
sx q[3];
rz(-0.62057453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9629843) q[0];
sx q[0];
rz(-1.4405788) q[0];
sx q[0];
rz(-2.0340023) q[0];
rz(1.3811318) q[1];
sx q[1];
rz(-1.3637204) q[1];
sx q[1];
rz(1.6345056) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0213623) q[0];
sx q[0];
rz(-1.7108546) q[0];
sx q[0];
rz(-0.39744486) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31839026) q[2];
sx q[2];
rz(-1.9620775) q[2];
sx q[2];
rz(-1.1423781) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9642051) q[1];
sx q[1];
rz(-1.2086165) q[1];
sx q[1];
rz(-1.7262579) q[1];
x q[2];
rz(-0.82689311) q[3];
sx q[3];
rz(-1.6312459) q[3];
sx q[3];
rz(-1.7884487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0570809) q[2];
sx q[2];
rz(-0.19597404) q[2];
sx q[2];
rz(-1.8692807) q[2];
rz(1.48014) q[3];
sx q[3];
rz(-2.0062165) q[3];
sx q[3];
rz(-2.541466) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24499527) q[0];
sx q[0];
rz(-0.56140459) q[0];
sx q[0];
rz(-1.2835314) q[0];
rz(3.1237579) q[1];
sx q[1];
rz(-2.5051038) q[1];
sx q[1];
rz(3.0160115) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3888729) q[0];
sx q[0];
rz(-2.5872091) q[0];
sx q[0];
rz(1.5977018) q[0];
rz(-pi) q[1];
rz(2.0415885) q[2];
sx q[2];
rz(-2.1934237) q[2];
sx q[2];
rz(2.560844) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8931742) q[1];
sx q[1];
rz(-1.9918973) q[1];
sx q[1];
rz(2.0401272) q[1];
x q[2];
rz(2.2108881) q[3];
sx q[3];
rz(-0.32880201) q[3];
sx q[3];
rz(0.14842515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.77832001) q[2];
sx q[2];
rz(-1.2691701) q[2];
sx q[2];
rz(-0.8030836) q[2];
rz(2.965029) q[3];
sx q[3];
rz(-1.8162138) q[3];
sx q[3];
rz(-0.77809063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9618027) q[0];
sx q[0];
rz(-2.640124) q[0];
sx q[0];
rz(-1.5928706) q[0];
rz(1.8190307) q[1];
sx q[1];
rz(-1.7772728) q[1];
sx q[1];
rz(-1.5900236) q[1];
rz(1.6336468) q[2];
sx q[2];
rz(-1.3286776) q[2];
sx q[2];
rz(1.5059289) q[2];
rz(2.8040721) q[3];
sx q[3];
rz(-1.6277085) q[3];
sx q[3];
rz(1.7276702) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
