OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.66981411) q[0];
sx q[0];
rz(-1.1416924) q[0];
sx q[0];
rz(-0.25622955) q[0];
rz(0.37707075) q[1];
sx q[1];
rz(-1.7989676) q[1];
sx q[1];
rz(1.0599729) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9794344) q[0];
sx q[0];
rz(-2.9720364) q[0];
sx q[0];
rz(-2.3441699) q[0];
x q[1];
rz(-0.04214123) q[2];
sx q[2];
rz(-2.4486447) q[2];
sx q[2];
rz(-2.9646089) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6457451) q[1];
sx q[1];
rz(-1.8105257) q[1];
sx q[1];
rz(1.354753) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6698631) q[3];
sx q[3];
rz(-1.4270253) q[3];
sx q[3];
rz(-3.037775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9609191) q[2];
sx q[2];
rz(-2.2283165) q[2];
sx q[2];
rz(1.712435) q[2];
rz(2.5299431) q[3];
sx q[3];
rz(-1.4432171) q[3];
sx q[3];
rz(-0.071368607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83577689) q[0];
sx q[0];
rz(-2.340305) q[0];
sx q[0];
rz(2.3961156) q[0];
rz(-1.9396793) q[1];
sx q[1];
rz(-1.8169836) q[1];
sx q[1];
rz(-1.036693) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3681524) q[0];
sx q[0];
rz(-2.7651909) q[0];
sx q[0];
rz(-1.8252116) q[0];
rz(-pi) q[1];
rz(1.5444504) q[2];
sx q[2];
rz(-2.3374632) q[2];
sx q[2];
rz(-0.94601099) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0526317) q[1];
sx q[1];
rz(-0.76238576) q[1];
sx q[1];
rz(0.41684581) q[1];
rz(-pi) q[2];
rz(-1.7880626) q[3];
sx q[3];
rz(-0.99743069) q[3];
sx q[3];
rz(-0.91106245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1408954) q[2];
sx q[2];
rz(-1.2709728) q[2];
sx q[2];
rz(-0.97266436) q[2];
rz(0.85727143) q[3];
sx q[3];
rz(-1.0664777) q[3];
sx q[3];
rz(1.8734141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5001517) q[0];
sx q[0];
rz(-2.2730136) q[0];
sx q[0];
rz(-0.80297536) q[0];
rz(0.650644) q[1];
sx q[1];
rz(-0.79901189) q[1];
sx q[1];
rz(-0.21779901) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7729491) q[0];
sx q[0];
rz(-1.9470707) q[0];
sx q[0];
rz(-0.89366389) q[0];
x q[1];
rz(2.5134263) q[2];
sx q[2];
rz(-1.3828613) q[2];
sx q[2];
rz(-0.70580259) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.21325363) q[1];
sx q[1];
rz(-1.6633185) q[1];
sx q[1];
rz(2.0615346) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59813114) q[3];
sx q[3];
rz(-2.7418828) q[3];
sx q[3];
rz(-2.3593115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.95152068) q[2];
sx q[2];
rz(-1.051544) q[2];
sx q[2];
rz(1.6311084) q[2];
rz(1.2269646) q[3];
sx q[3];
rz(-2.0881784) q[3];
sx q[3];
rz(2.1372883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25201061) q[0];
sx q[0];
rz(-0.70398206) q[0];
sx q[0];
rz(-0.66396436) q[0];
rz(-3.0162759) q[1];
sx q[1];
rz(-1.7071416) q[1];
sx q[1];
rz(2.1024316) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6923415) q[0];
sx q[0];
rz(-1.5806222) q[0];
sx q[0];
rz(1.5735606) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.98953668) q[2];
sx q[2];
rz(-0.46316812) q[2];
sx q[2];
rz(-2.4851967) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3287828) q[1];
sx q[1];
rz(-2.610958) q[1];
sx q[1];
rz(1.9932491) q[1];
rz(-pi) q[2];
rz(-1.4757206) q[3];
sx q[3];
rz(-2.0119442) q[3];
sx q[3];
rz(-0.26218647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8404954) q[2];
sx q[2];
rz(-1.7258464) q[2];
sx q[2];
rz(0.077979716) q[2];
rz(-1.1178499) q[3];
sx q[3];
rz(-2.100914) q[3];
sx q[3];
rz(2.1054721) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9341105) q[0];
sx q[0];
rz(-0.20620646) q[0];
sx q[0];
rz(0.41314405) q[0];
rz(-1.235599) q[1];
sx q[1];
rz(-1.3256336) q[1];
sx q[1];
rz(1.8280169) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8395555) q[0];
sx q[0];
rz(-1.3875204) q[0];
sx q[0];
rz(0.9367783) q[0];
rz(-pi) q[1];
rz(-0.92720534) q[2];
sx q[2];
rz(-2.5816133) q[2];
sx q[2];
rz(-0.72712979) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6690065) q[1];
sx q[1];
rz(-2.5372977) q[1];
sx q[1];
rz(-0.86493203) q[1];
rz(-2.026346) q[3];
sx q[3];
rz(-0.83720647) q[3];
sx q[3];
rz(2.7546492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1818587) q[2];
sx q[2];
rz(-1.2343957) q[2];
sx q[2];
rz(0.44446298) q[2];
rz(-1.9410761) q[3];
sx q[3];
rz(-1.8903172) q[3];
sx q[3];
rz(3.0357231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34894094) q[0];
sx q[0];
rz(-2.7466725) q[0];
sx q[0];
rz(-0.077202395) q[0];
rz(0.96241799) q[1];
sx q[1];
rz(-1.187477) q[1];
sx q[1];
rz(0.033800689) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1241959) q[0];
sx q[0];
rz(-2.5585045) q[0];
sx q[0];
rz(-0.65753196) q[0];
rz(-pi) q[1];
rz(0.54767139) q[2];
sx q[2];
rz(-1.593809) q[2];
sx q[2];
rz(-0.097824899) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5742366) q[1];
sx q[1];
rz(-2.1891382) q[1];
sx q[1];
rz(2.7537557) q[1];
rz(-2.364879) q[3];
sx q[3];
rz(-1.6459822) q[3];
sx q[3];
rz(3.1326339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.17322156) q[2];
sx q[2];
rz(-1.6049478) q[2];
sx q[2];
rz(2.8823493) q[2];
rz(-2.1974468) q[3];
sx q[3];
rz(-0.21374948) q[3];
sx q[3];
rz(-1.3313782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060870085) q[0];
sx q[0];
rz(-0.20651564) q[0];
sx q[0];
rz(2.5282705) q[0];
rz(0.90336409) q[1];
sx q[1];
rz(-0.84066835) q[1];
sx q[1];
rz(0.14796743) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67576236) q[0];
sx q[0];
rz(-1.0935385) q[0];
sx q[0];
rz(0.81264074) q[0];
x q[1];
rz(1.5388707) q[2];
sx q[2];
rz(-2.0096001) q[2];
sx q[2];
rz(1.9500458) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2206055) q[1];
sx q[1];
rz(-1.9927532) q[1];
sx q[1];
rz(0.76442952) q[1];
x q[2];
rz(-1.5336125) q[3];
sx q[3];
rz(-2.5761009) q[3];
sx q[3];
rz(-2.9786106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.938574) q[2];
sx q[2];
rz(-1.729894) q[2];
sx q[2];
rz(2.1417248) q[2];
rz(-1.3153007) q[3];
sx q[3];
rz(-0.2838997) q[3];
sx q[3];
rz(-0.6704754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1891747) q[0];
sx q[0];
rz(-2.2947831) q[0];
sx q[0];
rz(2.1536105) q[0];
rz(1.2692163) q[1];
sx q[1];
rz(-0.80943426) q[1];
sx q[1];
rz(-0.16407897) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1647427) q[0];
sx q[0];
rz(-1.4535507) q[0];
sx q[0];
rz(1.6567203) q[0];
x q[1];
rz(2.9226926) q[2];
sx q[2];
rz(-2.0563807) q[2];
sx q[2];
rz(1.4771493) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8089701) q[1];
sx q[1];
rz(-1.5738942) q[1];
sx q[1];
rz(1.078152) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3173728) q[3];
sx q[3];
rz(-2.4755423) q[3];
sx q[3];
rz(2.1061153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8998731) q[2];
sx q[2];
rz(-0.49979979) q[2];
sx q[2];
rz(-1.0618173) q[2];
rz(0.1344943) q[3];
sx q[3];
rz(-2.2420292) q[3];
sx q[3];
rz(3.0883446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6580842) q[0];
sx q[0];
rz(-0.71165076) q[0];
sx q[0];
rz(-0.2051951) q[0];
rz(2.9070053) q[1];
sx q[1];
rz(-1.4261475) q[1];
sx q[1];
rz(-0.1964143) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.043606) q[0];
sx q[0];
rz(-0.65904407) q[0];
sx q[0];
rz(2.090467) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8750849) q[2];
sx q[2];
rz(-0.54000137) q[2];
sx q[2];
rz(-2.2842479) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3383693) q[1];
sx q[1];
rz(-1.3937989) q[1];
sx q[1];
rz(1.8990252) q[1];
x q[2];
rz(0.78077448) q[3];
sx q[3];
rz(-1.9650243) q[3];
sx q[3];
rz(-3.0633283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3966169) q[2];
sx q[2];
rz(-1.1810415) q[2];
sx q[2];
rz(-1.9632001) q[2];
rz(-2.7922503) q[3];
sx q[3];
rz(-2.0125466) q[3];
sx q[3];
rz(-1.0497302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1402533) q[0];
sx q[0];
rz(-2.0297191) q[0];
sx q[0];
rz(-2.4358791) q[0];
rz(-2.4440675) q[1];
sx q[1];
rz(-1.3730647) q[1];
sx q[1];
rz(-0.90550214) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5030497) q[0];
sx q[0];
rz(-0.43404225) q[0];
sx q[0];
rz(0.040157138) q[0];
rz(-pi) q[1];
rz(-2.7780813) q[2];
sx q[2];
rz(-1.264252) q[2];
sx q[2];
rz(2.6078122) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9230629) q[1];
sx q[1];
rz(-2.8273812) q[1];
sx q[1];
rz(2.988222) q[1];
x q[2];
rz(-0.8373413) q[3];
sx q[3];
rz(-0.47347927) q[3];
sx q[3];
rz(-0.66094962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8991578) q[2];
sx q[2];
rz(-1.0589212) q[2];
sx q[2];
rz(2.4230797) q[2];
rz(2.4588623) q[3];
sx q[3];
rz(-1.1444789) q[3];
sx q[3];
rz(-0.1040641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.318442) q[0];
sx q[0];
rz(-2.5704076) q[0];
sx q[0];
rz(3.0388863) q[0];
rz(1.5009343) q[1];
sx q[1];
rz(-1.2702912) q[1];
sx q[1];
rz(2.695695) q[1];
rz(-0.6871625) q[2];
sx q[2];
rz(-2.9570815) q[2];
sx q[2];
rz(2.2312192) q[2];
rz(0.64651505) q[3];
sx q[3];
rz(-0.91520354) q[3];
sx q[3];
rz(-1.147816) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
