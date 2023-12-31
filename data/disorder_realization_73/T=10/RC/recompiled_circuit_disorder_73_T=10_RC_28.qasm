OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.9956545) q[0];
sx q[0];
rz(-0.50322682) q[0];
sx q[0];
rz(2.4174262) q[0];
rz(6.9231482) q[1];
sx q[1];
rz(5.7531113) q[1];
sx q[1];
rz(2.35676) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9871702) q[0];
sx q[0];
rz(-1.3599456) q[0];
sx q[0];
rz(-0.2984557) q[0];
rz(-0.41854026) q[2];
sx q[2];
rz(-1.4782018) q[2];
sx q[2];
rz(-1.4946403) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6341056) q[1];
sx q[1];
rz(-1.3480942) q[1];
sx q[1];
rz(1.0313862) q[1];
rz(-pi) q[2];
x q[2];
rz(1.647244) q[3];
sx q[3];
rz(-0.41562286) q[3];
sx q[3];
rz(-1.7474183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.589754) q[2];
sx q[2];
rz(-1.7244312) q[2];
sx q[2];
rz(0.067967728) q[2];
rz(-3.0170278) q[3];
sx q[3];
rz(-0.3228651) q[3];
sx q[3];
rz(1.386806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92000604) q[0];
sx q[0];
rz(-3.0060372) q[0];
sx q[0];
rz(-2.8979229) q[0];
rz(-2.5098353) q[1];
sx q[1];
rz(-1.4032204) q[1];
sx q[1];
rz(-1.3557281) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5403554) q[0];
sx q[0];
rz(-1.8260801) q[0];
sx q[0];
rz(0.028098696) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35661125) q[2];
sx q[2];
rz(-1.9372254) q[2];
sx q[2];
rz(-0.29417843) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1214952) q[1];
sx q[1];
rz(-1.1250245) q[1];
sx q[1];
rz(2.9989468) q[1];
rz(2.5494266) q[3];
sx q[3];
rz(-1.6093996) q[3];
sx q[3];
rz(-1.6174699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0791066) q[2];
sx q[2];
rz(-2.1495543) q[2];
sx q[2];
rz(2.8919354) q[2];
rz(0.50659531) q[3];
sx q[3];
rz(-1.6258312) q[3];
sx q[3];
rz(0.33199582) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8963985) q[0];
sx q[0];
rz(-1.2250552) q[0];
sx q[0];
rz(0.89843345) q[0];
rz(1.8067182) q[1];
sx q[1];
rz(-1.2355665) q[1];
sx q[1];
rz(1.2737087) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0262336) q[0];
sx q[0];
rz(-0.44239487) q[0];
sx q[0];
rz(0.94985234) q[0];
rz(-pi) q[1];
rz(-2.8286335) q[2];
sx q[2];
rz(-2.8319781) q[2];
sx q[2];
rz(1.5086053) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.102036) q[1];
sx q[1];
rz(-1.2219056) q[1];
sx q[1];
rz(-2.4400913) q[1];
rz(0.55176118) q[3];
sx q[3];
rz(-0.80358395) q[3];
sx q[3];
rz(1.4078275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0597824) q[2];
sx q[2];
rz(-2.5368097) q[2];
sx q[2];
rz(0.95345062) q[2];
rz(3.1070784) q[3];
sx q[3];
rz(-2.3551066) q[3];
sx q[3];
rz(0.22687337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-1.2801441) q[0];
sx q[0];
rz(-2.9252958) q[0];
sx q[0];
rz(-2.8934073) q[0];
rz(1.0379627) q[1];
sx q[1];
rz(-2.018441) q[1];
sx q[1];
rz(-0.074137069) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.603133) q[0];
sx q[0];
rz(-0.54766253) q[0];
sx q[0];
rz(1.0663701) q[0];
rz(-1.147869) q[2];
sx q[2];
rz(-0.72407702) q[2];
sx q[2];
rz(-0.99087447) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.292865) q[1];
sx q[1];
rz(-1.8855727) q[1];
sx q[1];
rz(-3.0175812) q[1];
rz(-pi) q[2];
rz(-1.7449964) q[3];
sx q[3];
rz(-1.3341691) q[3];
sx q[3];
rz(-0.06121204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2512102) q[2];
sx q[2];
rz(-2.7373098) q[2];
sx q[2];
rz(0.038643535) q[2];
rz(-0.97366992) q[3];
sx q[3];
rz(-0.49574167) q[3];
sx q[3];
rz(-0.27004778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6073109) q[0];
sx q[0];
rz(-1.6058291) q[0];
sx q[0];
rz(-1.3624396) q[0];
rz(-0.81659395) q[1];
sx q[1];
rz(-1.2885619) q[1];
sx q[1];
rz(1.978925) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029318854) q[0];
sx q[0];
rz(-3.0262039) q[0];
sx q[0];
rz(-1.3089887) q[0];
x q[1];
rz(0.16935279) q[2];
sx q[2];
rz(-1.2522109) q[2];
sx q[2];
rz(-2.7404075) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77046493) q[1];
sx q[1];
rz(-2.7264997) q[1];
sx q[1];
rz(-2.0626555) q[1];
rz(2.2546222) q[3];
sx q[3];
rz(-0.63347048) q[3];
sx q[3];
rz(-0.57852832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6267307) q[2];
sx q[2];
rz(-2.0613487) q[2];
sx q[2];
rz(-2.999372) q[2];
rz(0.90406117) q[3];
sx q[3];
rz(-1.3198493) q[3];
sx q[3];
rz(2.952125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11480039) q[0];
sx q[0];
rz(-2.8650706) q[0];
sx q[0];
rz(-1.4676771) q[0];
rz(2.5698075) q[1];
sx q[1];
rz(-0.3586868) q[1];
sx q[1];
rz(-0.30803672) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15333262) q[0];
sx q[0];
rz(-1.6969661) q[0];
sx q[0];
rz(-1.4754962) q[0];
x q[1];
rz(-2.3601989) q[2];
sx q[2];
rz(-0.77971824) q[2];
sx q[2];
rz(-1.7674854) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8923556) q[1];
sx q[1];
rz(-2.4294469) q[1];
sx q[1];
rz(-1.2675136) q[1];
x q[2];
rz(-1.4792535) q[3];
sx q[3];
rz(-1.9178101) q[3];
sx q[3];
rz(0.43859827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3623111) q[2];
sx q[2];
rz(-1.7411391) q[2];
sx q[2];
rz(1.1479088) q[2];
rz(2.4273196) q[3];
sx q[3];
rz(-0.89759421) q[3];
sx q[3];
rz(-0.64546293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4797453) q[0];
sx q[0];
rz(-2.3006738) q[0];
sx q[0];
rz(-3.0116144) q[0];
rz(-0.030844363) q[1];
sx q[1];
rz(-1.8519311) q[1];
sx q[1];
rz(-0.67108363) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0111611) q[0];
sx q[0];
rz(-1.5357619) q[0];
sx q[0];
rz(3.0879211) q[0];
rz(-pi) q[1];
rz(-0.030521557) q[2];
sx q[2];
rz(-1.754921) q[2];
sx q[2];
rz(0.39779278) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.56747251) q[1];
sx q[1];
rz(-1.6735055) q[1];
sx q[1];
rz(1.9436388) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89550771) q[3];
sx q[3];
rz(-1.607778) q[3];
sx q[3];
rz(-1.0469588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.122763) q[2];
sx q[2];
rz(-1.6100223) q[2];
sx q[2];
rz(-0.57787952) q[2];
rz(-3.1130062) q[3];
sx q[3];
rz(-1.8609906) q[3];
sx q[3];
rz(-1.2602497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9776483) q[0];
sx q[0];
rz(-0.90286911) q[0];
sx q[0];
rz(-0.41241616) q[0];
rz(-1.4498129) q[1];
sx q[1];
rz(-1.342536) q[1];
sx q[1];
rz(1.1669881) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1056571) q[0];
sx q[0];
rz(-1.3577537) q[0];
sx q[0];
rz(2.1980594) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27600482) q[2];
sx q[2];
rz(-0.88062421) q[2];
sx q[2];
rz(-1.3921757) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.2003157) q[1];
sx q[1];
rz(-1.3822767) q[1];
sx q[1];
rz(1.6016866) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8230209) q[3];
sx q[3];
rz(-0.54402292) q[3];
sx q[3];
rz(-2.9094484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2010487) q[2];
sx q[2];
rz(-2.2622435) q[2];
sx q[2];
rz(2.9525625) q[2];
rz(-0.14686251) q[3];
sx q[3];
rz(-0.18460128) q[3];
sx q[3];
rz(-1.3930901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1259595) q[0];
sx q[0];
rz(-1.8305612) q[0];
sx q[0];
rz(-2.2145859) q[0];
rz(-1.758763) q[1];
sx q[1];
rz(-0.60950509) q[1];
sx q[1];
rz(1.6519201) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.888436) q[0];
sx q[0];
rz(-0.90421593) q[0];
sx q[0];
rz(2.1783834) q[0];
rz(-3.0909806) q[2];
sx q[2];
rz(-1.0992556) q[2];
sx q[2];
rz(2.5052349) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.97396321) q[1];
sx q[1];
rz(-1.7427674) q[1];
sx q[1];
rz(-2.9457438) q[1];
rz(-pi) q[2];
rz(2.3143523) q[3];
sx q[3];
rz(-0.67875553) q[3];
sx q[3];
rz(2.7620706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5513409) q[2];
sx q[2];
rz(-2.6987023) q[2];
sx q[2];
rz(-1.7112973) q[2];
rz(2.5643505) q[3];
sx q[3];
rz(-0.8876628) q[3];
sx q[3];
rz(1.3841217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5532613) q[0];
sx q[0];
rz(-1.7734779) q[0];
sx q[0];
rz(-0.28840315) q[0];
rz(-0.53238955) q[1];
sx q[1];
rz(-0.45982292) q[1];
sx q[1];
rz(-0.14702252) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4924016) q[0];
sx q[0];
rz(-0.99897879) q[0];
sx q[0];
rz(-2.6834821) q[0];
rz(-pi) q[1];
rz(2.8617919) q[2];
sx q[2];
rz(-0.24926148) q[2];
sx q[2];
rz(1.793043) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.22132561) q[1];
sx q[1];
rz(-2.1888071) q[1];
sx q[1];
rz(3.0926535) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23707323) q[3];
sx q[3];
rz(-2.5519538) q[3];
sx q[3];
rz(-1.3414563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0599351) q[2];
sx q[2];
rz(-2.7066878) q[2];
sx q[2];
rz(0.74404136) q[2];
rz(-0.75731164) q[3];
sx q[3];
rz(-1.7777187) q[3];
sx q[3];
rz(1.3967167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.025678) q[0];
sx q[0];
rz(-1.0703351) q[0];
sx q[0];
rz(-1.0967789) q[0];
rz(2.3241282) q[1];
sx q[1];
rz(-1.2066963) q[1];
sx q[1];
rz(-0.6304601) q[1];
rz(1.5031917) q[2];
sx q[2];
rz(-2.0220145) q[2];
sx q[2];
rz(1.8215712) q[2];
rz(-2.1327303) q[3];
sx q[3];
rz(-1.6812134) q[3];
sx q[3];
rz(-2.5440661) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
