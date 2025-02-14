OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.661854) q[0];
sx q[0];
rz(5.148611) q[0];
sx q[0];
rz(10.892286) q[0];
rz(1.6788586) q[1];
sx q[1];
rz(-2.1958513) q[1];
sx q[1];
rz(0.50409627) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.534908) q[0];
sx q[0];
rz(-2.1143171) q[0];
sx q[0];
rz(-2.1786819) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7376704) q[2];
sx q[2];
rz(-1.9144399) q[2];
sx q[2];
rz(2.0640896) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2368184) q[1];
sx q[1];
rz(-0.95413757) q[1];
sx q[1];
rz(0.84994933) q[1];
rz(-pi) q[2];
rz(-2.8064742) q[3];
sx q[3];
rz(-1.6168645) q[3];
sx q[3];
rz(-0.81460458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.41245875) q[2];
sx q[2];
rz(-1.2428186) q[2];
sx q[2];
rz(-0.18623713) q[2];
rz(1.3650182) q[3];
sx q[3];
rz(-2.5297207) q[3];
sx q[3];
rz(1.2046825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3314826) q[0];
sx q[0];
rz(-0.45833603) q[0];
sx q[0];
rz(-3.0895184) q[0];
rz(2.1049818) q[1];
sx q[1];
rz(-0.60984817) q[1];
sx q[1];
rz(-2.5263272) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0947845) q[0];
sx q[0];
rz(-0.9082709) q[0];
sx q[0];
rz(-2.725968) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7024666) q[2];
sx q[2];
rz(-1.8063917) q[2];
sx q[2];
rz(2.5680096) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4897291) q[1];
sx q[1];
rz(-2.4861838) q[1];
sx q[1];
rz(-1.4784527) q[1];
rz(-pi) q[2];
rz(-1.3619945) q[3];
sx q[3];
rz(-2.3376645) q[3];
sx q[3];
rz(-1.3458136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.44876862) q[2];
sx q[2];
rz(-1.245446) q[2];
sx q[2];
rz(1.2358865) q[2];
rz(-2.0125194) q[3];
sx q[3];
rz(-2.2344507) q[3];
sx q[3];
rz(-0.83704078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6579984) q[0];
sx q[0];
rz(-0.75355419) q[0];
sx q[0];
rz(1.3760706) q[0];
rz(0.54028571) q[1];
sx q[1];
rz(-2.1018201) q[1];
sx q[1];
rz(1.4556063) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7175563) q[0];
sx q[0];
rz(-2.7323595) q[0];
sx q[0];
rz(0.17828973) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4810782) q[2];
sx q[2];
rz(-1.3680581) q[2];
sx q[2];
rz(2.3281472) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.27764717) q[1];
sx q[1];
rz(-1.9233812) q[1];
sx q[1];
rz(-0.90971281) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0946526) q[3];
sx q[3];
rz(-2.3848001) q[3];
sx q[3];
rz(1.688129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6471863) q[2];
sx q[2];
rz(-2.8072085) q[2];
sx q[2];
rz(-1.6953267) q[2];
rz(-1.9485731) q[3];
sx q[3];
rz(-2.1665067) q[3];
sx q[3];
rz(-1.4421991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0874264) q[0];
sx q[0];
rz(-2.8717201) q[0];
sx q[0];
rz(-0.42359459) q[0];
rz(-2.1845747) q[1];
sx q[1];
rz(-1.3514163) q[1];
sx q[1];
rz(-0.097600309) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97948697) q[0];
sx q[0];
rz(-0.66363664) q[0];
sx q[0];
rz(0.62941334) q[0];
x q[1];
rz(-0.1431485) q[2];
sx q[2];
rz(-2.2586056) q[2];
sx q[2];
rz(2.0542415) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7047119) q[1];
sx q[1];
rz(-1.2578674) q[1];
sx q[1];
rz(1.9347316) q[1];
rz(-pi) q[2];
x q[2];
rz(2.76582) q[3];
sx q[3];
rz(-0.76533459) q[3];
sx q[3];
rz(-0.96804141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.61509722) q[2];
sx q[2];
rz(-0.47250938) q[2];
sx q[2];
rz(0.16461593) q[2];
rz(-0.047867157) q[3];
sx q[3];
rz(-1.7367626) q[3];
sx q[3];
rz(-0.78711787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1614302) q[0];
sx q[0];
rz(-0.4466559) q[0];
sx q[0];
rz(1.584378) q[0];
rz(-0.36852512) q[1];
sx q[1];
rz(-1.3216113) q[1];
sx q[1];
rz(1.535086) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42035746) q[0];
sx q[0];
rz(-2.0518655) q[0];
sx q[0];
rz(0.98929437) q[0];
x q[1];
rz(-0.82442762) q[2];
sx q[2];
rz(-1.2807089) q[2];
sx q[2];
rz(0.96765358) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8884224) q[1];
sx q[1];
rz(-2.0520323) q[1];
sx q[1];
rz(2.8594244) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8312884) q[3];
sx q[3];
rz(-0.29803571) q[3];
sx q[3];
rz(-0.71825114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2313472) q[2];
sx q[2];
rz(-0.37800899) q[2];
sx q[2];
rz(-2.0965915) q[2];
rz(-0.16799489) q[3];
sx q[3];
rz(-1.955227) q[3];
sx q[3];
rz(2.7872938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1217693) q[0];
sx q[0];
rz(-1.0291809) q[0];
sx q[0];
rz(-1.9371012) q[0];
rz(2.3380741) q[1];
sx q[1];
rz(-2.3280227) q[1];
sx q[1];
rz(-2.4258851) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89838081) q[0];
sx q[0];
rz(-1.0168494) q[0];
sx q[0];
rz(0.57147632) q[0];
x q[1];
rz(-2.6936221) q[2];
sx q[2];
rz(-1.7501737) q[2];
sx q[2];
rz(2.568754) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0790076) q[1];
sx q[1];
rz(-0.91168303) q[1];
sx q[1];
rz(1.4277894) q[1];
x q[2];
rz(2.8479667) q[3];
sx q[3];
rz(-2.01247) q[3];
sx q[3];
rz(-2.3229118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7314926) q[2];
sx q[2];
rz(-2.4918753) q[2];
sx q[2];
rz(-1.043172) q[2];
rz(0.10797524) q[3];
sx q[3];
rz(-2.2004746) q[3];
sx q[3];
rz(-2.030355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1825948) q[0];
sx q[0];
rz(-1.3866871) q[0];
sx q[0];
rz(1.9166272) q[0];
rz(0.23172465) q[1];
sx q[1];
rz(-2.2544315) q[1];
sx q[1];
rz(-1.8909594) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5179493) q[0];
sx q[0];
rz(-0.8276437) q[0];
sx q[0];
rz(-1.736899) q[0];
rz(1.6804303) q[2];
sx q[2];
rz(-0.19936745) q[2];
sx q[2];
rz(-0.91857544) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2030484) q[1];
sx q[1];
rz(-0.37751679) q[1];
sx q[1];
rz(2.50752) q[1];
rz(1.5579079) q[3];
sx q[3];
rz(-2.1272214) q[3];
sx q[3];
rz(-0.22695676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.072824868) q[2];
sx q[2];
rz(-2.4710957) q[2];
sx q[2];
rz(-0.13216275) q[2];
rz(-0.69563785) q[3];
sx q[3];
rz(-1.1824824) q[3];
sx q[3];
rz(-2.3343991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.950133) q[0];
sx q[0];
rz(-1.5605254) q[0];
sx q[0];
rz(0.70924846) q[0];
rz(-1.5059772) q[1];
sx q[1];
rz(-1.154107) q[1];
sx q[1];
rz(2.0567315) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70299545) q[0];
sx q[0];
rz(-0.60351047) q[0];
sx q[0];
rz(0.10137697) q[0];
rz(2.3368767) q[2];
sx q[2];
rz(-1.448481) q[2];
sx q[2];
rz(1.027077) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0811396) q[1];
sx q[1];
rz(-2.6067002) q[1];
sx q[1];
rz(-0.51523955) q[1];
x q[2];
rz(-1.3832757) q[3];
sx q[3];
rz(-3.0802343) q[3];
sx q[3];
rz(0.96422577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6214577) q[2];
sx q[2];
rz(-2.1389213) q[2];
sx q[2];
rz(-1.8247068) q[2];
rz(-0.78553158) q[3];
sx q[3];
rz(-1.9436049) q[3];
sx q[3];
rz(2.7151916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21996466) q[0];
sx q[0];
rz(-1.6554609) q[0];
sx q[0];
rz(-0.40837902) q[0];
rz(-2.1829103) q[1];
sx q[1];
rz(-0.32902333) q[1];
sx q[1];
rz(1.6201409) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.085039728) q[0];
sx q[0];
rz(-1.3360177) q[0];
sx q[0];
rz(-1.1961557) q[0];
rz(-pi) q[1];
rz(2.5720308) q[2];
sx q[2];
rz(-1.6120835) q[2];
sx q[2];
rz(-2.2531367) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7359594) q[1];
sx q[1];
rz(-2.209747) q[1];
sx q[1];
rz(2.6239441) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19100325) q[3];
sx q[3];
rz(-1.2634209) q[3];
sx q[3];
rz(-1.6049811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.3756322) q[2];
sx q[2];
rz(-2.0291294) q[2];
sx q[2];
rz(-2.5755303) q[2];
rz(1.798897) q[3];
sx q[3];
rz(-0.415396) q[3];
sx q[3];
rz(0.28877637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5413496) q[0];
sx q[0];
rz(-0.039027795) q[0];
sx q[0];
rz(1.716123) q[0];
rz(2.0478981) q[1];
sx q[1];
rz(-0.92233557) q[1];
sx q[1];
rz(-2.7519382) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27586994) q[0];
sx q[0];
rz(-1.3148069) q[0];
sx q[0];
rz(-2.0098359) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1946743) q[2];
sx q[2];
rz(-1.3585261) q[2];
sx q[2];
rz(-0.52516261) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4466275) q[1];
sx q[1];
rz(-1.5942861) q[1];
sx q[1];
rz(1.9963677) q[1];
rz(-2.6171754) q[3];
sx q[3];
rz(-2.5811385) q[3];
sx q[3];
rz(-2.2956212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.23715544) q[2];
sx q[2];
rz(-0.5390141) q[2];
sx q[2];
rz(-1.1971486) q[2];
rz(-0.14687471) q[3];
sx q[3];
rz(-2.4683888) q[3];
sx q[3];
rz(-2.1541514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.420153) q[0];
sx q[0];
rz(-2.1099821) q[0];
sx q[0];
rz(-1.4954062) q[0];
rz(-2.738476) q[1];
sx q[1];
rz(-1.8164201) q[1];
sx q[1];
rz(1.4774189) q[1];
rz(1.6942799) q[2];
sx q[2];
rz(-0.99939838) q[2];
sx q[2];
rz(1.8746092) q[2];
rz(-0.38539376) q[3];
sx q[3];
rz(-1.7447532) q[3];
sx q[3];
rz(0.90739653) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
