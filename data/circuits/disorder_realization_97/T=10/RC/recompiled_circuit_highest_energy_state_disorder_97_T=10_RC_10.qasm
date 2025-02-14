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
rz(-1.1345743) q[0];
sx q[0];
rz(1.467508) q[0];
rz(1.6788586) q[1];
sx q[1];
rz(0.94574133) q[1];
sx q[1];
rz(8.9206817) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5228793) q[0];
sx q[0];
rz(-1.0600495) q[0];
sx q[0];
rz(0.63453959) q[0];
x q[1];
rz(0.43457793) q[2];
sx q[2];
rz(-2.761026) q[2];
sx q[2];
rz(1.6004882) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.91612736) q[1];
sx q[1];
rz(-2.230495) q[1];
sx q[1];
rz(0.74972357) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8064742) q[3];
sx q[3];
rz(-1.6168645) q[3];
sx q[3];
rz(-2.3269881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.41245875) q[2];
sx q[2];
rz(-1.8987741) q[2];
sx q[2];
rz(0.18623713) q[2];
rz(-1.3650182) q[3];
sx q[3];
rz(-2.5297207) q[3];
sx q[3];
rz(1.9369102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3314826) q[0];
sx q[0];
rz(-2.6832566) q[0];
sx q[0];
rz(-0.052074281) q[0];
rz(2.1049818) q[1];
sx q[1];
rz(-2.5317445) q[1];
sx q[1];
rz(-0.61526543) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8826854) q[0];
sx q[0];
rz(-1.2468161) q[0];
sx q[0];
rz(0.86470226) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7024666) q[2];
sx q[2];
rz(-1.3352009) q[2];
sx q[2];
rz(2.5680096) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.6518636) q[1];
sx q[1];
rz(-2.4861838) q[1];
sx q[1];
rz(-1.6631399) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3637487) q[3];
sx q[3];
rz(-1.7206216) q[3];
sx q[3];
rz(-2.7706287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.692824) q[2];
sx q[2];
rz(-1.8961467) q[2];
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
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6579984) q[0];
sx q[0];
rz(-2.3880385) q[0];
sx q[0];
rz(-1.765522) q[0];
rz(2.6013069) q[1];
sx q[1];
rz(-2.1018201) q[1];
sx q[1];
rz(-1.4556063) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7175563) q[0];
sx q[0];
rz(-0.40923318) q[0];
sx q[0];
rz(0.17828973) q[0];
rz(0.41103883) q[2];
sx q[2];
rz(-2.9201395) q[2];
sx q[2];
rz(-2.7483181) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.27764717) q[1];
sx q[1];
rz(-1.2182115) q[1];
sx q[1];
rz(2.2318798) q[1];
rz(-pi) q[2];
rz(-0.40850477) q[3];
sx q[3];
rz(-0.91445476) q[3];
sx q[3];
rz(1.0711627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6471863) q[2];
sx q[2];
rz(-0.33438412) q[2];
sx q[2];
rz(1.6953267) q[2];
rz(-1.9485731) q[3];
sx q[3];
rz(-0.97508591) q[3];
sx q[3];
rz(1.4421991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(1.0541662) q[0];
sx q[0];
rz(-0.26987258) q[0];
sx q[0];
rz(0.42359459) q[0];
rz(-2.1845747) q[1];
sx q[1];
rz(-1.7901763) q[1];
sx q[1];
rz(0.097600309) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4159629) q[0];
sx q[0];
rz(-1.0495674) q[0];
sx q[0];
rz(1.1394016) q[0];
rz(-pi) q[1];
rz(-2.2636534) q[2];
sx q[2];
rz(-1.6812455) q[2];
sx q[2];
rz(-0.3921961) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4368808) q[1];
sx q[1];
rz(-1.8837253) q[1];
sx q[1];
rz(-1.2068611) q[1];
rz(-2.41231) q[3];
sx q[3];
rz(-1.3137307) q[3];
sx q[3];
rz(2.2616539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5264954) q[2];
sx q[2];
rz(-2.6690833) q[2];
sx q[2];
rz(-2.9769767) q[2];
rz(-3.0937255) q[3];
sx q[3];
rz(-1.4048301) q[3];
sx q[3];
rz(2.3544748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1614302) q[0];
sx q[0];
rz(-0.4466559) q[0];
sx q[0];
rz(1.5572146) q[0];
rz(2.7730675) q[1];
sx q[1];
rz(-1.8199814) q[1];
sx q[1];
rz(-1.535086) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6958862) q[0];
sx q[0];
rz(-1.0622171) q[0];
sx q[0];
rz(2.5832547) q[0];
rz(2.317165) q[2];
sx q[2];
rz(-1.2807089) q[2];
sx q[2];
rz(-2.1739391) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6905744) q[1];
sx q[1];
rz(-1.8201882) q[1];
sx q[1];
rz(2.0687769) q[1];
x q[2];
rz(-1.3103043) q[3];
sx q[3];
rz(-0.29803571) q[3];
sx q[3];
rz(-0.71825114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0198233) q[0];
sx q[0];
rz(-2.1124117) q[0];
sx q[0];
rz(-1.2044915) q[0];
rz(-0.80351859) q[1];
sx q[1];
rz(-2.3280227) q[1];
sx q[1];
rz(0.71570754) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34622231) q[0];
sx q[0];
rz(-2.0487794) q[0];
sx q[0];
rz(-2.204889) q[0];
x q[1];
rz(-2.6936221) q[2];
sx q[2];
rz(-1.391419) q[2];
sx q[2];
rz(0.57283869) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3099311) q[1];
sx q[1];
rz(-0.67218053) q[1];
sx q[1];
rz(-0.18193717) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0296069) q[3];
sx q[3];
rz(-1.306064) q[3];
sx q[3];
rz(0.88065016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.41010007) q[2];
sx q[2];
rz(-0.64971739) q[2];
sx q[2];
rz(-2.0984207) q[2];
rz(-0.10797524) q[3];
sx q[3];
rz(-2.2004746) q[3];
sx q[3];
rz(-1.1112377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9589979) q[0];
sx q[0];
rz(-1.7549055) q[0];
sx q[0];
rz(-1.2249655) q[0];
rz(0.23172465) q[1];
sx q[1];
rz(-2.2544315) q[1];
sx q[1];
rz(1.2506332) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38076048) q[0];
sx q[0];
rz(-0.75801132) q[0];
sx q[0];
rz(-2.9635628) q[0];
x q[1];
rz(0.022103775) q[2];
sx q[2];
rz(-1.7689509) q[2];
sx q[2];
rz(2.1111859) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2030484) q[1];
sx q[1];
rz(-0.37751679) q[1];
sx q[1];
rz(0.63407268) q[1];
rz(-pi) q[2];
rz(2.5851303) q[3];
sx q[3];
rz(-1.5598522) q[3];
sx q[3];
rz(-1.8045604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0687678) q[2];
sx q[2];
rz(-2.4710957) q[2];
sx q[2];
rz(3.0094299) q[2];
rz(-0.69563785) q[3];
sx q[3];
rz(-1.1824824) q[3];
sx q[3];
rz(-2.3343991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19145963) q[0];
sx q[0];
rz(-1.5605254) q[0];
sx q[0];
rz(2.4323442) q[0];
rz(1.6356155) q[1];
sx q[1];
rz(-1.154107) q[1];
sx q[1];
rz(-1.0848612) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70299545) q[0];
sx q[0];
rz(-2.5380822) q[0];
sx q[0];
rz(0.10137697) q[0];
rz(-0.16895825) q[2];
sx q[2];
rz(-0.8118793) q[2];
sx q[2];
rz(2.4810227) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6425655) q[1];
sx q[1];
rz(-1.111218) q[1];
sx q[1];
rz(-1.8548555) q[1];
x q[2];
rz(-1.6310817) q[3];
sx q[3];
rz(-1.5593646) q[3];
sx q[3];
rz(2.722198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6214577) q[2];
sx q[2];
rz(-2.1389213) q[2];
sx q[2];
rz(1.3168859) q[2];
rz(-0.78553158) q[3];
sx q[3];
rz(-1.1979878) q[3];
sx q[3];
rz(-2.7151916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.921628) q[0];
sx q[0];
rz(-1.6554609) q[0];
sx q[0];
rz(0.40837902) q[0];
rz(-2.1829103) q[1];
sx q[1];
rz(-0.32902333) q[1];
sx q[1];
rz(1.6201409) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.085039728) q[0];
sx q[0];
rz(-1.3360177) q[0];
sx q[0];
rz(1.1961557) q[0];
rz(1.6198115) q[2];
sx q[2];
rz(-1.00178) q[2];
sx q[2];
rz(-2.4328277) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1679528) q[1];
sx q[1];
rz(-0.79881891) q[1];
sx q[1];
rz(-0.98320191) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9505894) q[3];
sx q[3];
rz(-1.8781717) q[3];
sx q[3];
rz(-1.6049811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7659605) q[2];
sx q[2];
rz(-2.0291294) q[2];
sx q[2];
rz(-0.56606236) q[2];
rz(1.798897) q[3];
sx q[3];
rz(-2.7261966) q[3];
sx q[3];
rz(-0.28877637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5413496) q[0];
sx q[0];
rz(-0.039027795) q[0];
sx q[0];
rz(-1.4254697) q[0];
rz(1.0936945) q[1];
sx q[1];
rz(-0.92233557) q[1];
sx q[1];
rz(2.7519382) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27586994) q[0];
sx q[0];
rz(-1.3148069) q[0];
sx q[0];
rz(-1.1317568) q[0];
x q[1];
rz(-2.8820428) q[2];
sx q[2];
rz(-0.96298157) q[2];
sx q[2];
rz(1.1961301) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4466275) q[1];
sx q[1];
rz(-1.5942861) q[1];
sx q[1];
rz(-1.9963677) q[1];
rz(-pi) q[2];
rz(-2.6439502) q[3];
sx q[3];
rz(-1.8402057) q[3];
sx q[3];
rz(-0.26925081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
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
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.7214397) q[0];
sx q[0];
rz(-1.0316105) q[0];
sx q[0];
rz(1.6461865) q[0];
rz(0.40311665) q[1];
sx q[1];
rz(-1.8164201) q[1];
sx q[1];
rz(1.4774189) q[1];
rz(-1.4473128) q[2];
sx q[2];
rz(-0.99939838) q[2];
sx q[2];
rz(1.8746092) q[2];
rz(0.38539376) q[3];
sx q[3];
rz(-1.3968395) q[3];
sx q[3];
rz(-2.2341961) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
