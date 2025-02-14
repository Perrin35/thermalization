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
rz(-2.1958513) q[1];
sx q[1];
rz(-2.6374964) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61871332) q[0];
sx q[0];
rz(-1.0600495) q[0];
sx q[0];
rz(2.5070531) q[0];
x q[1];
rz(-0.43457793) q[2];
sx q[2];
rz(-0.38056669) q[2];
sx q[2];
rz(-1.5411045) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0054746) q[1];
sx q[1];
rz(-2.1393993) q[1];
sx q[1];
rz(-0.75638812) q[1];
rz(0.13926718) q[3];
sx q[3];
rz(-0.33815171) q[3];
sx q[3];
rz(-0.6247181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7291339) q[2];
sx q[2];
rz(-1.8987741) q[2];
sx q[2];
rz(-0.18623713) q[2];
rz(1.7765744) q[3];
sx q[3];
rz(-2.5297207) q[3];
sx q[3];
rz(-1.2046825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81011009) q[0];
sx q[0];
rz(-2.6832566) q[0];
sx q[0];
rz(-3.0895184) q[0];
rz(-1.0366108) q[1];
sx q[1];
rz(-2.5317445) q[1];
sx q[1];
rz(2.5263272) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0947845) q[0];
sx q[0];
rz(-0.9082709) q[0];
sx q[0];
rz(-0.41562467) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9040163) q[2];
sx q[2];
rz(-1.698809) q[2];
sx q[2];
rz(1.0281171) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2959472) q[1];
sx q[1];
rz(-1.5145647) q[1];
sx q[1];
rz(-0.91744951) q[1];
x q[2];
rz(1.3619945) q[3];
sx q[3];
rz(-2.3376645) q[3];
sx q[3];
rz(1.3458136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.692824) q[2];
sx q[2];
rz(-1.245446) q[2];
sx q[2];
rz(1.2358865) q[2];
rz(2.0125194) q[3];
sx q[3];
rz(-0.90714199) q[3];
sx q[3];
rz(-0.83704078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4835943) q[0];
sx q[0];
rz(-2.3880385) q[0];
sx q[0];
rz(1.765522) q[0];
rz(-0.54028571) q[1];
sx q[1];
rz(-2.1018201) q[1];
sx q[1];
rz(-1.4556063) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7175563) q[0];
sx q[0];
rz(-0.40923318) q[0];
sx q[0];
rz(2.9633029) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7305538) q[2];
sx q[2];
rz(-0.22145311) q[2];
sx q[2];
rz(-2.7483181) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7109622) q[1];
sx q[1];
rz(-2.4050131) q[1];
sx q[1];
rz(2.1107091) q[1];
rz(-pi) q[2];
rz(0.40850477) q[3];
sx q[3];
rz(-0.91445476) q[3];
sx q[3];
rz(2.07043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.49440631) q[2];
sx q[2];
rz(-2.8072085) q[2];
sx q[2];
rz(1.6953267) q[2];
rz(-1.1930195) q[3];
sx q[3];
rz(-2.1665067) q[3];
sx q[3];
rz(-1.6993935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0874264) q[0];
sx q[0];
rz(-2.8717201) q[0];
sx q[0];
rz(-2.7179981) q[0];
rz(2.1845747) q[1];
sx q[1];
rz(-1.3514163) q[1];
sx q[1];
rz(0.097600309) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1621057) q[0];
sx q[0];
rz(-2.477956) q[0];
sx q[0];
rz(-0.62941334) q[0];
x q[1];
rz(2.2636534) q[2];
sx q[2];
rz(-1.4603472) q[2];
sx q[2];
rz(-0.3921961) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4368808) q[1];
sx q[1];
rz(-1.2578674) q[1];
sx q[1];
rz(1.9347316) q[1];
rz(-pi) q[2];
rz(-1.9097435) q[3];
sx q[3];
rz(-2.2710861) q[3];
sx q[3];
rz(-0.46745121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5264954) q[2];
sx q[2];
rz(-0.47250938) q[2];
sx q[2];
rz(0.16461593) q[2];
rz(3.0937255) q[3];
sx q[3];
rz(-1.7367626) q[3];
sx q[3];
rz(2.3544748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1614302) q[0];
sx q[0];
rz(-2.6949368) q[0];
sx q[0];
rz(1.5572146) q[0];
rz(-2.7730675) q[1];
sx q[1];
rz(-1.3216113) q[1];
sx q[1];
rz(-1.535086) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7212352) q[0];
sx q[0];
rz(-2.0518655) q[0];
sx q[0];
rz(-2.1522983) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82442762) q[2];
sx q[2];
rz(-1.2807089) q[2];
sx q[2];
rz(2.1739391) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2531703) q[1];
sx q[1];
rz(-2.0520323) q[1];
sx q[1];
rz(-0.28216823) q[1];
x q[2];
rz(1.3103043) q[3];
sx q[3];
rz(-0.29803571) q[3];
sx q[3];
rz(-2.4233415) q[3];
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
rz(0.16799489) q[3];
sx q[3];
rz(-1.1863656) q[3];
sx q[3];
rz(2.7872938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.0198233) q[0];
sx q[0];
rz(-1.0291809) q[0];
sx q[0];
rz(1.9371012) q[0];
rz(0.80351859) q[1];
sx q[1];
rz(-2.3280227) q[1];
sx q[1];
rz(-0.71570754) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7835282) q[0];
sx q[0];
rz(-0.77371374) q[0];
sx q[0];
rz(0.85229243) q[0];
rz(-pi) q[1];
rz(-1.3722695) q[2];
sx q[2];
rz(-1.1305222) q[2];
sx q[2];
rz(1.0834875) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.83166158) q[1];
sx q[1];
rz(-0.67218053) q[1];
sx q[1];
rz(0.18193717) q[1];
rz(1.1119858) q[3];
sx q[3];
rz(-1.8355287) q[3];
sx q[3];
rz(-2.2609425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7314926) q[2];
sx q[2];
rz(-2.4918753) q[2];
sx q[2];
rz(2.0984207) q[2];
rz(3.0336174) q[3];
sx q[3];
rz(-0.94111809) q[3];
sx q[3];
rz(-2.030355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9589979) q[0];
sx q[0];
rz(-1.7549055) q[0];
sx q[0];
rz(1.9166272) q[0];
rz(2.909868) q[1];
sx q[1];
rz(-2.2544315) q[1];
sx q[1];
rz(1.8909594) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5179493) q[0];
sx q[0];
rz(-0.8276437) q[0];
sx q[0];
rz(1.736899) q[0];
rz(-3.1194889) q[2];
sx q[2];
rz(-1.3726418) q[2];
sx q[2];
rz(1.0304067) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2030484) q[1];
sx q[1];
rz(-0.37751679) q[1];
sx q[1];
rz(-0.63407268) q[1];
rz(-pi) q[2];
rz(0.55646236) q[3];
sx q[3];
rz(-1.5817405) q[3];
sx q[3];
rz(1.3370322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0687678) q[2];
sx q[2];
rz(-0.67049694) q[2];
sx q[2];
rz(3.0094299) q[2];
rz(-0.69563785) q[3];
sx q[3];
rz(-1.1824824) q[3];
sx q[3];
rz(0.80719358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.950133) q[0];
sx q[0];
rz(-1.5810672) q[0];
sx q[0];
rz(-0.70924846) q[0];
rz(-1.6356155) q[1];
sx q[1];
rz(-1.154107) q[1];
sx q[1];
rz(-2.0567315) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1902311) q[0];
sx q[0];
rz(-1.6282646) q[0];
sx q[0];
rz(2.5404853) q[0];
rz(2.9726344) q[2];
sx q[2];
rz(-2.3297133) q[2];
sx q[2];
rz(0.66056992) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.198609) q[1];
sx q[1];
rz(-1.8246905) q[1];
sx q[1];
rz(2.6655546) q[1];
x q[2];
rz(1.510511) q[3];
sx q[3];
rz(-1.5593646) q[3];
sx q[3];
rz(-0.41939467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6214577) q[2];
sx q[2];
rz(-2.1389213) q[2];
sx q[2];
rz(1.8247068) q[2];
rz(2.3560611) q[3];
sx q[3];
rz(-1.9436049) q[3];
sx q[3];
rz(2.7151916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21996466) q[0];
sx q[0];
rz(-1.4861318) q[0];
sx q[0];
rz(-0.40837902) q[0];
rz(-2.1829103) q[1];
sx q[1];
rz(-2.8125693) q[1];
sx q[1];
rz(-1.6201409) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0565529) q[0];
sx q[0];
rz(-1.805575) q[0];
sx q[0];
rz(-1.1961557) q[0];
rz(1.5217811) q[2];
sx q[2];
rz(-1.00178) q[2];
sx q[2];
rz(-0.70876497) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4056333) q[1];
sx q[1];
rz(-2.209747) q[1];
sx q[1];
rz(-2.6239441) q[1];
rz(-1.2580885) q[3];
sx q[3];
rz(-1.3888479) q[3];
sx q[3];
rz(0.092620919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.3756322) q[2];
sx q[2];
rz(-1.1124632) q[2];
sx q[2];
rz(-0.56606236) q[2];
rz(1.3426956) q[3];
sx q[3];
rz(-2.7261966) q[3];
sx q[3];
rz(0.28877637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5413496) q[0];
sx q[0];
rz(-3.1025649) q[0];
sx q[0];
rz(1.4254697) q[0];
rz(2.0478981) q[1];
sx q[1];
rz(-0.92233557) q[1];
sx q[1];
rz(0.38965449) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3411639) q[0];
sx q[0];
rz(-0.50398705) q[0];
sx q[0];
rz(-1.0188854) q[0];
x q[1];
rz(2.1946743) q[2];
sx q[2];
rz(-1.7830666) q[2];
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
rz(pi/2) q[0];
rz(0.07239711) q[1];
sx q[1];
rz(-2.715413) q[1];
sx q[1];
rz(-1.5139493) q[1];
rz(0.49764244) q[3];
sx q[3];
rz(-1.3013869) q[3];
sx q[3];
rz(-2.8723418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9044372) q[2];
sx q[2];
rz(-0.5390141) q[2];
sx q[2];
rz(-1.944444) q[2];
rz(0.14687471) q[3];
sx q[3];
rz(-2.4683888) q[3];
sx q[3];
rz(2.1541514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7214397) q[0];
sx q[0];
rz(-1.0316105) q[0];
sx q[0];
rz(1.6461865) q[0];
rz(2.738476) q[1];
sx q[1];
rz(-1.3251726) q[1];
sx q[1];
rz(-1.6641738) q[1];
rz(0.57488048) q[2];
sx q[2];
rz(-1.6745865) q[2];
sx q[2];
rz(-2.9048017) q[2];
rz(2.7561989) q[3];
sx q[3];
rz(-1.7447532) q[3];
sx q[3];
rz(0.90739653) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
