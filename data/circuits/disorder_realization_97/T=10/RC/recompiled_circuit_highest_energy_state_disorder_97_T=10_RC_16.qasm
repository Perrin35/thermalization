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
rz(0.94574133) q[1];
sx q[1];
rz(8.9206817) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.534908) q[0];
sx q[0];
rz(-1.0272756) q[0];
sx q[0];
rz(-2.1786819) q[0];
rz(-pi) q[1];
rz(0.43457793) q[2];
sx q[2];
rz(-0.38056669) q[2];
sx q[2];
rz(1.5411045) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2368184) q[1];
sx q[1];
rz(-2.1874551) q[1];
sx q[1];
rz(-2.2916433) q[1];
x q[2];
rz(-0.13926718) q[3];
sx q[3];
rz(-0.33815171) q[3];
sx q[3];
rz(0.6247181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7291339) q[2];
sx q[2];
rz(-1.2428186) q[2];
sx q[2];
rz(2.9553555) q[2];
rz(-1.7765744) q[3];
sx q[3];
rz(-2.5297207) q[3];
sx q[3];
rz(-1.9369102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81011009) q[0];
sx q[0];
rz(-0.45833603) q[0];
sx q[0];
rz(3.0895184) q[0];
rz(2.1049818) q[1];
sx q[1];
rz(-2.5317445) q[1];
sx q[1];
rz(2.5263272) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2589073) q[0];
sx q[0];
rz(-1.8947766) q[0];
sx q[0];
rz(2.2768904) q[0];
rz(1.7024666) q[2];
sx q[2];
rz(-1.8063917) q[2];
sx q[2];
rz(-0.57358303) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.76814684) q[1];
sx q[1];
rz(-2.2229338) q[1];
sx q[1];
rz(3.0708205) q[1];
rz(-pi) q[2];
rz(-2.9297057) q[3];
sx q[3];
rz(-0.78915262) q[3];
sx q[3];
rz(2.0921538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.692824) q[2];
sx q[2];
rz(-1.245446) q[2];
sx q[2];
rz(-1.9057062) q[2];
rz(1.1290733) q[3];
sx q[3];
rz(-0.90714199) q[3];
sx q[3];
rz(0.83704078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6579984) q[0];
sx q[0];
rz(-2.3880385) q[0];
sx q[0];
rz(-1.765522) q[0];
rz(-2.6013069) q[1];
sx q[1];
rz(-1.0397725) q[1];
sx q[1];
rz(1.6859863) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8309915) q[0];
sx q[0];
rz(-1.5001703) q[0];
sx q[0];
rz(-2.7381606) q[0];
rz(-pi) q[1];
rz(0.41103883) q[2];
sx q[2];
rz(-2.9201395) q[2];
sx q[2];
rz(-2.7483181) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0307349) q[1];
sx q[1];
rz(-0.95672119) q[1];
sx q[1];
rz(-0.43621896) q[1];
rz(-pi) q[2];
rz(2.0469401) q[3];
sx q[3];
rz(-2.3848001) q[3];
sx q[3];
rz(1.4534637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.49440631) q[2];
sx q[2];
rz(-2.8072085) q[2];
sx q[2];
rz(-1.4462659) q[2];
rz(1.9485731) q[3];
sx q[3];
rz(-0.97508591) q[3];
sx q[3];
rz(-1.4421991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0874264) q[0];
sx q[0];
rz(-0.26987258) q[0];
sx q[0];
rz(0.42359459) q[0];
rz(-2.1845747) q[1];
sx q[1];
rz(-1.3514163) q[1];
sx q[1];
rz(-0.097600309) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1621057) q[0];
sx q[0];
rz(-0.66363664) q[0];
sx q[0];
rz(2.5121793) q[0];
x q[1];
rz(1.7427069) q[2];
sx q[2];
rz(-2.4414276) q[2];
sx q[2];
rz(-1.3106048) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3279539) q[1];
sx q[1];
rz(-0.47537741) q[1];
sx q[1];
rz(-0.83303501) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.72928263) q[3];
sx q[3];
rz(-1.8278619) q[3];
sx q[3];
rz(2.2616539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5264954) q[2];
sx q[2];
rz(-0.47250938) q[2];
sx q[2];
rz(-0.16461593) q[2];
rz(-3.0937255) q[3];
sx q[3];
rz(-1.4048301) q[3];
sx q[3];
rz(2.3544748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1614302) q[0];
sx q[0];
rz(-0.4466559) q[0];
sx q[0];
rz(-1.5572146) q[0];
rz(-2.7730675) q[1];
sx q[1];
rz(-1.8199814) q[1];
sx q[1];
rz(1.535086) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6958862) q[0];
sx q[0];
rz(-2.0793756) q[0];
sx q[0];
rz(-2.5832547) q[0];
x q[1];
rz(-2.755411) q[2];
sx q[2];
rz(-0.86241041) q[2];
sx q[2];
rz(0.86175534) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.45101825) q[1];
sx q[1];
rz(-1.8201882) q[1];
sx q[1];
rz(1.0728157) q[1];
rz(-pi) q[2];
rz(1.2822578) q[3];
sx q[3];
rz(-1.6464982) q[3];
sx q[3];
rz(-2.5385365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2313472) q[2];
sx q[2];
rz(-2.7635837) q[2];
sx q[2];
rz(2.0965915) q[2];
rz(2.9735978) q[3];
sx q[3];
rz(-1.1863656) q[3];
sx q[3];
rz(-2.7872938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1217693) q[0];
sx q[0];
rz(-1.0291809) q[0];
sx q[0];
rz(1.9371012) q[0];
rz(-2.3380741) q[1];
sx q[1];
rz(-0.81356994) q[1];
sx q[1];
rz(-2.4258851) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7953703) q[0];
sx q[0];
rz(-1.0928133) q[0];
sx q[0];
rz(-0.93670364) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7693232) q[2];
sx q[2];
rz(-1.1305222) q[2];
sx q[2];
rz(1.0834875) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.83166158) q[1];
sx q[1];
rz(-2.4694121) q[1];
sx q[1];
rz(-0.18193717) q[1];
rz(-0.29362595) q[3];
sx q[3];
rz(-2.01247) q[3];
sx q[3];
rz(0.81868086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.41010007) q[2];
sx q[2];
rz(-0.64971739) q[2];
sx q[2];
rz(2.0984207) q[2];
rz(-3.0336174) q[3];
sx q[3];
rz(-0.94111809) q[3];
sx q[3];
rz(-1.1112377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1825948) q[0];
sx q[0];
rz(-1.7549055) q[0];
sx q[0];
rz(-1.9166272) q[0];
rz(-2.909868) q[1];
sx q[1];
rz(-2.2544315) q[1];
sx q[1];
rz(1.2506332) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5179493) q[0];
sx q[0];
rz(-2.3139489) q[0];
sx q[0];
rz(1.4046937) q[0];
rz(-pi) q[1];
rz(-1.6804303) q[2];
sx q[2];
rz(-2.9422252) q[2];
sx q[2];
rz(2.2230172) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.53377235) q[1];
sx q[1];
rz(-1.8723066) q[1];
sx q[1];
rz(1.34006) q[1];
rz(-0.020718109) q[3];
sx q[3];
rz(-2.585034) q[3];
sx q[3];
rz(-2.8902365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.072824868) q[2];
sx q[2];
rz(-0.67049694) q[2];
sx q[2];
rz(3.0094299) q[2];
rz(0.69563785) q[3];
sx q[3];
rz(-1.1824824) q[3];
sx q[3];
rz(2.3343991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19145963) q[0];
sx q[0];
rz(-1.5605254) q[0];
sx q[0];
rz(0.70924846) q[0];
rz(-1.5059772) q[1];
sx q[1];
rz(-1.9874856) q[1];
sx q[1];
rz(-2.0567315) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1902311) q[0];
sx q[0];
rz(-1.5133281) q[0];
sx q[0];
rz(2.5404853) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9726344) q[2];
sx q[2];
rz(-2.3297133) q[2];
sx q[2];
rz(2.4810227) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.9429837) q[1];
sx q[1];
rz(-1.3169022) q[1];
sx q[1];
rz(0.47603807) q[1];
x q[2];
rz(-3.1301401) q[3];
sx q[3];
rz(-1.6310777) q[3];
sx q[3];
rz(-1.1520916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.52013493) q[2];
sx q[2];
rz(-1.0026714) q[2];
sx q[2];
rz(1.3168859) q[2];
rz(0.78553158) q[3];
sx q[3];
rz(-1.9436049) q[3];
sx q[3];
rz(0.42640105) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.921628) q[0];
sx q[0];
rz(-1.4861318) q[0];
sx q[0];
rz(-2.7332136) q[0];
rz(0.95868239) q[1];
sx q[1];
rz(-2.8125693) q[1];
sx q[1];
rz(-1.6201409) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7470541) q[0];
sx q[0];
rz(-1.9346721) q[0];
sx q[0];
rz(0.25157148) q[0];
rz(-pi) q[1];
rz(-1.6198115) q[2];
sx q[2];
rz(-1.00178) q[2];
sx q[2];
rz(-0.70876497) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.97363981) q[1];
sx q[1];
rz(-0.79881891) q[1];
sx q[1];
rz(-0.98320191) q[1];
x q[2];
rz(-1.8835041) q[3];
sx q[3];
rz(-1.3888479) q[3];
sx q[3];
rz(-0.092620919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.3756322) q[2];
sx q[2];
rz(-2.0291294) q[2];
sx q[2];
rz(0.56606236) q[2];
rz(-1.3426956) q[3];
sx q[3];
rz(-0.415396) q[3];
sx q[3];
rz(-2.8528163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.600243) q[0];
sx q[0];
rz(-3.1025649) q[0];
sx q[0];
rz(-1.716123) q[0];
rz(2.0478981) q[1];
sx q[1];
rz(-0.92233557) q[1];
sx q[1];
rz(0.38965449) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80042875) q[0];
sx q[0];
rz(-0.50398705) q[0];
sx q[0];
rz(-1.0188854) q[0];
rz(-1.924224) q[2];
sx q[2];
rz(-2.4871748) q[2];
sx q[2];
rz(-0.76088727) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.13481465) q[1];
sx q[1];
rz(-1.9962427) q[1];
sx q[1];
rz(3.1158034) q[1];
rz(-2.6439502) q[3];
sx q[3];
rz(-1.3013869) q[3];
sx q[3];
rz(0.26925081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23715544) q[2];
sx q[2];
rz(-2.6025786) q[2];
sx q[2];
rz(1.1971486) q[2];
rz(0.14687471) q[3];
sx q[3];
rz(-2.4683888) q[3];
sx q[3];
rz(-0.98744121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.420153) q[0];
sx q[0];
rz(-2.1099821) q[0];
sx q[0];
rz(-1.4954062) q[0];
rz(2.738476) q[1];
sx q[1];
rz(-1.3251726) q[1];
sx q[1];
rz(-1.6641738) q[1];
rz(-0.57488048) q[2];
sx q[2];
rz(-1.4670062) q[2];
sx q[2];
rz(0.23679096) q[2];
rz(0.4372863) q[3];
sx q[3];
rz(-0.42103817) q[3];
sx q[3];
rz(-0.26013817) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
