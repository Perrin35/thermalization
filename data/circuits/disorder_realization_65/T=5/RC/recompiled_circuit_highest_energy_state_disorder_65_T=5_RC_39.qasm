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
rz(-1.3574358) q[0];
sx q[0];
rz(-2.5166002) q[0];
sx q[0];
rz(-1.5224737) q[0];
rz(-1.9982665) q[1];
sx q[1];
rz(-0.63653094) q[1];
sx q[1];
rz(-1.7826537) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1899043) q[0];
sx q[0];
rz(-2.1225327) q[0];
sx q[0];
rz(2.8008528) q[0];
x q[1];
rz(0.29523103) q[2];
sx q[2];
rz(-1.9257987) q[2];
sx q[2];
rz(-0.4685185) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.93329079) q[1];
sx q[1];
rz(-1.33889) q[1];
sx q[1];
rz(1.307322) q[1];
x q[2];
rz(3.0380958) q[3];
sx q[3];
rz(-2.0465133) q[3];
sx q[3];
rz(2.5535943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6978567) q[2];
sx q[2];
rz(-2.3190658) q[2];
sx q[2];
rz(-2.6653384) q[2];
rz(-0.13992986) q[3];
sx q[3];
rz(-1.6401451) q[3];
sx q[3];
rz(-3.068315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78854617) q[0];
sx q[0];
rz(-1.6206425) q[0];
sx q[0];
rz(0.36343685) q[0];
rz(0.67047554) q[1];
sx q[1];
rz(-1.3251708) q[1];
sx q[1];
rz(-1.0796116) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.712942) q[0];
sx q[0];
rz(-1.043937) q[0];
sx q[0];
rz(-2.8654979) q[0];
rz(-pi) q[1];
rz(2.1035467) q[2];
sx q[2];
rz(-1.4584625) q[2];
sx q[2];
rz(2.7137604) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0091925) q[1];
sx q[1];
rz(-1.3295439) q[1];
sx q[1];
rz(-3.0962837) q[1];
rz(-pi) q[2];
rz(0.3877181) q[3];
sx q[3];
rz(-0.61655515) q[3];
sx q[3];
rz(-1.316837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4662027) q[2];
sx q[2];
rz(-0.73489302) q[2];
sx q[2];
rz(2.7030763) q[2];
rz(1.0112666) q[3];
sx q[3];
rz(-0.68532419) q[3];
sx q[3];
rz(1.4810168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48280516) q[0];
sx q[0];
rz(-2.6130982) q[0];
sx q[0];
rz(-0.82861376) q[0];
rz(1.2476745) q[1];
sx q[1];
rz(-1.9745461) q[1];
sx q[1];
rz(2.8819328) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22406507) q[0];
sx q[0];
rz(-0.7389141) q[0];
sx q[0];
rz(-1.9784443) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9405878) q[2];
sx q[2];
rz(-1.9281053) q[2];
sx q[2];
rz(1.0543038) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3500525) q[1];
sx q[1];
rz(-2.2660221) q[1];
sx q[1];
rz(-0.94757845) q[1];
x q[2];
rz(-2.8494495) q[3];
sx q[3];
rz(-0.93860561) q[3];
sx q[3];
rz(1.8504253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9344249) q[2];
sx q[2];
rz(-1.577689) q[2];
sx q[2];
rz(-2.3210607) q[2];
rz(2.3151243) q[3];
sx q[3];
rz(-0.99244899) q[3];
sx q[3];
rz(-1.3124527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8823223) q[0];
sx q[0];
rz(-2.3311908) q[0];
sx q[0];
rz(-0.93233863) q[0];
rz(1.8238292) q[1];
sx q[1];
rz(-1.1612786) q[1];
sx q[1];
rz(-3.0188149) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30194382) q[0];
sx q[0];
rz(-2.1084063) q[0];
sx q[0];
rz(-1.4728236) q[0];
rz(2.8343646) q[2];
sx q[2];
rz(-2.3005552) q[2];
sx q[2];
rz(-2.407069) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.40298324) q[1];
sx q[1];
rz(-1.5697641) q[1];
sx q[1];
rz(-2.3139075) q[1];
x q[2];
rz(0.56468954) q[3];
sx q[3];
rz(-1.3211234) q[3];
sx q[3];
rz(2.3752729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.54470283) q[2];
sx q[2];
rz(-0.92146102) q[2];
sx q[2];
rz(-0.60453647) q[2];
rz(-0.70945159) q[3];
sx q[3];
rz(-2.3716898) q[3];
sx q[3];
rz(-2.3179222) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9619047) q[0];
sx q[0];
rz(-3.0178495) q[0];
sx q[0];
rz(1.7748348) q[0];
rz(-2.1429515) q[1];
sx q[1];
rz(-2.3256681) q[1];
sx q[1];
rz(2.6008115) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3546412) q[0];
sx q[0];
rz(-1.4313414) q[0];
sx q[0];
rz(2.9324173) q[0];
rz(-2.3930757) q[2];
sx q[2];
rz(-2.243804) q[2];
sx q[2];
rz(1.7130913) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0980254) q[1];
sx q[1];
rz(-0.35570733) q[1];
sx q[1];
rz(0.58009899) q[1];
rz(2.5752047) q[3];
sx q[3];
rz(-1.4175804) q[3];
sx q[3];
rz(-0.18463102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46405408) q[2];
sx q[2];
rz(-0.81391922) q[2];
sx q[2];
rz(0.82826725) q[2];
rz(-1.6895435) q[3];
sx q[3];
rz(-1.4938846) q[3];
sx q[3];
rz(2.233708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0238817) q[0];
sx q[0];
rz(-1.9458867) q[0];
sx q[0];
rz(-1.9697795) q[0];
rz(-0.68012971) q[1];
sx q[1];
rz(-1.9155733) q[1];
sx q[1];
rz(2.5853058) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9754575) q[0];
sx q[0];
rz(-1.8248471) q[0];
sx q[0];
rz(1.9587166) q[0];
rz(-3.1388144) q[2];
sx q[2];
rz(-1.4242715) q[2];
sx q[2];
rz(0.79039449) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.57842931) q[1];
sx q[1];
rz(-2.142504) q[1];
sx q[1];
rz(1.494638) q[1];
rz(-1.0917967) q[3];
sx q[3];
rz(-0.99078876) q[3];
sx q[3];
rz(0.73778462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3351626) q[2];
sx q[2];
rz(-1.8589636) q[2];
sx q[2];
rz(0.20509091) q[2];
rz(-2.3388376) q[3];
sx q[3];
rz(-1.1057248) q[3];
sx q[3];
rz(0.55454379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4270808) q[0];
sx q[0];
rz(-1.0733805) q[0];
sx q[0];
rz(-0.22739534) q[0];
rz(2.3511476) q[1];
sx q[1];
rz(-0.77722725) q[1];
sx q[1];
rz(0.8367742) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6444708) q[0];
sx q[0];
rz(-2.5724412) q[0];
sx q[0];
rz(-0.41438024) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.927194) q[2];
sx q[2];
rz(-2.1292392) q[2];
sx q[2];
rz(1.3791549) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0322276) q[1];
sx q[1];
rz(-1.9739474) q[1];
sx q[1];
rz(-1.228986) q[1];
rz(-1.1985336) q[3];
sx q[3];
rz(-2.296631) q[3];
sx q[3];
rz(1.6282631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5922015) q[2];
sx q[2];
rz(-1.4791919) q[2];
sx q[2];
rz(0.57360348) q[2];
rz(2.8163689) q[3];
sx q[3];
rz(-0.89541382) q[3];
sx q[3];
rz(0.4044683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0322872) q[0];
sx q[0];
rz(-2.9635297) q[0];
sx q[0];
rz(-0.67894116) q[0];
rz(-3.0067054) q[1];
sx q[1];
rz(-2.0811681) q[1];
sx q[1];
rz(-1.7178242) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8979864) q[0];
sx q[0];
rz(-2.4084414) q[0];
sx q[0];
rz(2.2871458) q[0];
rz(0.94332327) q[2];
sx q[2];
rz(-2.3296142) q[2];
sx q[2];
rz(1.2710149) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9805124) q[1];
sx q[1];
rz(-0.80787611) q[1];
sx q[1];
rz(0.041461583) q[1];
rz(-pi) q[2];
rz(-2.1072949) q[3];
sx q[3];
rz(-0.52742672) q[3];
sx q[3];
rz(2.3303383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5474825) q[2];
sx q[2];
rz(-2.2308733) q[2];
sx q[2];
rz(-0.2317079) q[2];
rz(2.1314651) q[3];
sx q[3];
rz(-1.2332656) q[3];
sx q[3];
rz(-0.89099425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4657087) q[0];
sx q[0];
rz(-0.75513419) q[0];
sx q[0];
rz(2.6278507) q[0];
rz(-1.4328009) q[1];
sx q[1];
rz(-1.865973) q[1];
sx q[1];
rz(1.6339711) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9378273) q[0];
sx q[0];
rz(-1.4368334) q[0];
sx q[0];
rz(-1.6197657) q[0];
x q[1];
rz(2.8880226) q[2];
sx q[2];
rz(-1.7268983) q[2];
sx q[2];
rz(1.9954322) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.53879879) q[1];
sx q[1];
rz(-1.4414296) q[1];
sx q[1];
rz(-1.4461229) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8225372) q[3];
sx q[3];
rz(-1.5466403) q[3];
sx q[3];
rz(0.59059737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.12019176) q[2];
sx q[2];
rz(-0.84453619) q[2];
sx q[2];
rz(0.25944844) q[2];
rz(-1.9868959) q[3];
sx q[3];
rz(-0.76609937) q[3];
sx q[3];
rz(1.8416789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65086377) q[0];
sx q[0];
rz(-1.6614953) q[0];
sx q[0];
rz(-1.4890626) q[0];
rz(2.5015855) q[1];
sx q[1];
rz(-1.0970486) q[1];
sx q[1];
rz(-2.1515813) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31371597) q[0];
sx q[0];
rz(-1.4826688) q[0];
sx q[0];
rz(0.14493305) q[0];
rz(1.1009388) q[2];
sx q[2];
rz(-2.2820916) q[2];
sx q[2];
rz(-2.2180722) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9089936) q[1];
sx q[1];
rz(-0.4505583) q[1];
sx q[1];
rz(-1.1934936) q[1];
x q[2];
rz(3.1166273) q[3];
sx q[3];
rz(-1.7933266) q[3];
sx q[3];
rz(-0.74786598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2085569) q[2];
sx q[2];
rz(-1.8427589) q[2];
sx q[2];
rz(-2.9162858) q[2];
rz(2.2202668) q[3];
sx q[3];
rz(-0.39042979) q[3];
sx q[3];
rz(3.0909753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0113572) q[0];
sx q[0];
rz(-2.7738032) q[0];
sx q[0];
rz(-2.2875447) q[0];
rz(1.8232952) q[1];
sx q[1];
rz(-1.6164936) q[1];
sx q[1];
rz(2.2093538) q[1];
rz(-1.9226045) q[2];
sx q[2];
rz(-2.8657037) q[2];
sx q[2];
rz(-2.5421179) q[2];
rz(-0.18410889) q[3];
sx q[3];
rz(-1.0297903) q[3];
sx q[3];
rz(-0.31576706) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
