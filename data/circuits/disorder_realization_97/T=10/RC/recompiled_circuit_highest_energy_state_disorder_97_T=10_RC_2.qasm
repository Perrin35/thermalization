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
rz(-0.47973862) q[0];
sx q[0];
rz(-2.0070183) q[0];
sx q[0];
rz(-1.467508) q[0];
rz(1.6788586) q[1];
sx q[1];
rz(0.94574133) q[1];
sx q[1];
rz(8.9206817) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.534908) q[0];
sx q[0];
rz(-1.0272756) q[0];
sx q[0];
rz(-0.96291079) q[0];
x q[1];
rz(1.4039223) q[2];
sx q[2];
rz(-1.9144399) q[2];
sx q[2];
rz(2.0640896) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0054746) q[1];
sx q[1];
rz(-2.1393993) q[1];
sx q[1];
rz(-2.3852045) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8064742) q[3];
sx q[3];
rz(-1.5247282) q[3];
sx q[3];
rz(-2.3269881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.41245875) q[2];
sx q[2];
rz(-1.2428186) q[2];
sx q[2];
rz(2.9553555) q[2];
rz(-1.3650182) q[3];
sx q[3];
rz(-2.5297207) q[3];
sx q[3];
rz(1.9369102) q[3];
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
rz(2.3314826) q[0];
sx q[0];
rz(-0.45833603) q[0];
sx q[0];
rz(3.0895184) q[0];
rz(-1.0366108) q[1];
sx q[1];
rz(-2.5317445) q[1];
sx q[1];
rz(2.5263272) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2589073) q[0];
sx q[0];
rz(-1.8947766) q[0];
sx q[0];
rz(-0.86470226) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6411166) q[2];
sx q[2];
rz(-2.8722974) q[2];
sx q[2];
rz(-3.0840906) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.6518636) q[1];
sx q[1];
rz(-2.4861838) q[1];
sx q[1];
rz(1.6631399) q[1];
rz(2.9297057) q[3];
sx q[3];
rz(-2.35244) q[3];
sx q[3];
rz(2.0921538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.692824) q[2];
sx q[2];
rz(-1.245446) q[2];
sx q[2];
rz(-1.2358865) q[2];
rz(-2.0125194) q[3];
sx q[3];
rz(-2.2344507) q[3];
sx q[3];
rz(-0.83704078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4835943) q[0];
sx q[0];
rz(-0.75355419) q[0];
sx q[0];
rz(-1.765522) q[0];
rz(-0.54028571) q[1];
sx q[1];
rz(-2.1018201) q[1];
sx q[1];
rz(1.6859863) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9115096) q[0];
sx q[0];
rz(-1.1684281) q[0];
sx q[0];
rz(-1.6475639) q[0];
x q[1];
rz(0.20353453) q[2];
sx q[2];
rz(-1.4829205) q[2];
sx q[2];
rz(-2.3661302) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4306304) q[1];
sx q[1];
rz(-2.4050131) q[1];
sx q[1];
rz(1.0308835) q[1];
rz(2.0469401) q[3];
sx q[3];
rz(-0.75679251) q[3];
sx q[3];
rz(1.688129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6471863) q[2];
sx q[2];
rz(-2.8072085) q[2];
sx q[2];
rz(-1.4462659) q[2];
rz(-1.9485731) q[3];
sx q[3];
rz(-2.1665067) q[3];
sx q[3];
rz(1.6993935) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0541662) q[0];
sx q[0];
rz(-2.8717201) q[0];
sx q[0];
rz(-2.7179981) q[0];
rz(-2.1845747) q[1];
sx q[1];
rz(-1.3514163) q[1];
sx q[1];
rz(3.0439923) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4159629) q[0];
sx q[0];
rz(-2.0920252) q[0];
sx q[0];
rz(1.1394016) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9984442) q[2];
sx q[2];
rz(-0.88298702) q[2];
sx q[2];
rz(2.0542415) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7047119) q[1];
sx q[1];
rz(-1.2578674) q[1];
sx q[1];
rz(-1.2068611) q[1];
rz(2.76582) q[3];
sx q[3];
rz(-2.3762581) q[3];
sx q[3];
rz(0.96804141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.61509722) q[2];
sx q[2];
rz(-0.47250938) q[2];
sx q[2];
rz(-2.9769767) q[2];
rz(0.047867157) q[3];
sx q[3];
rz(-1.4048301) q[3];
sx q[3];
rz(2.3544748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98016244) q[0];
sx q[0];
rz(-0.4466559) q[0];
sx q[0];
rz(1.5572146) q[0];
rz(0.36852512) q[1];
sx q[1];
rz(-1.8199814) q[1];
sx q[1];
rz(-1.6065067) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53705789) q[0];
sx q[0];
rz(-2.4050483) q[0];
sx q[0];
rz(-0.81088541) q[0];
rz(-pi) q[1];
rz(-2.317165) q[2];
sx q[2];
rz(-1.8608837) q[2];
sx q[2];
rz(-2.1739391) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6905744) q[1];
sx q[1];
rz(-1.8201882) q[1];
sx q[1];
rz(-2.0687769) q[1];
rz(-pi) q[2];
rz(1.8593349) q[3];
sx q[3];
rz(-1.4950945) q[3];
sx q[3];
rz(0.60305616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9102455) q[2];
sx q[2];
rz(-2.7635837) q[2];
sx q[2];
rz(-2.0965915) q[2];
rz(-0.16799489) q[3];
sx q[3];
rz(-1.955227) q[3];
sx q[3];
rz(-0.35429889) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0198233) q[0];
sx q[0];
rz(-2.1124117) q[0];
sx q[0];
rz(1.9371012) q[0];
rz(-0.80351859) q[1];
sx q[1];
rz(-2.3280227) q[1];
sx q[1];
rz(-2.4258851) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7953703) q[0];
sx q[0];
rz(-2.0487794) q[0];
sx q[0];
rz(0.93670364) q[0];
rz(-pi) q[1];
x q[1];
rz(0.44797051) q[2];
sx q[2];
rz(-1.391419) q[2];
sx q[2];
rz(0.57283869) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3099311) q[1];
sx q[1];
rz(-2.4694121) q[1];
sx q[1];
rz(0.18193717) q[1];
rz(0.29362595) q[3];
sx q[3];
rz(-1.1291227) q[3];
sx q[3];
rz(0.81868086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.41010007) q[2];
sx q[2];
rz(-0.64971739) q[2];
sx q[2];
rz(-2.0984207) q[2];
rz(0.10797524) q[3];
sx q[3];
rz(-2.2004746) q[3];
sx q[3];
rz(1.1112377) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1825948) q[0];
sx q[0];
rz(-1.3866871) q[0];
sx q[0];
rz(-1.9166272) q[0];
rz(2.909868) q[1];
sx q[1];
rz(-2.2544315) q[1];
sx q[1];
rz(1.8909594) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0601025) q[0];
sx q[0];
rz(-1.4487473) q[0];
sx q[0];
rz(0.7500612) q[0];
x q[1];
rz(-1.768998) q[2];
sx q[2];
rz(-1.5491252) q[2];
sx q[2];
rz(-0.54474165) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2030484) q[1];
sx q[1];
rz(-2.7640759) q[1];
sx q[1];
rz(2.50752) q[1];
rz(-pi) q[2];
x q[2];
rz(0.55646236) q[3];
sx q[3];
rz(-1.5817405) q[3];
sx q[3];
rz(-1.8045604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0687678) q[2];
sx q[2];
rz(-2.4710957) q[2];
sx q[2];
rz(-0.13216275) q[2];
rz(-0.69563785) q[3];
sx q[3];
rz(-1.9591103) q[3];
sx q[3];
rz(-0.80719358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.950133) q[0];
sx q[0];
rz(-1.5810672) q[0];
sx q[0];
rz(-2.4323442) q[0];
rz(1.5059772) q[1];
sx q[1];
rz(-1.9874856) q[1];
sx q[1];
rz(-1.0848612) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5615254) q[0];
sx q[0];
rz(-0.97082061) q[0];
sx q[0];
rz(-1.6404433) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16895825) q[2];
sx q[2];
rz(-2.3297133) q[2];
sx q[2];
rz(0.66056992) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.9429837) q[1];
sx q[1];
rz(-1.3169022) q[1];
sx q[1];
rz(-2.6655546) q[1];
rz(-3.1301401) q[3];
sx q[3];
rz(-1.6310777) q[3];
sx q[3];
rz(-1.1520916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.52013493) q[2];
sx q[2];
rz(-1.0026714) q[2];
sx q[2];
rz(1.3168859) q[2];
rz(0.78553158) q[3];
sx q[3];
rz(-1.9436049) q[3];
sx q[3];
rz(-2.7151916) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.921628) q[0];
sx q[0];
rz(-1.6554609) q[0];
sx q[0];
rz(-2.7332136) q[0];
rz(0.95868239) q[1];
sx q[1];
rz(-2.8125693) q[1];
sx q[1];
rz(1.5214517) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.085039728) q[0];
sx q[0];
rz(-1.3360177) q[0];
sx q[0];
rz(1.9454369) q[0];
x q[1];
rz(0.076456377) q[2];
sx q[2];
rz(-0.57089134) q[2];
sx q[2];
rz(-0.61797188) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.16220763) q[1];
sx q[1];
rz(-1.9792792) q[1];
sx q[1];
rz(-0.86345203) q[1];
x q[2];
rz(-0.19100325) q[3];
sx q[3];
rz(-1.2634209) q[3];
sx q[3];
rz(1.6049811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.3756322) q[2];
sx q[2];
rz(-1.1124632) q[2];
sx q[2];
rz(-2.5755303) q[2];
rz(1.3426956) q[3];
sx q[3];
rz(-2.7261966) q[3];
sx q[3];
rz(-2.8528163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5413496) q[0];
sx q[0];
rz(-3.1025649) q[0];
sx q[0];
rz(-1.716123) q[0];
rz(-2.0478981) q[1];
sx q[1];
rz(-2.2192571) q[1];
sx q[1];
rz(0.38965449) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4132767) q[0];
sx q[0];
rz(-1.9945869) q[0];
sx q[0];
rz(-2.8601147) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8820428) q[2];
sx q[2];
rz(-0.96298157) q[2];
sx q[2];
rz(1.1961301) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.07239711) q[1];
sx q[1];
rz(-0.42617961) q[1];
sx q[1];
rz(-1.6276433) q[1];
x q[2];
rz(-1.2663307) q[3];
sx q[3];
rz(-2.048945) q[3];
sx q[3];
rz(1.696451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.23715544) q[2];
sx q[2];
rz(-0.5390141) q[2];
sx q[2];
rz(1.944444) q[2];
rz(-2.9947179) q[3];
sx q[3];
rz(-2.4683888) q[3];
sx q[3];
rz(2.1541514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(2.9523136) q[2];
sx q[2];
rz(-0.58313417) q[2];
sx q[2];
rz(-1.492576) q[2];
rz(1.3833787) q[3];
sx q[3];
rz(-1.1915177) q[3];
sx q[3];
rz(2.4080924) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
