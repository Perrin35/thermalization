OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5698009) q[0];
sx q[0];
rz(4.0507841) q[0];
sx q[0];
rz(7.9180766) q[0];
rz(1.2070967) q[1];
sx q[1];
rz(-2.6488882) q[1];
sx q[1];
rz(-0.55941137) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6772422) q[0];
sx q[0];
rz(-1.4046245) q[0];
sx q[0];
rz(-2.8234923) q[0];
rz(-pi) q[1];
rz(0.71601358) q[2];
sx q[2];
rz(-0.56173827) q[2];
sx q[2];
rz(-1.5519976) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.2845498) q[1];
sx q[1];
rz(-0.27978292) q[1];
sx q[1];
rz(-1.5955052) q[1];
x q[2];
rz(2.6443566) q[3];
sx q[3];
rz(-2.1057893) q[3];
sx q[3];
rz(2.3930156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.59014615) q[2];
sx q[2];
rz(-0.85942736) q[2];
sx q[2];
rz(2.0476511) q[2];
rz(1.9048196) q[3];
sx q[3];
rz(-1.1636795) q[3];
sx q[3];
rz(2.8804603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.3742438) q[0];
sx q[0];
rz(-1.2232895) q[0];
sx q[0];
rz(-2.4236524) q[0];
rz(-1.9473437) q[1];
sx q[1];
rz(-0.4078882) q[1];
sx q[1];
rz(1.6494707) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46280086) q[0];
sx q[0];
rz(-2.959842) q[0];
sx q[0];
rz(2.3857855) q[0];
rz(-pi) q[1];
rz(0.43955438) q[2];
sx q[2];
rz(-1.2258523) q[2];
sx q[2];
rz(-2.3203691) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.9886943) q[1];
sx q[1];
rz(-1.1225268) q[1];
sx q[1];
rz(0.50181915) q[1];
x q[2];
rz(2.9255052) q[3];
sx q[3];
rz(-0.83744811) q[3];
sx q[3];
rz(2.2186389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0481999) q[2];
sx q[2];
rz(-0.77646774) q[2];
sx q[2];
rz(2.7788567) q[2];
rz(2.2705966) q[3];
sx q[3];
rz(-1.3716776) q[3];
sx q[3];
rz(-1.3236275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094630346) q[0];
sx q[0];
rz(-1.8383263) q[0];
sx q[0];
rz(-2.5808425) q[0];
rz(-2.1669855) q[1];
sx q[1];
rz(-2.2798996) q[1];
sx q[1];
rz(0.21751705) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64338486) q[0];
sx q[0];
rz(-1.5899993) q[0];
sx q[0];
rz(-0.58450825) q[0];
rz(-1.0571805) q[2];
sx q[2];
rz(-1.4808146) q[2];
sx q[2];
rz(-2.1291901) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.88918658) q[1];
sx q[1];
rz(-1.6853309) q[1];
sx q[1];
rz(1.8327528) q[1];
rz(-pi) q[2];
rz(-2.5325696) q[3];
sx q[3];
rz(-1.6224738) q[3];
sx q[3];
rz(0.9059815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4983623) q[2];
sx q[2];
rz(-2.3537894) q[2];
sx q[2];
rz(2.1882449) q[2];
rz(-1.0777473) q[3];
sx q[3];
rz(-1.4606303) q[3];
sx q[3];
rz(-0.47422153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
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
rz(1.0039907) q[0];
sx q[0];
rz(-0.1739665) q[0];
sx q[0];
rz(2.2250788) q[0];
rz(-2.7169531) q[1];
sx q[1];
rz(-1.3820796) q[1];
sx q[1];
rz(2.0133846) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5461499) q[0];
sx q[0];
rz(-0.91888035) q[0];
sx q[0];
rz(2.2007887) q[0];
x q[1];
rz(1.3638956) q[2];
sx q[2];
rz(-0.82447663) q[2];
sx q[2];
rz(1.7306223) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.27552106) q[1];
sx q[1];
rz(-1.0539319) q[1];
sx q[1];
rz(-1.2244774) q[1];
rz(-pi) q[2];
rz(0.22596328) q[3];
sx q[3];
rz(-0.98528242) q[3];
sx q[3];
rz(2.8087392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.26587036) q[2];
sx q[2];
rz(-1.46572) q[2];
sx q[2];
rz(1.1992559) q[2];
rz(1.1572329) q[3];
sx q[3];
rz(-0.80409378) q[3];
sx q[3];
rz(-0.51586241) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1542094) q[0];
sx q[0];
rz(-3.0976384) q[0];
sx q[0];
rz(-0.92535812) q[0];
rz(0.56272733) q[1];
sx q[1];
rz(-2.1268763) q[1];
sx q[1];
rz(-2.7664807) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1031418) q[0];
sx q[0];
rz(-1.4618357) q[0];
sx q[0];
rz(0.54619638) q[0];
x q[1];
rz(3.0882224) q[2];
sx q[2];
rz(-2.0056021) q[2];
sx q[2];
rz(-2.82719) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.28175345) q[1];
sx q[1];
rz(-1.5956164) q[1];
sx q[1];
rz(0.16902216) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6869223) q[3];
sx q[3];
rz(-1.934762) q[3];
sx q[3];
rz(0.92586799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.23842266) q[2];
sx q[2];
rz(-2.2245202) q[2];
sx q[2];
rz(-0.49752107) q[2];
rz(-2.7072952) q[3];
sx q[3];
rz(-1.6853761) q[3];
sx q[3];
rz(0.98880497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.3797746) q[0];
sx q[0];
rz(-0.86123818) q[0];
sx q[0];
rz(-2.8616943) q[0];
rz(2.1474536) q[1];
sx q[1];
rz(-0.86052624) q[1];
sx q[1];
rz(-2.5631189) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6743603) q[0];
sx q[0];
rz(-1.064309) q[0];
sx q[0];
rz(-1.1367528) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4388588) q[2];
sx q[2];
rz(-1.8763855) q[2];
sx q[2];
rz(-1.0440799) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7565365) q[1];
sx q[1];
rz(-2.2900683) q[1];
sx q[1];
rz(1.8315218) q[1];
rz(-pi) q[2];
rz(-1.6719477) q[3];
sx q[3];
rz(-1.3037852) q[3];
sx q[3];
rz(1.284541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0591187) q[2];
sx q[2];
rz(-2.5422577) q[2];
sx q[2];
rz(-2.2590051) q[2];
rz(-2.5146218) q[3];
sx q[3];
rz(-0.65492237) q[3];
sx q[3];
rz(-2.2466808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51666981) q[0];
sx q[0];
rz(-2.1997917) q[0];
sx q[0];
rz(-0.00061568419) q[0];
rz(-2.2387538) q[1];
sx q[1];
rz(-0.87264624) q[1];
sx q[1];
rz(3.0965064) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49360291) q[0];
sx q[0];
rz(-1.6479074) q[0];
sx q[0];
rz(-0.48992975) q[0];
rz(-pi) q[1];
rz(-0.22901625) q[2];
sx q[2];
rz(-0.75521246) q[2];
sx q[2];
rz(1.7984185) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5645091) q[1];
sx q[1];
rz(-2.2840743) q[1];
sx q[1];
rz(1.1611931) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5302913) q[3];
sx q[3];
rz(-1.0795648) q[3];
sx q[3];
rz(-2.4424551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8467329) q[2];
sx q[2];
rz(-0.555987) q[2];
sx q[2];
rz(2.2247458) q[2];
rz(0.89546853) q[3];
sx q[3];
rz(-1.5119036) q[3];
sx q[3];
rz(-1.2112613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.004892) q[0];
sx q[0];
rz(-3.0558375) q[0];
sx q[0];
rz(0.46491796) q[0];
rz(1.9526941) q[1];
sx q[1];
rz(-1.3465954) q[1];
sx q[1];
rz(-0.37110034) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8086149) q[0];
sx q[0];
rz(-0.59691256) q[0];
sx q[0];
rz(0.78947584) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6109054) q[2];
sx q[2];
rz(-2.3903987) q[2];
sx q[2];
rz(1.7263279) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7830989) q[1];
sx q[1];
rz(-0.91121735) q[1];
sx q[1];
rz(-2.5765319) q[1];
x q[2];
rz(-0.38406541) q[3];
sx q[3];
rz(-1.3017941) q[3];
sx q[3];
rz(-0.88685461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4798639) q[2];
sx q[2];
rz(-0.78833818) q[2];
sx q[2];
rz(0.41024497) q[2];
rz(-1.5069626) q[3];
sx q[3];
rz(-2.7362636) q[3];
sx q[3];
rz(-3.098367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0931382) q[0];
sx q[0];
rz(-0.52670902) q[0];
sx q[0];
rz(-0.17573892) q[0];
rz(-2.3161092) q[1];
sx q[1];
rz(-1.8800507) q[1];
sx q[1];
rz(0.82297355) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1639481) q[0];
sx q[0];
rz(-0.86714449) q[0];
sx q[0];
rz(2.4022341) q[0];
rz(0.31655471) q[2];
sx q[2];
rz(-1.9783894) q[2];
sx q[2];
rz(2.0616693) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.6187669) q[1];
sx q[1];
rz(-0.78729388) q[1];
sx q[1];
rz(-2.6403422) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4863344) q[3];
sx q[3];
rz(-0.41197398) q[3];
sx q[3];
rz(-0.99124817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2814111) q[2];
sx q[2];
rz(-1.9965636) q[2];
sx q[2];
rz(0.54110503) q[2];
rz(2.7152854) q[3];
sx q[3];
rz(-2.6796894) q[3];
sx q[3];
rz(-2.2947252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.68778872) q[0];
sx q[0];
rz(-2.332088) q[0];
sx q[0];
rz(0.33568207) q[0];
rz(3.1188534) q[1];
sx q[1];
rz(-0.72240654) q[1];
sx q[1];
rz(2.4679599) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.861872) q[0];
sx q[0];
rz(-0.69929823) q[0];
sx q[0];
rz(2.5847816) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6725476) q[2];
sx q[2];
rz(-1.4802336) q[2];
sx q[2];
rz(-3.0770214) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2362263) q[1];
sx q[1];
rz(-2.1520808) q[1];
sx q[1];
rz(-0.24317113) q[1];
rz(0.38506759) q[3];
sx q[3];
rz(-0.35849471) q[3];
sx q[3];
rz(0.26824238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6966887) q[2];
sx q[2];
rz(-0.99385571) q[2];
sx q[2];
rz(-1.0059086) q[2];
rz(0.78399793) q[3];
sx q[3];
rz(-0.32117143) q[3];
sx q[3];
rz(1.4348437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.0189331) q[0];
sx q[0];
rz(-1.7187332) q[0];
sx q[0];
rz(-1.1359489) q[0];
rz(0.38416531) q[1];
sx q[1];
rz(-1.1687678) q[1];
sx q[1];
rz(1.6484177) q[1];
rz(0.21239077) q[2];
sx q[2];
rz(-1.7707386) q[2];
sx q[2];
rz(2.7223827) q[2];
rz(0.85687153) q[3];
sx q[3];
rz(-1.8663434) q[3];
sx q[3];
rz(0.28874884) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
