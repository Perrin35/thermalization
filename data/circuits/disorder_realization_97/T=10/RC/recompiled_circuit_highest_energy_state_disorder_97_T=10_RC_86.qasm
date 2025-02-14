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
rz(-1.462734) q[1];
sx q[1];
rz(-0.94574133) q[1];
sx q[1];
rz(2.6374964) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60668463) q[0];
sx q[0];
rz(-2.1143171) q[0];
sx q[0];
rz(0.96291079) q[0];
rz(-pi) q[1];
rz(-2.7934876) q[2];
sx q[2];
rz(-1.7278302) q[2];
sx q[2];
rz(0.43660313) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.13611804) q[1];
sx q[1];
rz(-2.1393993) q[1];
sx q[1];
rz(2.3852045) q[1];
rz(0.33511843) q[3];
sx q[3];
rz(-1.6168645) q[3];
sx q[3];
rz(-0.81460458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.41245875) q[2];
sx q[2];
rz(-1.8987741) q[2];
sx q[2];
rz(-2.9553555) q[2];
rz(1.3650182) q[3];
sx q[3];
rz(-0.61187196) q[3];
sx q[3];
rz(1.9369102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3314826) q[0];
sx q[0];
rz(-0.45833603) q[0];
sx q[0];
rz(-0.052074281) q[0];
rz(2.1049818) q[1];
sx q[1];
rz(-2.5317445) q[1];
sx q[1];
rz(2.5263272) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4724169) q[0];
sx q[0];
rz(-2.376498) q[0];
sx q[0];
rz(-1.0932176) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6411166) q[2];
sx q[2];
rz(-2.8722974) q[2];
sx q[2];
rz(-0.057502086) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.8456455) q[1];
sx q[1];
rz(-1.5145647) q[1];
sx q[1];
rz(2.2241431) q[1];
x q[2];
rz(-0.21188696) q[3];
sx q[3];
rz(-2.35244) q[3];
sx q[3];
rz(-1.0494389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.692824) q[2];
sx q[2];
rz(-1.8961467) q[2];
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
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6579984) q[0];
sx q[0];
rz(-2.3880385) q[0];
sx q[0];
rz(-1.765522) q[0];
rz(0.54028571) q[1];
sx q[1];
rz(-1.0397725) q[1];
sx q[1];
rz(-1.4556063) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8309915) q[0];
sx q[0];
rz(-1.5001703) q[0];
sx q[0];
rz(-2.7381606) q[0];
rz(2.7305538) q[2];
sx q[2];
rz(-0.22145311) q[2];
sx q[2];
rz(-2.7483181) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4306304) q[1];
sx q[1];
rz(-0.73657958) q[1];
sx q[1];
rz(2.1107091) q[1];
x q[2];
rz(2.2690587) q[3];
sx q[3];
rz(-1.8909406) q[3];
sx q[3];
rz(0.24139377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.49440631) q[2];
sx q[2];
rz(-0.33438412) q[2];
sx q[2];
rz(1.6953267) q[2];
rz(1.9485731) q[3];
sx q[3];
rz(-2.1665067) q[3];
sx q[3];
rz(1.4421991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0541662) q[0];
sx q[0];
rz(-2.8717201) q[0];
sx q[0];
rz(-2.7179981) q[0];
rz(2.1845747) q[1];
sx q[1];
rz(-1.7901763) q[1];
sx q[1];
rz(-0.097600309) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1621057) q[0];
sx q[0];
rz(-0.66363664) q[0];
sx q[0];
rz(2.5121793) q[0];
rz(-pi) q[1];
rz(-2.9984442) q[2];
sx q[2];
rz(-0.88298702) q[2];
sx q[2];
rz(2.0542415) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7047119) q[1];
sx q[1];
rz(-1.8837253) q[1];
sx q[1];
rz(-1.2068611) q[1];
rz(-2.76582) q[3];
sx q[3];
rz(-2.3762581) q[3];
sx q[3];
rz(-0.96804141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5264954) q[2];
sx q[2];
rz(-0.47250938) q[2];
sx q[2];
rz(-0.16461593) q[2];
rz(-0.047867157) q[3];
sx q[3];
rz(-1.4048301) q[3];
sx q[3];
rz(0.78711787) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1614302) q[0];
sx q[0];
rz(-0.4466559) q[0];
sx q[0];
rz(-1.584378) q[0];
rz(0.36852512) q[1];
sx q[1];
rz(-1.8199814) q[1];
sx q[1];
rz(-1.6065067) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4457065) q[0];
sx q[0];
rz(-2.0793756) q[0];
sx q[0];
rz(-0.558338) q[0];
rz(1.1565882) q[2];
sx q[2];
rz(-2.3510859) q[2];
sx q[2];
rz(-0.3031916) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.69362946) q[1];
sx q[1];
rz(-2.5893859) q[1];
sx q[1];
rz(-2.0606661) q[1];
rz(-pi) q[2];
rz(1.8312884) q[3];
sx q[3];
rz(-0.29803571) q[3];
sx q[3];
rz(2.4233415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9102455) q[2];
sx q[2];
rz(-2.7635837) q[2];
sx q[2];
rz(1.0450012) q[2];
rz(-0.16799489) q[3];
sx q[3];
rz(-1.955227) q[3];
sx q[3];
rz(2.7872938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
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
rz(-1.2044915) q[0];
rz(-0.80351859) q[1];
sx q[1];
rz(-2.3280227) q[1];
sx q[1];
rz(0.71570754) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2432118) q[0];
sx q[0];
rz(-2.1247433) q[0];
sx q[0];
rz(-2.5701163) q[0];
x q[1];
rz(-1.3722695) q[2];
sx q[2];
rz(-2.0110705) q[2];
sx q[2];
rz(2.0581051) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5454272) q[1];
sx q[1];
rz(-1.6837032) q[1];
sx q[1];
rz(0.66408709) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1200772) q[3];
sx q[3];
rz(-0.52495723) q[3];
sx q[3];
rz(0.20302816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.41010007) q[2];
sx q[2];
rz(-2.4918753) q[2];
sx q[2];
rz(-2.0984207) q[2];
rz(0.10797524) q[3];
sx q[3];
rz(-0.94111809) q[3];
sx q[3];
rz(-1.1112377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.9589979) q[0];
sx q[0];
rz(-1.7549055) q[0];
sx q[0];
rz(-1.2249655) q[0];
rz(0.23172465) q[1];
sx q[1];
rz(-2.2544315) q[1];
sx q[1];
rz(-1.8909594) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0601025) q[0];
sx q[0];
rz(-1.6928453) q[0];
sx q[0];
rz(2.3915315) q[0];
x q[1];
rz(1.3725946) q[2];
sx q[2];
rz(-1.5924675) q[2];
sx q[2];
rz(-2.596851) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.96737427) q[1];
sx q[1];
rz(-1.7909482) q[1];
sx q[1];
rz(-0.30921127) q[1];
x q[2];
rz(3.1208745) q[3];
sx q[3];
rz(-0.55655863) q[3];
sx q[3];
rz(-0.25135612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0687678) q[2];
sx q[2];
rz(-0.67049694) q[2];
sx q[2];
rz(3.0094299) q[2];
rz(-0.69563785) q[3];
sx q[3];
rz(-1.9591103) q[3];
sx q[3];
rz(-0.80719358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.950133) q[0];
sx q[0];
rz(-1.5810672) q[0];
sx q[0];
rz(0.70924846) q[0];
rz(1.5059772) q[1];
sx q[1];
rz(-1.154107) q[1];
sx q[1];
rz(1.0848612) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4385972) q[0];
sx q[0];
rz(-2.5380822) q[0];
sx q[0];
rz(0.10137697) q[0];
rz(-pi) q[1];
rz(0.80471595) q[2];
sx q[2];
rz(-1.6931117) q[2];
sx q[2];
rz(-2.1145156) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0604531) q[1];
sx q[1];
rz(-2.6067002) q[1];
sx q[1];
rz(-2.6263531) q[1];
rz(-pi) q[2];
rz(1.6310817) q[3];
sx q[3];
rz(-1.5593646) q[3];
sx q[3];
rz(0.41939467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6214577) q[2];
sx q[2];
rz(-1.0026714) q[2];
sx q[2];
rz(-1.3168859) q[2];
rz(0.78553158) q[3];
sx q[3];
rz(-1.1979878) q[3];
sx q[3];
rz(2.7151916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21996466) q[0];
sx q[0];
rz(-1.4861318) q[0];
sx q[0];
rz(0.40837902) q[0];
rz(-0.95868239) q[1];
sx q[1];
rz(-2.8125693) q[1];
sx q[1];
rz(-1.5214517) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1216089) q[0];
sx q[0];
rz(-0.43916288) q[0];
sx q[0];
rz(0.99187054) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6198115) q[2];
sx q[2];
rz(-2.1398126) q[2];
sx q[2];
rz(-0.70876497) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.16220763) q[1];
sx q[1];
rz(-1.1623135) q[1];
sx q[1];
rz(-0.86345203) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9505894) q[3];
sx q[3];
rz(-1.2634209) q[3];
sx q[3];
rz(1.5366116) q[3];
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
rz(-1.3426956) q[3];
sx q[3];
rz(-2.7261966) q[3];
sx q[3];
rz(2.8528163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.5413496) q[0];
sx q[0];
rz(-3.1025649) q[0];
sx q[0];
rz(1.4254697) q[0];
rz(1.0936945) q[1];
sx q[1];
rz(-2.2192571) q[1];
sx q[1];
rz(-2.7519382) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.728316) q[0];
sx q[0];
rz(-1.9945869) q[0];
sx q[0];
rz(-2.8601147) q[0];
rz(-pi) q[1];
x q[1];
rz(0.25954982) q[2];
sx q[2];
rz(-2.1786111) q[2];
sx q[2];
rz(-1.1961301) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0691955) q[1];
sx q[1];
rz(-0.42617961) q[1];
sx q[1];
rz(-1.5139493) q[1];
rz(-pi) q[2];
rz(0.52441729) q[3];
sx q[3];
rz(-2.5811385) q[3];
sx q[3];
rz(0.84597142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.23715544) q[2];
sx q[2];
rz(-0.5390141) q[2];
sx q[2];
rz(1.1971486) q[2];
rz(-0.14687471) q[3];
sx q[3];
rz(-2.4683888) q[3];
sx q[3];
rz(0.98744121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7214397) q[0];
sx q[0];
rz(-2.1099821) q[0];
sx q[0];
rz(-1.4954062) q[0];
rz(0.40311665) q[1];
sx q[1];
rz(-1.8164201) q[1];
sx q[1];
rz(1.4774189) q[1];
rz(-2.9523136) q[2];
sx q[2];
rz(-2.5584585) q[2];
sx q[2];
rz(1.6490166) q[2];
rz(-2.7561989) q[3];
sx q[3];
rz(-1.3968395) q[3];
sx q[3];
rz(-2.2341961) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
