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
rz(-2.2351216) q[0];
sx q[0];
rz(-1.6989166) q[0];
sx q[0];
rz(-0.28191167) q[0];
rz(0.52892041) q[1];
sx q[1];
rz(4.79098) q[1];
sx q[1];
rz(11.00287) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.179047) q[0];
sx q[0];
rz(-1.3907278) q[0];
sx q[0];
rz(-1.942304) q[0];
rz(-pi) q[1];
rz(-0.018997832) q[2];
sx q[2];
rz(-0.27488959) q[2];
sx q[2];
rz(-1.2541447) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8488449) q[1];
sx q[1];
rz(-2.1969921) q[1];
sx q[1];
rz(1.4153773) q[1];
rz(-pi) q[2];
rz(-0.6880349) q[3];
sx q[3];
rz(-1.8092938) q[3];
sx q[3];
rz(-1.1900589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6866744) q[2];
sx q[2];
rz(-0.29759559) q[2];
sx q[2];
rz(-2.5285524) q[2];
rz(2.6679299) q[3];
sx q[3];
rz(-1.9434171) q[3];
sx q[3];
rz(-1.4325498) q[3];
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
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4050201) q[0];
sx q[0];
rz(-2.9614145) q[0];
sx q[0];
rz(-0.75102425) q[0];
rz(0.48149064) q[1];
sx q[1];
rz(-1.0844237) q[1];
sx q[1];
rz(-2.1717333) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6365627) q[0];
sx q[0];
rz(-1.9410994) q[0];
sx q[0];
rz(-0.67815336) q[0];
x q[1];
rz(0.18377797) q[2];
sx q[2];
rz(-0.2193976) q[2];
sx q[2];
rz(0.32599005) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1456993) q[1];
sx q[1];
rz(-2.4384192) q[1];
sx q[1];
rz(1.9759167) q[1];
rz(-pi) q[2];
rz(-2.5447846) q[3];
sx q[3];
rz(-1.9442026) q[3];
sx q[3];
rz(-1.663409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.91886175) q[2];
sx q[2];
rz(-1.6657882) q[2];
sx q[2];
rz(-2.6027423) q[2];
rz(0.017596267) q[3];
sx q[3];
rz(-0.16418695) q[3];
sx q[3];
rz(0.9084475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4410412) q[0];
sx q[0];
rz(-2.7563162) q[0];
sx q[0];
rz(-3.0803296) q[0];
rz(1.5047269) q[1];
sx q[1];
rz(-0.69925362) q[1];
sx q[1];
rz(-1.9452852) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17249566) q[0];
sx q[0];
rz(-2.3845256) q[0];
sx q[0];
rz(-1.9226546) q[0];
rz(-pi) q[1];
rz(-0.2537276) q[2];
sx q[2];
rz(-2.8710033) q[2];
sx q[2];
rz(-1.21278) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.34431008) q[1];
sx q[1];
rz(-0.6629262) q[1];
sx q[1];
rz(2.1144889) q[1];
rz(-1.9798093) q[3];
sx q[3];
rz(-0.73506415) q[3];
sx q[3];
rz(0.72003698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.68941826) q[2];
sx q[2];
rz(-1.8345366) q[2];
sx q[2];
rz(-2.7090731) q[2];
rz(-2.050926) q[3];
sx q[3];
rz(-0.64041036) q[3];
sx q[3];
rz(-1.0961016) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5890305) q[0];
sx q[0];
rz(-0.93695372) q[0];
sx q[0];
rz(-2.9402148) q[0];
rz(0.51678139) q[1];
sx q[1];
rz(-0.38025451) q[1];
sx q[1];
rz(-2.1929599) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.059221642) q[0];
sx q[0];
rz(-1.9192524) q[0];
sx q[0];
rz(-0.90014622) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29251137) q[2];
sx q[2];
rz(-0.88380948) q[2];
sx q[2];
rz(-2.3623737) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.779114) q[1];
sx q[1];
rz(-2.6370905) q[1];
sx q[1];
rz(-3.1140559) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52093769) q[3];
sx q[3];
rz(-1.4864941) q[3];
sx q[3];
rz(-1.6430294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2739233) q[2];
sx q[2];
rz(-2.7322768) q[2];
sx q[2];
rz(-0.55740994) q[2];
rz(-3.0607767) q[3];
sx q[3];
rz(-1.9010474) q[3];
sx q[3];
rz(-1.2062937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4163365) q[0];
sx q[0];
rz(-1.6660322) q[0];
sx q[0];
rz(-2.1112554) q[0];
rz(-0.15580767) q[1];
sx q[1];
rz(-2.0741597) q[1];
sx q[1];
rz(0.20419289) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0263173) q[0];
sx q[0];
rz(-1.6395307) q[0];
sx q[0];
rz(2.8990549) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5956018) q[2];
sx q[2];
rz(-0.24790774) q[2];
sx q[2];
rz(1.1746097) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9928007) q[1];
sx q[1];
rz(-0.53314236) q[1];
sx q[1];
rz(-0.41457446) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.810012) q[3];
sx q[3];
rz(-0.81221928) q[3];
sx q[3];
rz(1.9977457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.92945176) q[2];
sx q[2];
rz(-1.6484304) q[2];
sx q[2];
rz(-0.41352752) q[2];
rz(0.84189576) q[3];
sx q[3];
rz(-1.3827518) q[3];
sx q[3];
rz(0.97257417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0464762) q[0];
sx q[0];
rz(-2.0501417) q[0];
sx q[0];
rz(0.0055775642) q[0];
rz(-1.7976409) q[1];
sx q[1];
rz(-0.19897142) q[1];
sx q[1];
rz(0.70095789) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0049135) q[0];
sx q[0];
rz(-2.4304996) q[0];
sx q[0];
rz(-2.2836766) q[0];
rz(-pi) q[1];
rz(2.7024038) q[2];
sx q[2];
rz(-2.1648295) q[2];
sx q[2];
rz(-2.8043945) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8474104) q[1];
sx q[1];
rz(-2.8632616) q[1];
sx q[1];
rz(0.14974447) q[1];
rz(2.088955) q[3];
sx q[3];
rz(-2.0817882) q[3];
sx q[3];
rz(0.68447733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.31412101) q[2];
sx q[2];
rz(-2.8077329) q[2];
sx q[2];
rz(1.9742879) q[2];
rz(0.88268924) q[3];
sx q[3];
rz(-2.3736931) q[3];
sx q[3];
rz(0.32514969) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90671396) q[0];
sx q[0];
rz(-1.119708) q[0];
sx q[0];
rz(3.0562905) q[0];
rz(2.3872497) q[1];
sx q[1];
rz(-2.6383548) q[1];
sx q[1];
rz(1.0275966) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9011544) q[0];
sx q[0];
rz(-1.580997) q[0];
sx q[0];
rz(-2.9290694) q[0];
x q[1];
rz(1.5761887) q[2];
sx q[2];
rz(-2.1674334) q[2];
sx q[2];
rz(-0.83534681) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.8366521) q[1];
sx q[1];
rz(-1.8936367) q[1];
sx q[1];
rz(3.0129827) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47691741) q[3];
sx q[3];
rz(-2.0774042) q[3];
sx q[3];
rz(1.5445386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9895642) q[2];
sx q[2];
rz(-0.98735183) q[2];
sx q[2];
rz(0.41934553) q[2];
rz(2.5069405) q[3];
sx q[3];
rz(-2.3366163) q[3];
sx q[3];
rz(2.0161207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.80970508) q[0];
sx q[0];
rz(-2.2655847) q[0];
sx q[0];
rz(2.4494655) q[0];
rz(-2.5634815) q[1];
sx q[1];
rz(-0.41936857) q[1];
sx q[1];
rz(-1.7587761) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.064474) q[0];
sx q[0];
rz(-3.0713284) q[0];
sx q[0];
rz(-1.2901359) q[0];
rz(-pi) q[1];
rz(-0.77246364) q[2];
sx q[2];
rz(-0.83784762) q[2];
sx q[2];
rz(-0.39847429) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9585413) q[1];
sx q[1];
rz(-1.019283) q[1];
sx q[1];
rz(-2.8666878) q[1];
x q[2];
rz(2.6666497) q[3];
sx q[3];
rz(-0.90427665) q[3];
sx q[3];
rz(-0.24196821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.77070037) q[2];
sx q[2];
rz(-2.6110677) q[2];
sx q[2];
rz(1.4250866) q[2];
rz(0.51618451) q[3];
sx q[3];
rz(-2.1631212) q[3];
sx q[3];
rz(-1.2396575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4952963) q[0];
sx q[0];
rz(-1.4343867) q[0];
sx q[0];
rz(-2.7700951) q[0];
rz(2.3067572) q[1];
sx q[1];
rz(-0.8131665) q[1];
sx q[1];
rz(0.017597839) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4478701) q[0];
sx q[0];
rz(-1.2945064) q[0];
sx q[0];
rz(-2.1135065) q[0];
x q[1];
rz(1.8210635) q[2];
sx q[2];
rz(-2.2907718) q[2];
sx q[2];
rz(-1.7257476) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0549893) q[1];
sx q[1];
rz(-1.9412244) q[1];
sx q[1];
rz(0.61062529) q[1];
rz(1.2491717) q[3];
sx q[3];
rz(-0.98708803) q[3];
sx q[3];
rz(-0.58313384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1060433) q[2];
sx q[2];
rz(-2.273166) q[2];
sx q[2];
rz(-3.0412728) q[2];
rz(0.22942461) q[3];
sx q[3];
rz(-2.094163) q[3];
sx q[3];
rz(2.6917698) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4556274) q[0];
sx q[0];
rz(-2.8534511) q[0];
sx q[0];
rz(-0.57383865) q[0];
rz(-1.0539184) q[1];
sx q[1];
rz(-2.0951447) q[1];
sx q[1];
rz(-0.67646772) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5887816) q[0];
sx q[0];
rz(-1.8491321) q[0];
sx q[0];
rz(-1.3180483) q[0];
x q[1];
rz(-1.9755115) q[2];
sx q[2];
rz(-1.3339443) q[2];
sx q[2];
rz(-0.39654695) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0645731) q[1];
sx q[1];
rz(-1.6592245) q[1];
sx q[1];
rz(-0.52845533) q[1];
x q[2];
rz(-2.1031688) q[3];
sx q[3];
rz(-1.7422973) q[3];
sx q[3];
rz(2.4166783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.68710589) q[2];
sx q[2];
rz(-2.3837619) q[2];
sx q[2];
rz(1.494361) q[2];
rz(2.2729661) q[3];
sx q[3];
rz(-2.4354911) q[3];
sx q[3];
rz(-0.47006616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0710707) q[0];
sx q[0];
rz(-1.949953) q[0];
sx q[0];
rz(-1.9198887) q[0];
rz(1.8078049) q[1];
sx q[1];
rz(-1.8232657) q[1];
sx q[1];
rz(2.5008536) q[1];
rz(0.76539466) q[2];
sx q[2];
rz(-1.5049878) q[2];
sx q[2];
rz(2.4743248) q[2];
rz(0.044338772) q[3];
sx q[3];
rz(-1.7462742) q[3];
sx q[3];
rz(-2.3901226) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
