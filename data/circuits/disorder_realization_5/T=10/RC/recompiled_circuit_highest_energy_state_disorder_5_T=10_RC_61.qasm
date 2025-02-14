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
rz(1.16601) q[0];
sx q[0];
rz(-0.87870413) q[0];
sx q[0];
rz(0.19600828) q[0];
rz(-1.4930383) q[1];
sx q[1];
rz(-0.9382481) q[1];
sx q[1];
rz(0.87747639) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2169133) q[0];
sx q[0];
rz(-1.3558421) q[0];
sx q[0];
rz(2.0314905) q[0];
rz(-pi) q[1];
rz(2.4634667) q[2];
sx q[2];
rz(-2.6161751) q[2];
sx q[2];
rz(-0.12066387) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0044549) q[1];
sx q[1];
rz(-2.6124705) q[1];
sx q[1];
rz(-1.0164539) q[1];
x q[2];
rz(-2.0480754) q[3];
sx q[3];
rz(-1.358489) q[3];
sx q[3];
rz(-2.4119854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.33285546) q[2];
sx q[2];
rz(-0.70781195) q[2];
sx q[2];
rz(-0.037192496) q[2];
rz(-0.91341364) q[3];
sx q[3];
rz(-0.98228684) q[3];
sx q[3];
rz(1.6325525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7971802) q[0];
sx q[0];
rz(-2.0672815) q[0];
sx q[0];
rz(-0.2221701) q[0];
rz(0.19368681) q[1];
sx q[1];
rz(-1.2682468) q[1];
sx q[1];
rz(2.8817315) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65529275) q[0];
sx q[0];
rz(-1.813714) q[0];
sx q[0];
rz(-1.1134558) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.1078306) q[2];
sx q[2];
rz(-1.3930151) q[2];
sx q[2];
rz(-1.1441292) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4603468) q[1];
sx q[1];
rz(-2.859194) q[1];
sx q[1];
rz(0.71747924) q[1];
x q[2];
rz(-1.0845029) q[3];
sx q[3];
rz(-2.2688365) q[3];
sx q[3];
rz(1.2459038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3746609) q[2];
sx q[2];
rz(-2.5587475) q[2];
sx q[2];
rz(2.984821) q[2];
rz(0.80592704) q[3];
sx q[3];
rz(-1.0904049) q[3];
sx q[3];
rz(0.52483112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90461516) q[0];
sx q[0];
rz(-0.1122864) q[0];
sx q[0];
rz(1.9554546) q[0];
rz(-0.063749464) q[1];
sx q[1];
rz(-1.5260162) q[1];
sx q[1];
rz(0.41248163) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97226364) q[0];
sx q[0];
rz(-2.0893851) q[0];
sx q[0];
rz(-0.56904582) q[0];
rz(-pi) q[1];
rz(-1.3774728) q[2];
sx q[2];
rz(-0.68916048) q[2];
sx q[2];
rz(-1.2528407) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3220249) q[1];
sx q[1];
rz(-1.7734999) q[1];
sx q[1];
rz(-2.8611819) q[1];
rz(-pi) q[2];
rz(0.58150141) q[3];
sx q[3];
rz(-0.99101725) q[3];
sx q[3];
rz(1.2691154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5731262) q[2];
sx q[2];
rz(-0.44034475) q[2];
sx q[2];
rz(-1.9754515) q[2];
rz(2.3447573) q[3];
sx q[3];
rz(-1.7723869) q[3];
sx q[3];
rz(-0.78127512) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2260988) q[0];
sx q[0];
rz(-1.508536) q[0];
sx q[0];
rz(-1.728212) q[0];
rz(0.49890292) q[1];
sx q[1];
rz(-1.3830802) q[1];
sx q[1];
rz(-1.2938719) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2436566) q[0];
sx q[0];
rz(-1.5439646) q[0];
sx q[0];
rz(-1.4856337) q[0];
x q[1];
rz(1.7912553) q[2];
sx q[2];
rz(-0.33500698) q[2];
sx q[2];
rz(-2.9398244) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9905967) q[1];
sx q[1];
rz(-1.5676204) q[1];
sx q[1];
rz(3.0795829) q[1];
x q[2];
rz(1.5086725) q[3];
sx q[3];
rz(-1.5034564) q[3];
sx q[3];
rz(1.289866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6125907) q[2];
sx q[2];
rz(-2.4606073) q[2];
sx q[2];
rz(0.16981086) q[2];
rz(2.7134907) q[3];
sx q[3];
rz(-1.7120818) q[3];
sx q[3];
rz(-0.16038173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.661628) q[0];
sx q[0];
rz(-2.7492838) q[0];
sx q[0];
rz(-2.3062134) q[0];
rz(0.63940489) q[1];
sx q[1];
rz(-2.4765922) q[1];
sx q[1];
rz(-2.0506052) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3314455) q[0];
sx q[0];
rz(-1.5584599) q[0];
sx q[0];
rz(2.4452565) q[0];
rz(2.8418737) q[2];
sx q[2];
rz(-1.7309566) q[2];
sx q[2];
rz(1.0417787) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.94746339) q[1];
sx q[1];
rz(-0.76912472) q[1];
sx q[1];
rz(-0.98712732) q[1];
rz(-2.1017081) q[3];
sx q[3];
rz(-0.58947488) q[3];
sx q[3];
rz(0.34661759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.368448) q[2];
sx q[2];
rz(-1.1870563) q[2];
sx q[2];
rz(3.1375569) q[2];
rz(2.0648547) q[3];
sx q[3];
rz(-2.4439947) q[3];
sx q[3];
rz(-2.9546886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3183597) q[0];
sx q[0];
rz(-1.3613181) q[0];
sx q[0];
rz(2.5079492) q[0];
rz(2.7404495) q[1];
sx q[1];
rz(-2.432343) q[1];
sx q[1];
rz(2.4493682) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5244103) q[0];
sx q[0];
rz(-1.1600736) q[0];
sx q[0];
rz(2.3401712) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2622617) q[2];
sx q[2];
rz(-0.82595982) q[2];
sx q[2];
rz(0.334563) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0226188) q[1];
sx q[1];
rz(-1.8198208) q[1];
sx q[1];
rz(2.157446) q[1];
x q[2];
rz(-0.68128392) q[3];
sx q[3];
rz(-2.9046106) q[3];
sx q[3];
rz(2.9221688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9358518) q[2];
sx q[2];
rz(-0.78448272) q[2];
sx q[2];
rz(-1.9007696) q[2];
rz(1.9567018) q[3];
sx q[3];
rz(-1.840206) q[3];
sx q[3];
rz(1.544781) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.009509) q[0];
sx q[0];
rz(-1.3030095) q[0];
sx q[0];
rz(-1.3707772) q[0];
rz(-0.41807434) q[1];
sx q[1];
rz(-1.6590174) q[1];
sx q[1];
rz(-0.48666993) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3330227) q[0];
sx q[0];
rz(-1.4736703) q[0];
sx q[0];
rz(2.1311231) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1407632) q[2];
sx q[2];
rz(-2.5116133) q[2];
sx q[2];
rz(0.85106817) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.44341893) q[1];
sx q[1];
rz(-1.3935745) q[1];
sx q[1];
rz(-1.366057) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6177931) q[3];
sx q[3];
rz(-1.1011657) q[3];
sx q[3];
rz(-2.3835973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.70249867) q[2];
sx q[2];
rz(-1.0489901) q[2];
sx q[2];
rz(0.19719633) q[2];
rz(-2.6321453) q[3];
sx q[3];
rz(-2.8809437) q[3];
sx q[3];
rz(2.313224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2800804) q[0];
sx q[0];
rz(-2.0674288) q[0];
sx q[0];
rz(1.7751088) q[0];
rz(-0.81661433) q[1];
sx q[1];
rz(-1.0895224) q[1];
sx q[1];
rz(-2.8588967) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6940281) q[0];
sx q[0];
rz(-0.97937939) q[0];
sx q[0];
rz(-2.924218) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.28699283) q[2];
sx q[2];
rz(-0.32397917) q[2];
sx q[2];
rz(-0.92008495) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7636488) q[1];
sx q[1];
rz(-0.70893439) q[1];
sx q[1];
rz(1.8590742) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52146179) q[3];
sx q[3];
rz(-1.9127426) q[3];
sx q[3];
rz(-1.6638883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1064328) q[2];
sx q[2];
rz(-2.1356434) q[2];
sx q[2];
rz(2.1584568) q[2];
rz(1.2784917) q[3];
sx q[3];
rz(-0.92517868) q[3];
sx q[3];
rz(2.2127051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4814257) q[0];
sx q[0];
rz(-1.3496512) q[0];
sx q[0];
rz(-0.31326681) q[0];
rz(-0.52472862) q[1];
sx q[1];
rz(-0.26291651) q[1];
sx q[1];
rz(-1.8892586) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1850644) q[0];
sx q[0];
rz(-1.4502429) q[0];
sx q[0];
rz(-0.53114364) q[0];
x q[1];
rz(-2.5400794) q[2];
sx q[2];
rz(-0.85962112) q[2];
sx q[2];
rz(1.892923) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.59555) q[1];
sx q[1];
rz(-2.3417406) q[1];
sx q[1];
rz(-2.355666) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81635565) q[3];
sx q[3];
rz(-1.6301148) q[3];
sx q[3];
rz(0.99494464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.19920443) q[2];
sx q[2];
rz(-0.62493268) q[2];
sx q[2];
rz(-1.7601298) q[2];
rz(0.30113164) q[3];
sx q[3];
rz(-2.1589409) q[3];
sx q[3];
rz(0.38100955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9075266) q[0];
sx q[0];
rz(-2.1387687) q[0];
sx q[0];
rz(-1.7940849) q[0];
rz(-2.2566336) q[1];
sx q[1];
rz(-0.62565175) q[1];
sx q[1];
rz(-1.7431097) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2725582) q[0];
sx q[0];
rz(-1.5814651) q[0];
sx q[0];
rz(1.19633) q[0];
rz(-pi) q[1];
rz(1.3560881) q[2];
sx q[2];
rz(-2.4248059) q[2];
sx q[2];
rz(1.4625664) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8902379) q[1];
sx q[1];
rz(-1.6511464) q[1];
sx q[1];
rz(-1.608196) q[1];
rz(-pi) q[2];
rz(2.1330256) q[3];
sx q[3];
rz(-1.3763714) q[3];
sx q[3];
rz(-2.090522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0884023) q[2];
sx q[2];
rz(-1.3206427) q[2];
sx q[2];
rz(-2.1957446) q[2];
rz(-2.451402) q[3];
sx q[3];
rz(-1.2226356) q[3];
sx q[3];
rz(-1.8405731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90203862) q[0];
sx q[0];
rz(-1.6033462) q[0];
sx q[0];
rz(-0.70815804) q[0];
rz(-3.0802849) q[1];
sx q[1];
rz(-0.61534449) q[1];
sx q[1];
rz(-1.4475488) q[1];
rz(1.4881143) q[2];
sx q[2];
rz(-1.2808242) q[2];
sx q[2];
rz(-2.6416953) q[2];
rz(1.7508909) q[3];
sx q[3];
rz(-1.4994499) q[3];
sx q[3];
rz(-2.3105619) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
