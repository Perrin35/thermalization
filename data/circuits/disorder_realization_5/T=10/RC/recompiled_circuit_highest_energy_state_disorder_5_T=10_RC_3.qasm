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
rz(2.2628885) q[0];
sx q[0];
rz(9.2287697) q[0];
rz(1.6485543) q[1];
sx q[1];
rz(-2.2033446) q[1];
sx q[1];
rz(-0.87747639) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75979027) q[0];
sx q[0];
rz(-0.50509278) q[0];
sx q[0];
rz(2.0272966) q[0];
x q[1];
rz(-2.4634667) q[2];
sx q[2];
rz(-0.52541753) q[2];
sx q[2];
rz(-0.12066387) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6266105) q[1];
sx q[1];
rz(-2.0143854) q[1];
sx q[1];
rz(-0.29859467) q[1];
rz(-2.9035197) q[3];
sx q[3];
rz(-2.0364982) q[3];
sx q[3];
rz(-2.4089486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8087372) q[2];
sx q[2];
rz(-2.4337807) q[2];
sx q[2];
rz(0.037192496) q[2];
rz(-2.228179) q[3];
sx q[3];
rz(-2.1593058) q[3];
sx q[3];
rz(-1.5090401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7971802) q[0];
sx q[0];
rz(-2.0672815) q[0];
sx q[0];
rz(-2.9194226) q[0];
rz(0.19368681) q[1];
sx q[1];
rz(-1.8733459) q[1];
sx q[1];
rz(-2.8817315) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65529275) q[0];
sx q[0];
rz(-1.813714) q[0];
sx q[0];
rz(2.0281369) q[0];
x q[1];
rz(2.1104576) q[2];
sx q[2];
rz(-0.20763131) q[2];
sx q[2];
rz(-0.59484824) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7227962) q[1];
sx q[1];
rz(-1.3592615) q[1];
sx q[1];
rz(-1.7593033) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50825676) q[3];
sx q[3];
rz(-0.82672366) q[3];
sx q[3];
rz(-1.9342157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.76693177) q[2];
sx q[2];
rz(-0.58284512) q[2];
sx q[2];
rz(-0.15677162) q[2];
rz(2.3356656) q[3];
sx q[3];
rz(-1.0904049) q[3];
sx q[3];
rz(2.6167615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2369775) q[0];
sx q[0];
rz(-0.1122864) q[0];
sx q[0];
rz(-1.9554546) q[0];
rz(0.063749464) q[1];
sx q[1];
rz(-1.6155764) q[1];
sx q[1];
rz(-2.729111) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97226364) q[0];
sx q[0];
rz(-2.0893851) q[0];
sx q[0];
rz(-2.5725468) q[0];
rz(0.89084741) q[2];
sx q[2];
rz(-1.6932704) q[2];
sx q[2];
rz(-2.9735931) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3220249) q[1];
sx q[1];
rz(-1.7734999) q[1];
sx q[1];
rz(0.2804108) q[1];
x q[2];
rz(-2.2354911) q[3];
sx q[3];
rz(-1.0933439) q[3];
sx q[3];
rz(2.4942644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5731262) q[2];
sx q[2];
rz(-2.7012479) q[2];
sx q[2];
rz(1.1661412) q[2];
rz(-0.79683534) q[3];
sx q[3];
rz(-1.7723869) q[3];
sx q[3];
rz(-0.78127512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
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
rz(1.2260988) q[0];
sx q[0];
rz(-1.6330566) q[0];
sx q[0];
rz(1.728212) q[0];
rz(-0.49890292) q[1];
sx q[1];
rz(-1.3830802) q[1];
sx q[1];
rz(1.2938719) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.897936) q[0];
sx q[0];
rz(-1.597628) q[0];
sx q[0];
rz(1.4856337) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8982688) q[2];
sx q[2];
rz(-1.6427543) q[2];
sx q[2];
rz(1.5775934) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7728987) q[1];
sx q[1];
rz(-3.0795018) q[1];
sx q[1];
rz(0.05120488) q[1];
rz(-2.3975375) q[3];
sx q[3];
rz(-3.0500055) q[3];
sx q[3];
rz(-2.0360144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6125907) q[2];
sx q[2];
rz(-2.4606073) q[2];
sx q[2];
rz(-0.16981086) q[2];
rz(0.42810193) q[3];
sx q[3];
rz(-1.7120818) q[3];
sx q[3];
rz(-2.9812109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47996461) q[0];
sx q[0];
rz(-0.39230883) q[0];
sx q[0];
rz(0.83537927) q[0];
rz(-0.63940489) q[1];
sx q[1];
rz(-2.4765922) q[1];
sx q[1];
rz(2.0506052) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3314455) q[0];
sx q[0];
rz(-1.5831328) q[0];
sx q[0];
rz(-0.69633616) q[0];
rz(0.29971896) q[2];
sx q[2];
rz(-1.4106361) q[2];
sx q[2];
rz(-2.0998139) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.450836) q[1];
sx q[1];
rz(-0.95162205) q[1];
sx q[1];
rz(-2.6515533) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0398846) q[3];
sx q[3];
rz(-0.58947488) q[3];
sx q[3];
rz(2.7949751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7731446) q[2];
sx q[2];
rz(-1.9545363) q[2];
sx q[2];
rz(-3.1375569) q[2];
rz(-2.0648547) q[3];
sx q[3];
rz(-0.69759798) q[3];
sx q[3];
rz(0.18690404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-0.70924962) q[1];
sx q[1];
rz(-2.4493682) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56256912) q[0];
sx q[0];
rz(-2.2897567) q[0];
sx q[0];
rz(-1.0114874) q[0];
rz(-2.8234286) q[2];
sx q[2];
rz(-2.3468694) q[2];
sx q[2];
rz(-0.77407167) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0226188) q[1];
sx q[1];
rz(-1.8198208) q[1];
sx q[1];
rz(0.98414661) q[1];
rz(-pi) q[2];
rz(2.4603087) q[3];
sx q[3];
rz(-2.9046106) q[3];
sx q[3];
rz(2.9221688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.20574084) q[2];
sx q[2];
rz(-0.78448272) q[2];
sx q[2];
rz(-1.9007696) q[2];
rz(-1.1848909) q[3];
sx q[3];
rz(-1.840206) q[3];
sx q[3];
rz(-1.5968116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1320837) q[0];
sx q[0];
rz(-1.8385831) q[0];
sx q[0];
rz(-1.7708154) q[0];
rz(0.41807434) q[1];
sx q[1];
rz(-1.6590174) q[1];
sx q[1];
rz(-2.6549227) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.226017) q[0];
sx q[0];
rz(-2.5737983) q[0];
sx q[0];
rz(-1.752112) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7667747) q[2];
sx q[2];
rz(-2.0897802) q[2];
sx q[2];
rz(-1.6199552) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8312953) q[1];
sx q[1];
rz(-0.26997176) q[1];
sx q[1];
rz(0.84862535) q[1];
rz(3.0492856) q[3];
sx q[3];
rz(-2.669791) q[3];
sx q[3];
rz(-2.2800454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.70249867) q[2];
sx q[2];
rz(-2.0926026) q[2];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8615123) q[0];
sx q[0];
rz(-2.0674288) q[0];
sx q[0];
rz(1.7751088) q[0];
rz(-2.3249783) q[1];
sx q[1];
rz(-1.0895224) q[1];
sx q[1];
rz(-0.28269592) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3168516) q[0];
sx q[0];
rz(-0.62559375) q[0];
sx q[0];
rz(-1.2600785) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4760232) q[2];
sx q[2];
rz(-1.2605209) q[2];
sx q[2];
rz(-0.6182593) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0055089) q[1];
sx q[1];
rz(-2.2448531) q[1];
sx q[1];
rz(2.9024209) q[1];
rz(-pi) q[2];
rz(-1.9603086) q[3];
sx q[3];
rz(-1.0822902) q[3];
sx q[3];
rz(2.8581885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.035159811) q[2];
sx q[2];
rz(-1.0059493) q[2];
sx q[2];
rz(-0.98313588) q[2];
rz(1.2784917) q[3];
sx q[3];
rz(-0.92517868) q[3];
sx q[3];
rz(-0.92888752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4814257) q[0];
sx q[0];
rz(-1.7919414) q[0];
sx q[0];
rz(2.8283258) q[0];
rz(2.616864) q[1];
sx q[1];
rz(-2.8786761) q[1];
sx q[1];
rz(1.8892586) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9577873) q[0];
sx q[0];
rz(-2.5982214) q[0];
sx q[0];
rz(-0.23475351) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60151327) q[2];
sx q[2];
rz(-0.85962112) q[2];
sx q[2];
rz(-1.892923) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.63284167) q[1];
sx q[1];
rz(-2.1023884) q[1];
sx q[1];
rz(0.94137194) q[1];
rz(-pi) q[2];
rz(2.325237) q[3];
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
rz(-2.840461) q[3];
sx q[3];
rz(-0.98265177) q[3];
sx q[3];
rz(2.7605831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9075266) q[0];
sx q[0];
rz(-1.0028239) q[0];
sx q[0];
rz(-1.3475077) q[0];
rz(-2.2566336) q[1];
sx q[1];
rz(-0.62565175) q[1];
sx q[1];
rz(-1.7431097) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32537726) q[0];
sx q[0];
rz(-2.7669816) q[0];
sx q[0];
rz(1.5999567) q[0];
x q[1];
rz(1.3560881) q[2];
sx q[2];
rz(-0.71678679) q[2];
sx q[2];
rz(-1.4625664) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.32244476) q[1];
sx q[1];
rz(-1.5335173) q[1];
sx q[1];
rz(-0.080406043) q[1];
rz(-2.9129254) q[3];
sx q[3];
rz(-2.1211984) q[3];
sx q[3];
rz(-0.64085863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0884023) q[2];
sx q[2];
rz(-1.3206427) q[2];
sx q[2];
rz(-0.94584805) q[2];
rz(0.69019067) q[3];
sx q[3];
rz(-1.2226356) q[3];
sx q[3];
rz(1.3010196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.90203862) q[0];
sx q[0];
rz(-1.6033462) q[0];
sx q[0];
rz(-0.70815804) q[0];
rz(0.061307727) q[1];
sx q[1];
rz(-0.61534449) q[1];
sx q[1];
rz(-1.4475488) q[1];
rz(1.6534783) q[2];
sx q[2];
rz(-1.8607685) q[2];
sx q[2];
rz(0.49989732) q[2];
rz(-1.9504299) q[3];
sx q[3];
rz(-2.9480231) q[3];
sx q[3];
rz(2.775016) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
