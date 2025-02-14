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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24841079) q[0];
sx q[0];
rz(-2.0201004) q[0];
sx q[0];
rz(0.23907678) q[0];
rz(-pi) q[1];
rz(-1.9196366) q[2];
sx q[2];
rz(-1.1695122) q[2];
sx q[2];
rz(2.2711585) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.5149821) q[1];
sx q[1];
rz(-1.1272073) q[1];
sx q[1];
rz(2.842998) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[1];
rz(2.8087372) q[2];
sx q[2];
rz(-0.70781195) q[2];
sx q[2];
rz(0.037192496) q[2];
rz(-2.228179) q[3];
sx q[3];
rz(-2.1593058) q[3];
sx q[3];
rz(1.6325525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34441242) q[0];
sx q[0];
rz(-2.0672815) q[0];
sx q[0];
rz(2.9194226) q[0];
rz(0.19368681) q[1];
sx q[1];
rz(-1.8733459) q[1];
sx q[1];
rz(-2.8817315) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4862999) q[0];
sx q[0];
rz(-1.3278786) q[0];
sx q[0];
rz(-1.1134558) q[0];
x q[1];
rz(-1.3919984) q[2];
sx q[2];
rz(-1.4646718) q[2];
sx q[2];
rz(0.44580844) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4187964) q[1];
sx q[1];
rz(-1.3592615) q[1];
sx q[1];
rz(1.3822894) q[1];
rz(-0.75921311) q[3];
sx q[3];
rz(-1.9369643) q[3];
sx q[3];
rz(0.002634332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.76693177) q[2];
sx q[2];
rz(-2.5587475) q[2];
sx q[2];
rz(-2.984821) q[2];
rz(-2.3356656) q[3];
sx q[3];
rz(-1.0904049) q[3];
sx q[3];
rz(0.52483112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2369775) q[0];
sx q[0];
rz(-3.0293063) q[0];
sx q[0];
rz(1.9554546) q[0];
rz(-3.0778432) q[1];
sx q[1];
rz(-1.6155764) q[1];
sx q[1];
rz(0.41248163) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29152399) q[0];
sx q[0];
rz(-2.0578034) q[0];
sx q[0];
rz(2.1662232) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3774728) q[2];
sx q[2];
rz(-2.4524322) q[2];
sx q[2];
rz(-1.8887519) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8349065) q[1];
sx q[1];
rz(-1.2962771) q[1];
sx q[1];
rz(1.7815018) q[1];
rz(-pi) q[2];
rz(-0.8728506) q[3];
sx q[3];
rz(-0.79668364) q[3];
sx q[3];
rz(2.7484675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5731262) q[2];
sx q[2];
rz(-2.7012479) q[2];
sx q[2];
rz(1.9754515) q[2];
rz(-2.3447573) q[3];
sx q[3];
rz(-1.3692057) q[3];
sx q[3];
rz(2.3603175) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9154938) q[0];
sx q[0];
rz(-1.6330566) q[0];
sx q[0];
rz(-1.728212) q[0];
rz(0.49890292) q[1];
sx q[1];
rz(-1.7585124) q[1];
sx q[1];
rz(-1.8477207) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.897936) q[0];
sx q[0];
rz(-1.597628) q[0];
sx q[0];
rz(1.6559589) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2433238) q[2];
sx q[2];
rz(-1.4988384) q[2];
sx q[2];
rz(-1.5775934) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1509959) q[1];
sx q[1];
rz(-1.5739723) q[1];
sx q[1];
rz(-0.062009723) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0741229) q[3];
sx q[3];
rz(-1.6327792) q[3];
sx q[3];
rz(2.8648479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.52900195) q[2];
sx q[2];
rz(-2.4606073) q[2];
sx q[2];
rz(0.16981086) q[2];
rz(-0.42810193) q[3];
sx q[3];
rz(-1.4295108) q[3];
sx q[3];
rz(0.16038173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.661628) q[0];
sx q[0];
rz(-2.7492838) q[0];
sx q[0];
rz(-0.83537927) q[0];
rz(2.5021878) q[1];
sx q[1];
rz(-2.4765922) q[1];
sx q[1];
rz(-1.0909874) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22459659) q[0];
sx q[0];
rz(-2.4451655) q[0];
sx q[0];
rz(0.019231872) q[0];
x q[1];
rz(1.4032989) q[2];
sx q[2];
rz(-1.27503) q[2];
sx q[2];
rz(2.5633321) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1941293) q[1];
sx q[1];
rz(-2.3724679) q[1];
sx q[1];
rz(-2.1544653) q[1];
rz(-pi) q[2];
rz(1.0476607) q[3];
sx q[3];
rz(-1.8561279) q[3];
sx q[3];
rz(0.77013515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.368448) q[2];
sx q[2];
rz(-1.1870563) q[2];
sx q[2];
rz(-3.1375569) q[2];
rz(-2.0648547) q[3];
sx q[3];
rz(-0.69759798) q[3];
sx q[3];
rz(-2.9546886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3183597) q[0];
sx q[0];
rz(-1.7802745) q[0];
sx q[0];
rz(-2.5079492) q[0];
rz(2.7404495) q[1];
sx q[1];
rz(-0.70924962) q[1];
sx q[1];
rz(0.69222442) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.322583) q[0];
sx q[0];
rz(-2.2624709) q[0];
sx q[0];
rz(-2.5965967) q[0];
rz(-1.2622617) q[2];
sx q[2];
rz(-2.3156328) q[2];
sx q[2];
rz(0.334563) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8521523) q[1];
sx q[1];
rz(-1.0045144) q[1];
sx q[1];
rz(0.29636611) q[1];
rz(-2.9561437) q[3];
sx q[3];
rz(-1.4223961) q[3];
sx q[3];
rz(0.68391358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.20574084) q[2];
sx q[2];
rz(-2.3571099) q[2];
sx q[2];
rz(1.240823) q[2];
rz(-1.1848909) q[3];
sx q[3];
rz(-1.3013867) q[3];
sx q[3];
rz(-1.544781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1320837) q[0];
sx q[0];
rz(-1.8385831) q[0];
sx q[0];
rz(1.7708154) q[0];
rz(2.7235183) q[1];
sx q[1];
rz(-1.4825753) q[1];
sx q[1];
rz(0.48666993) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91557568) q[0];
sx q[0];
rz(-0.56779438) q[0];
sx q[0];
rz(1.3894807) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1407632) q[2];
sx q[2];
rz(-0.62997937) q[2];
sx q[2];
rz(-0.85106817) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8312953) q[1];
sx q[1];
rz(-0.26997176) q[1];
sx q[1];
rz(-0.84862535) q[1];
rz(-pi) q[2];
rz(0.47007665) q[3];
sx q[3];
rz(-1.6127018) q[3];
sx q[3];
rz(-2.3500729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.70249867) q[2];
sx q[2];
rz(-2.0926026) q[2];
sx q[2];
rz(-2.9443963) q[2];
rz(-0.50944734) q[3];
sx q[3];
rz(-0.26064894) q[3];
sx q[3];
rz(-0.82836866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8615123) q[0];
sx q[0];
rz(-2.0674288) q[0];
sx q[0];
rz(1.3664838) q[0];
rz(-0.81661433) q[1];
sx q[1];
rz(-2.0520703) q[1];
sx q[1];
rz(2.8588967) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3168516) q[0];
sx q[0];
rz(-2.5159989) q[0];
sx q[0];
rz(1.8815142) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31158547) q[2];
sx q[2];
rz(-1.6610314) q[2];
sx q[2];
rz(2.2180706) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.41425175) q[1];
sx q[1];
rz(-1.7569572) q[1];
sx q[1];
rz(0.8826137) q[1];
rz(2.5211996) q[3];
sx q[3];
rz(-0.61479688) q[3];
sx q[3];
rz(2.7063587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1064328) q[2];
sx q[2];
rz(-1.0059493) q[2];
sx q[2];
rz(-0.98313588) q[2];
rz(1.8631009) q[3];
sx q[3];
rz(-2.216414) q[3];
sx q[3];
rz(2.2127051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66016692) q[0];
sx q[0];
rz(-1.7919414) q[0];
sx q[0];
rz(0.31326681) q[0];
rz(-0.52472862) q[1];
sx q[1];
rz(-0.26291651) q[1];
sx q[1];
rz(1.252334) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1850644) q[0];
sx q[0];
rz(-1.6913497) q[0];
sx q[0];
rz(-0.53114364) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76340126) q[2];
sx q[2];
rz(-2.0138676) q[2];
sx q[2];
rz(0.099066548) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5460427) q[1];
sx q[1];
rz(-2.3417406) q[1];
sx q[1];
rz(-0.78592664) q[1];
x q[2];
rz(0.081324025) q[3];
sx q[3];
rz(-2.3235851) q[3];
sx q[3];
rz(0.52018702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9423882) q[2];
sx q[2];
rz(-0.62493268) q[2];
sx q[2];
rz(1.3814629) q[2];
rz(-0.30113164) q[3];
sx q[3];
rz(-2.1589409) q[3];
sx q[3];
rz(-0.38100955) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23406601) q[0];
sx q[0];
rz(-1.0028239) q[0];
sx q[0];
rz(1.3475077) q[0];
rz(2.2566336) q[1];
sx q[1];
rz(-0.62565175) q[1];
sx q[1];
rz(1.7431097) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8475474) q[0];
sx q[0];
rz(-1.1963524) q[0];
sx q[0];
rz(-3.1301296) q[0];
rz(-pi) q[1];
x q[1];
rz(0.18357205) q[2];
sx q[2];
rz(-2.2677448) q[2];
sx q[2];
rz(1.1810034) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2513548) q[1];
sx q[1];
rz(-1.6511464) q[1];
sx q[1];
rz(-1.608196) q[1];
x q[2];
rz(2.9129254) q[3];
sx q[3];
rz(-2.1211984) q[3];
sx q[3];
rz(-2.500734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.053190319) q[2];
sx q[2];
rz(-1.3206427) q[2];
sx q[2];
rz(0.94584805) q[2];
rz(2.451402) q[3];
sx q[3];
rz(-1.2226356) q[3];
sx q[3];
rz(-1.3010196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90203862) q[0];
sx q[0];
rz(-1.5382465) q[0];
sx q[0];
rz(2.4334346) q[0];
rz(-3.0802849) q[1];
sx q[1];
rz(-0.61534449) q[1];
sx q[1];
rz(-1.4475488) q[1];
rz(2.8715677) q[2];
sx q[2];
rz(-2.8403828) q[2];
sx q[2];
rz(-2.3595911) q[2];
rz(0.072515247) q[3];
sx q[3];
rz(-1.7504277) q[3];
sx q[3];
rz(-0.75274368) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
