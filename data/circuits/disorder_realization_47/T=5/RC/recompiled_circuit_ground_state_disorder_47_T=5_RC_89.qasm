OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1563675) q[0];
sx q[0];
rz(-1.2824143) q[0];
sx q[0];
rz(3.0106944) q[0];
rz(-2.6137597) q[1];
sx q[1];
rz(-0.37017828) q[1];
sx q[1];
rz(0.17732492) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41577121) q[0];
sx q[0];
rz(-2.3292589) q[0];
sx q[0];
rz(-1.2052631) q[0];
rz(-pi) q[1];
rz(-1.2700291) q[2];
sx q[2];
rz(-1.352515) q[2];
sx q[2];
rz(1.0839562) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0541815) q[1];
sx q[1];
rz(-2.0435395) q[1];
sx q[1];
rz(-2.8768538) q[1];
rz(-1.9779284) q[3];
sx q[3];
rz(-0.61591776) q[3];
sx q[3];
rz(0.27995279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6191972) q[2];
sx q[2];
rz(-2.0530687) q[2];
sx q[2];
rz(1.8001451) q[2];
rz(-1.7893192) q[3];
sx q[3];
rz(-1.6100581) q[3];
sx q[3];
rz(-1.2927829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.751048) q[0];
sx q[0];
rz(-1.9123257) q[0];
sx q[0];
rz(-2.6265662) q[0];
rz(0.36188778) q[1];
sx q[1];
rz(-1.3970951) q[1];
sx q[1];
rz(-0.99641689) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5549042) q[0];
sx q[0];
rz(-1.3350272) q[0];
sx q[0];
rz(0.20637189) q[0];
rz(-0.17536945) q[2];
sx q[2];
rz(-1.9992454) q[2];
sx q[2];
rz(-2.1074949) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.60153466) q[1];
sx q[1];
rz(-0.52943474) q[1];
sx q[1];
rz(1.387818) q[1];
x q[2];
rz(2.6985136) q[3];
sx q[3];
rz(-2.0650117) q[3];
sx q[3];
rz(0.48393238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4126052) q[2];
sx q[2];
rz(-0.94634405) q[2];
sx q[2];
rz(-1.6607355) q[2];
rz(2.6442773) q[3];
sx q[3];
rz(-0.9484843) q[3];
sx q[3];
rz(-2.4174387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6052674) q[0];
sx q[0];
rz(-1.241642) q[0];
sx q[0];
rz(1.9507677) q[0];
rz(-0.39769998) q[1];
sx q[1];
rz(-1.5813446) q[1];
sx q[1];
rz(-0.89888987) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9346823) q[0];
sx q[0];
rz(-2.795462) q[0];
sx q[0];
rz(-1.3735361) q[0];
rz(-pi) q[1];
rz(-1.6562881) q[2];
sx q[2];
rz(-1.7608661) q[2];
sx q[2];
rz(-0.23641931) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7601493) q[1];
sx q[1];
rz(-0.53434125) q[1];
sx q[1];
rz(1.4830073) q[1];
x q[2];
rz(1.2736788) q[3];
sx q[3];
rz(-2.6032902) q[3];
sx q[3];
rz(0.75002128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1206104) q[2];
sx q[2];
rz(-1.1889428) q[2];
sx q[2];
rz(-2.9494542) q[2];
rz(-0.7297248) q[3];
sx q[3];
rz(-0.14484043) q[3];
sx q[3];
rz(0.63989583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14908734) q[0];
sx q[0];
rz(-2.3893116) q[0];
sx q[0];
rz(-1.8335023) q[0];
rz(-0.84450841) q[1];
sx q[1];
rz(-1.6788071) q[1];
sx q[1];
rz(0.62215296) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7956411) q[0];
sx q[0];
rz(-1.1591732) q[0];
sx q[0];
rz(0.53070416) q[0];
rz(-1.5160558) q[2];
sx q[2];
rz(-1.3856882) q[2];
sx q[2];
rz(2.6814044) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0627985) q[1];
sx q[1];
rz(-1.5192966) q[1];
sx q[1];
rz(-2.9777793) q[1];
x q[2];
rz(-0.31720576) q[3];
sx q[3];
rz(-1.460808) q[3];
sx q[3];
rz(-3.0134137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1239803) q[2];
sx q[2];
rz(-1.7698741) q[2];
sx q[2];
rz(-0.89517108) q[2];
rz(0.34919843) q[3];
sx q[3];
rz(-2.5694191) q[3];
sx q[3];
rz(-0.23381843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2772086) q[0];
sx q[0];
rz(-1.9338436) q[0];
sx q[0];
rz(2.5982017) q[0];
rz(-3.0282989) q[1];
sx q[1];
rz(-1.7430867) q[1];
sx q[1];
rz(-1.8494122) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1975511) q[0];
sx q[0];
rz(-0.8055976) q[0];
sx q[0];
rz(2.1155691) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2207556) q[2];
sx q[2];
rz(-1.2652768) q[2];
sx q[2];
rz(-0.68543514) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0185011) q[1];
sx q[1];
rz(-0.95884174) q[1];
sx q[1];
rz(-0.83142821) q[1];
rz(-pi) q[2];
rz(0.73907799) q[3];
sx q[3];
rz(-1.0024286) q[3];
sx q[3];
rz(-3.0773602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8285404) q[2];
sx q[2];
rz(-1.7834168) q[2];
sx q[2];
rz(-2.6521315) q[2];
rz(0.66568565) q[3];
sx q[3];
rz(-2.5299215) q[3];
sx q[3];
rz(2.3556975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2419484) q[0];
sx q[0];
rz(-0.90270942) q[0];
sx q[0];
rz(-0.06074252) q[0];
rz(-1.2784866) q[1];
sx q[1];
rz(-2.2272031) q[1];
sx q[1];
rz(-1.0260822) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7728108) q[0];
sx q[0];
rz(-2.368481) q[0];
sx q[0];
rz(-0.85526534) q[0];
rz(0.32569484) q[2];
sx q[2];
rz(-2.1148588) q[2];
sx q[2];
rz(-2.0607255) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7674305) q[1];
sx q[1];
rz(-2.8547127) q[1];
sx q[1];
rz(0.5021023) q[1];
rz(-pi) q[2];
rz(1.2450663) q[3];
sx q[3];
rz(-1.0661849) q[3];
sx q[3];
rz(0.65629279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.25980276) q[2];
sx q[2];
rz(-1.5564352) q[2];
sx q[2];
rz(-3.0628824) q[2];
rz(1.4824661) q[3];
sx q[3];
rz(-0.31646287) q[3];
sx q[3];
rz(0.44221529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.215613) q[0];
sx q[0];
rz(-2.7613566) q[0];
sx q[0];
rz(-1.6967787) q[0];
rz(-2.3346057) q[1];
sx q[1];
rz(-1.746256) q[1];
sx q[1];
rz(-0.011946202) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0204457) q[0];
sx q[0];
rz(-0.63948005) q[0];
sx q[0];
rz(1.5259652) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19334601) q[2];
sx q[2];
rz(-2.6259661) q[2];
sx q[2];
rz(0.28179291) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.641549) q[1];
sx q[1];
rz(-2.1807359) q[1];
sx q[1];
rz(-0.088492101) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9254382) q[3];
sx q[3];
rz(-0.4222479) q[3];
sx q[3];
rz(-0.99340445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5936467) q[2];
sx q[2];
rz(-2.839489) q[2];
sx q[2];
rz(-0.83615237) q[2];
rz(-0.052308403) q[3];
sx q[3];
rz(-1.6736504) q[3];
sx q[3];
rz(0.93594319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1702105) q[0];
sx q[0];
rz(-2.2844071) q[0];
sx q[0];
rz(-0.49825391) q[0];
rz(-1.0055297) q[1];
sx q[1];
rz(-1.4509095) q[1];
sx q[1];
rz(3.0301869) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.446602) q[0];
sx q[0];
rz(-1.1518307) q[0];
sx q[0];
rz(-0.10563157) q[0];
rz(-pi) q[1];
rz(-2.7132052) q[2];
sx q[2];
rz(-1.7991788) q[2];
sx q[2];
rz(-1.0699748) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1131683) q[1];
sx q[1];
rz(-2.3481124) q[1];
sx q[1];
rz(1.9774578) q[1];
rz(-2.8008055) q[3];
sx q[3];
rz(-2.7772745) q[3];
sx q[3];
rz(-2.7076569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.19693836) q[2];
sx q[2];
rz(-1.3613746) q[2];
sx q[2];
rz(1.4062175) q[2];
rz(1.1574636) q[3];
sx q[3];
rz(-2.1486053) q[3];
sx q[3];
rz(-1.5895313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0208397) q[0];
sx q[0];
rz(-0.73796219) q[0];
sx q[0];
rz(-1.4388168) q[0];
rz(0.84398794) q[1];
sx q[1];
rz(-1.8344717) q[1];
sx q[1];
rz(0.2058952) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1263329) q[0];
sx q[0];
rz(-1.6113937) q[0];
sx q[0];
rz(0.6313398) q[0];
x q[1];
rz(1.5175784) q[2];
sx q[2];
rz(-0.88376617) q[2];
sx q[2];
rz(-1.991697) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3033324) q[1];
sx q[1];
rz(-0.15619677) q[1];
sx q[1];
rz(2.6172943) q[1];
rz(2.6456553) q[3];
sx q[3];
rz(-1.933799) q[3];
sx q[3];
rz(-2.7285226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0042808) q[2];
sx q[2];
rz(-1.0575123) q[2];
sx q[2];
rz(1.6005969) q[2];
rz(2.6565523) q[3];
sx q[3];
rz(-1.78616) q[3];
sx q[3];
rz(0.11390991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.4892905) q[0];
sx q[0];
rz(-0.052209608) q[0];
sx q[0];
rz(1.2963649) q[0];
rz(1.0507874) q[1];
sx q[1];
rz(-1.6810345) q[1];
sx q[1];
rz(-2.5349862) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8826128) q[0];
sx q[0];
rz(-1.4073696) q[0];
sx q[0];
rz(0.033074065) q[0];
x q[1];
rz(1.6982618) q[2];
sx q[2];
rz(-1.0147139) q[2];
sx q[2];
rz(0.13917222) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.21571017) q[1];
sx q[1];
rz(-1.6318738) q[1];
sx q[1];
rz(-0.30558773) q[1];
rz(-1.4539155) q[3];
sx q[3];
rz(-1.28834) q[3];
sx q[3];
rz(-2.4931537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7217094) q[2];
sx q[2];
rz(-2.0100644) q[2];
sx q[2];
rz(-0.69880542) q[2];
rz(1.4612259) q[3];
sx q[3];
rz(-2.1278087) q[3];
sx q[3];
rz(0.47718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2496495) q[0];
sx q[0];
rz(-1.6310545) q[0];
sx q[0];
rz(-1.7153836) q[0];
rz(-0.37829788) q[1];
sx q[1];
rz(-0.62722798) q[1];
sx q[1];
rz(2.41082) q[1];
rz(2.8332491) q[2];
sx q[2];
rz(-2.0258122) q[2];
sx q[2];
rz(-3.0292055) q[2];
rz(-3.0918269) q[3];
sx q[3];
rz(-0.55703041) q[3];
sx q[3];
rz(1.218474) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
