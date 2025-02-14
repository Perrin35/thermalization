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
rz(1.5647178) q[0];
sx q[0];
rz(-1.325542) q[0];
sx q[0];
rz(2.5817885) q[0];
rz(-1.7436104) q[1];
sx q[1];
rz(-1.3275194) q[1];
sx q[1];
rz(-1.3865857) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5450534) q[0];
sx q[0];
rz(-0.60205108) q[0];
sx q[0];
rz(-2.028378) q[0];
rz(2.8459889) q[2];
sx q[2];
rz(-1.6787375) q[2];
sx q[2];
rz(0.88376494) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.4999285) q[1];
sx q[1];
rz(-1.4870502) q[1];
sx q[1];
rz(-0.20831627) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5941526) q[3];
sx q[3];
rz(-2.8508715) q[3];
sx q[3];
rz(-0.61222044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6036512) q[2];
sx q[2];
rz(-2.5141022) q[2];
sx q[2];
rz(0.44169912) q[2];
rz(2.3407827) q[3];
sx q[3];
rz(-1.8947314) q[3];
sx q[3];
rz(2.7504564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2342728) q[0];
sx q[0];
rz(-1.3751625) q[0];
sx q[0];
rz(-0.32670879) q[0];
rz(-1.7960583) q[1];
sx q[1];
rz(-1.6163454) q[1];
sx q[1];
rz(-0.84652841) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1382111) q[0];
sx q[0];
rz(-1.4633006) q[0];
sx q[0];
rz(1.18446) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88794215) q[2];
sx q[2];
rz(-1.789621) q[2];
sx q[2];
rz(1.3252979) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.1607779) q[1];
sx q[1];
rz(-2.6123568) q[1];
sx q[1];
rz(1.5917626) q[1];
x q[2];
rz(-2.8277581) q[3];
sx q[3];
rz(-2.3646128) q[3];
sx q[3];
rz(-1.8557567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9091866) q[2];
sx q[2];
rz(-2.0441983) q[2];
sx q[2];
rz(-0.48420134) q[2];
rz(1.3062612) q[3];
sx q[3];
rz(-0.56070119) q[3];
sx q[3];
rz(-1.9561214) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66121286) q[0];
sx q[0];
rz(-0.73037195) q[0];
sx q[0];
rz(2.6015306) q[0];
rz(-1.9909319) q[1];
sx q[1];
rz(-1.2723424) q[1];
sx q[1];
rz(2.5475492) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6912508) q[0];
sx q[0];
rz(-1.461554) q[0];
sx q[0];
rz(1.7418377) q[0];
rz(-0.39806367) q[2];
sx q[2];
rz(-0.69483611) q[2];
sx q[2];
rz(0.24206012) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.6642528) q[1];
sx q[1];
rz(-0.76271933) q[1];
sx q[1];
rz(-1.4422756) q[1];
x q[2];
rz(1.698231) q[3];
sx q[3];
rz(-2.6563247) q[3];
sx q[3];
rz(1.8125774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8192886) q[2];
sx q[2];
rz(-2.1257336) q[2];
sx q[2];
rz(0.6907531) q[2];
rz(0.47762075) q[3];
sx q[3];
rz(-2.728929) q[3];
sx q[3];
rz(-2.7121108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-1.8884884) q[0];
sx q[0];
rz(-0.34585837) q[0];
sx q[0];
rz(-0.75230569) q[0];
rz(2.2350156) q[1];
sx q[1];
rz(-0.38175672) q[1];
sx q[1];
rz(-2.9125772) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2370214) q[0];
sx q[0];
rz(-1.8638049) q[0];
sx q[0];
rz(1.8051534) q[0];
rz(-pi) q[1];
rz(-0.53900881) q[2];
sx q[2];
rz(-1.3752364) q[2];
sx q[2];
rz(0.53949088) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6743857) q[1];
sx q[1];
rz(-1.3027281) q[1];
sx q[1];
rz(1.3246791) q[1];
rz(-pi) q[2];
rz(1.9360156) q[3];
sx q[3];
rz(-2.1250279) q[3];
sx q[3];
rz(-0.54868719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0166246) q[2];
sx q[2];
rz(-0.11398537) q[2];
sx q[2];
rz(-1.0329186) q[2];
rz(2.0756663) q[3];
sx q[3];
rz(-1.67098) q[3];
sx q[3];
rz(1.3037995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8835939) q[0];
sx q[0];
rz(-3.0553387) q[0];
sx q[0];
rz(-1.8121383) q[0];
rz(-2.6138002) q[1];
sx q[1];
rz(-1.8254435) q[1];
sx q[1];
rz(-0.39906183) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85241991) q[0];
sx q[0];
rz(-1.6483083) q[0];
sx q[0];
rz(1.2184675) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1588604) q[2];
sx q[2];
rz(-2.2592681) q[2];
sx q[2];
rz(0.43718034) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.1414772) q[1];
sx q[1];
rz(-2.2908604) q[1];
sx q[1];
rz(-1.9898371) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9513017) q[3];
sx q[3];
rz(-2.2566593) q[3];
sx q[3];
rz(-3.0113285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1479147) q[2];
sx q[2];
rz(-0.75239158) q[2];
sx q[2];
rz(-1.0047151) q[2];
rz(-2.8435454) q[3];
sx q[3];
rz(-1.7883251) q[3];
sx q[3];
rz(-2.0731488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7816958) q[0];
sx q[0];
rz(-1.6258465) q[0];
sx q[0];
rz(-1.9500649) q[0];
rz(-2.6021992) q[1];
sx q[1];
rz(-1.0788147) q[1];
sx q[1];
rz(2.3753812) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3193885) q[0];
sx q[0];
rz(-0.66492101) q[0];
sx q[0];
rz(0.98916905) q[0];
x q[1];
rz(-0.64935342) q[2];
sx q[2];
rz(-1.0217862) q[2];
sx q[2];
rz(1.8446814) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.898452) q[1];
sx q[1];
rz(-2.2808373) q[1];
sx q[1];
rz(0.30177994) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2003464) q[3];
sx q[3];
rz(-1.660957) q[3];
sx q[3];
rz(2.3925635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3300276) q[2];
sx q[2];
rz(-1.4469688) q[2];
sx q[2];
rz(-2.2583029) q[2];
rz(-2.6106994) q[3];
sx q[3];
rz(-0.58404946) q[3];
sx q[3];
rz(-0.77597165) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74060488) q[0];
sx q[0];
rz(-2.8477836) q[0];
sx q[0];
rz(-2.4921932) q[0];
rz(0.3736639) q[1];
sx q[1];
rz(-2.3039736) q[1];
sx q[1];
rz(-1.1423473) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72066618) q[0];
sx q[0];
rz(-2.775179) q[0];
sx q[0];
rz(2.0271676) q[0];
x q[1];
rz(-2.9317103) q[2];
sx q[2];
rz(-1.58605) q[2];
sx q[2];
rz(-1.5535056) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5264633) q[1];
sx q[1];
rz(-2.3285236) q[1];
sx q[1];
rz(2.710756) q[1];
rz(-pi) q[2];
rz(-1.6616529) q[3];
sx q[3];
rz(-2.6808219) q[3];
sx q[3];
rz(2.6449869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.00804) q[2];
sx q[2];
rz(-0.69591659) q[2];
sx q[2];
rz(-0.35888654) q[2];
rz(2.533203) q[3];
sx q[3];
rz(-1.5865934) q[3];
sx q[3];
rz(-2.2801094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1818039) q[0];
sx q[0];
rz(-0.53878468) q[0];
sx q[0];
rz(2.1754919) q[0];
rz(-2.6995662) q[1];
sx q[1];
rz(-1.288488) q[1];
sx q[1];
rz(-0.40837049) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3672542) q[0];
sx q[0];
rz(-1.5440732) q[0];
sx q[0];
rz(1.7033808) q[0];
x q[1];
rz(-0.14773519) q[2];
sx q[2];
rz(-1.6134173) q[2];
sx q[2];
rz(3.0377667) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.30701559) q[1];
sx q[1];
rz(-0.6041827) q[1];
sx q[1];
rz(-1.8240364) q[1];
rz(-pi) q[2];
rz(0.28520198) q[3];
sx q[3];
rz(-0.75213471) q[3];
sx q[3];
rz(1.7908781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8926706) q[2];
sx q[2];
rz(-0.20888027) q[2];
sx q[2];
rz(-1.199031) q[2];
rz(1.2976973) q[3];
sx q[3];
rz(-1.0782995) q[3];
sx q[3];
rz(3.1395932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2099828) q[0];
sx q[0];
rz(-1.7117806) q[0];
sx q[0];
rz(-1.1071052) q[0];
rz(2.9529086) q[1];
sx q[1];
rz(-2.76177) q[1];
sx q[1];
rz(1.8170549) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0640776) q[0];
sx q[0];
rz(-2.0821838) q[0];
sx q[0];
rz(-1.8989424) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1679384) q[2];
sx q[2];
rz(-1.5533548) q[2];
sx q[2];
rz(-0.23250015) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.75809385) q[1];
sx q[1];
rz(-2.1111087) q[1];
sx q[1];
rz(3.1385149) q[1];
rz(0.56745133) q[3];
sx q[3];
rz(-2.3559915) q[3];
sx q[3];
rz(-0.42450617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.3482427) q[2];
sx q[2];
rz(-0.57955727) q[2];
sx q[2];
rz(-2.9938475) q[2];
rz(0.82792264) q[3];
sx q[3];
rz(-1.7235651) q[3];
sx q[3];
rz(-0.9399606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70656908) q[0];
sx q[0];
rz(-2.7913385) q[0];
sx q[0];
rz(-1.6178004) q[0];
rz(-1.891547) q[1];
sx q[1];
rz(-0.83386546) q[1];
sx q[1];
rz(-0.42125431) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6235879) q[0];
sx q[0];
rz(-0.54749085) q[0];
sx q[0];
rz(2.3836552) q[0];
x q[1];
rz(1.4827316) q[2];
sx q[2];
rz(-2.1906664) q[2];
sx q[2];
rz(-1.6234835) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7819216) q[1];
sx q[1];
rz(-0.087347833) q[1];
sx q[1];
rz(-3.0328301) q[1];
rz(-pi) q[2];
rz(-0.45205595) q[3];
sx q[3];
rz(-1.3516507) q[3];
sx q[3];
rz(-2.5170086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0201575) q[2];
sx q[2];
rz(-0.73548663) q[2];
sx q[2];
rz(0.45041931) q[2];
rz(-1.2823229) q[3];
sx q[3];
rz(-2.4560792) q[3];
sx q[3];
rz(0.57257563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12238518) q[0];
sx q[0];
rz(-1.5337802) q[0];
sx q[0];
rz(1.5564729) q[0];
rz(-2.4284594) q[1];
sx q[1];
rz(-1.5596371) q[1];
sx q[1];
rz(0.023991931) q[1];
rz(-2.6789011) q[2];
sx q[2];
rz(-2.0491055) q[2];
sx q[2];
rz(-2.9841205) q[2];
rz(0.080056277) q[3];
sx q[3];
rz(-1.4984984) q[3];
sx q[3];
rz(0.62582918) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
