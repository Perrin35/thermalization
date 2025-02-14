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
rz(-0.5598042) q[0];
rz(1.3979823) q[1];
sx q[1];
rz(-1.8140732) q[1];
sx q[1];
rz(-1.7550069) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35980269) q[0];
sx q[0];
rz(-1.823678) q[0];
sx q[0];
rz(-2.1232312) q[0];
rz(-pi) q[1];
rz(-0.29560372) q[2];
sx q[2];
rz(-1.6787375) q[2];
sx q[2];
rz(0.88376494) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0530449) q[1];
sx q[1];
rz(-1.7783718) q[1];
sx q[1];
rz(1.6563841) q[1];
rz(-1.5941526) q[3];
sx q[3];
rz(-0.29072116) q[3];
sx q[3];
rz(-0.61222044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5379415) q[2];
sx q[2];
rz(-2.5141022) q[2];
sx q[2];
rz(-0.44169912) q[2];
rz(-0.80080992) q[3];
sx q[3];
rz(-1.8947314) q[3];
sx q[3];
rz(2.7504564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2342728) q[0];
sx q[0];
rz(-1.3751625) q[0];
sx q[0];
rz(0.32670879) q[0];
rz(1.3455343) q[1];
sx q[1];
rz(-1.6163454) q[1];
sx q[1];
rz(-0.84652841) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5237997) q[0];
sx q[0];
rz(-1.1868068) q[0];
sx q[0];
rz(3.0256173) q[0];
rz(-pi) q[1];
rz(-0.2791762) q[2];
sx q[2];
rz(-2.2344032) q[2];
sx q[2];
rz(-0.070726591) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4281171) q[1];
sx q[1];
rz(-1.5813811) q[1];
sx q[1];
rz(1.0416563) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75196619) q[3];
sx q[3];
rz(-1.3526254) q[3];
sx q[3];
rz(3.0840517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9091866) q[2];
sx q[2];
rz(-1.0973944) q[2];
sx q[2];
rz(-2.6573913) q[2];
rz(-1.8353315) q[3];
sx q[3];
rz(-0.56070119) q[3];
sx q[3];
rz(1.1854712) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66121286) q[0];
sx q[0];
rz(-2.4112207) q[0];
sx q[0];
rz(2.6015306) q[0];
rz(1.1506608) q[1];
sx q[1];
rz(-1.2723424) q[1];
sx q[1];
rz(2.5475492) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4580989) q[0];
sx q[0];
rz(-2.9389295) q[0];
sx q[0];
rz(0.99838241) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65512983) q[2];
sx q[2];
rz(-1.821604) q[2];
sx q[2];
rz(2.1252968) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.81338367) q[1];
sx q[1];
rz(-1.4821308) q[1];
sx q[1];
rz(0.81221272) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0888799) q[3];
sx q[3];
rz(-1.6301117) q[3];
sx q[3];
rz(2.7869567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8192886) q[2];
sx q[2];
rz(-1.015859) q[2];
sx q[2];
rz(0.6907531) q[2];
rz(2.6639719) q[3];
sx q[3];
rz(-0.4126637) q[3];
sx q[3];
rz(0.42948183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.59531) q[0];
sx q[0];
rz(-2.7685099) q[0];
sx q[0];
rz(2.4855916) q[0];
rz(-pi) q[1];
rz(-0.53900881) q[2];
sx q[2];
rz(-1.3752364) q[2];
sx q[2];
rz(-2.6021018) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4333231) q[1];
sx q[1];
rz(-2.7796942) q[1];
sx q[1];
rz(-0.72558484) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5563779) q[3];
sx q[3];
rz(-1.8793595) q[3];
sx q[3];
rz(-1.2206961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1249681) q[2];
sx q[2];
rz(-0.11398537) q[2];
sx q[2];
rz(-1.0329186) q[2];
rz(-2.0756663) q[3];
sx q[3];
rz(-1.4706127) q[3];
sx q[3];
rz(-1.8377931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8835939) q[0];
sx q[0];
rz(-3.0553387) q[0];
sx q[0];
rz(1.8121383) q[0];
rz(2.6138002) q[1];
sx q[1];
rz(-1.3161491) q[1];
sx q[1];
rz(2.7425308) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.746839) q[0];
sx q[0];
rz(-1.9220214) q[0];
sx q[0];
rz(-0.082562692) q[0];
rz(2.40995) q[2];
sx q[2];
rz(-1.8851309) q[2];
sx q[2];
rz(-0.86282496) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.45259055) q[1];
sx q[1];
rz(-2.3277644) q[1];
sx q[1];
rz(-0.43431538) q[1];
rz(-pi) q[2];
rz(-0.87598975) q[3];
sx q[3];
rz(-1.717698) q[3];
sx q[3];
rz(1.5796652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.993678) q[2];
sx q[2];
rz(-2.3892011) q[2];
sx q[2];
rz(-1.0047151) q[2];
rz(-2.8435454) q[3];
sx q[3];
rz(-1.3532676) q[3];
sx q[3];
rz(-1.0684439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7816958) q[0];
sx q[0];
rz(-1.6258465) q[0];
sx q[0];
rz(-1.1915278) q[0];
rz(2.6021992) q[1];
sx q[1];
rz(-2.0627779) q[1];
sx q[1];
rz(-0.76621145) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82220413) q[0];
sx q[0];
rz(-0.66492101) q[0];
sx q[0];
rz(2.1524236) q[0];
x q[1];
rz(0.64935342) q[2];
sx q[2];
rz(-2.1198065) q[2];
sx q[2];
rz(1.8446814) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6886544) q[1];
sx q[1];
rz(-2.3804745) q[1];
sx q[1];
rz(1.2378986) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2003464) q[3];
sx q[3];
rz(-1.660957) q[3];
sx q[3];
rz(0.7490292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3300276) q[2];
sx q[2];
rz(-1.4469688) q[2];
sx q[2];
rz(-0.88328973) q[2];
rz(2.6106994) q[3];
sx q[3];
rz(-0.58404946) q[3];
sx q[3];
rz(-2.365621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74060488) q[0];
sx q[0];
rz(-2.8477836) q[0];
sx q[0];
rz(0.64939943) q[0];
rz(-2.7679288) q[1];
sx q[1];
rz(-0.83761907) q[1];
sx q[1];
rz(-1.9992453) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72066618) q[0];
sx q[0];
rz(-2.775179) q[0];
sx q[0];
rz(1.1144251) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9317103) q[2];
sx q[2];
rz(-1.58605) q[2];
sx q[2];
rz(1.588087) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.025665034) q[1];
sx q[1];
rz(-2.2916404) q[1];
sx q[1];
rz(1.1551108) q[1];
rz(-pi) q[2];
rz(1.1116696) q[3];
sx q[3];
rz(-1.5304426) q[3];
sx q[3];
rz(-0.99276517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.00804) q[2];
sx q[2];
rz(-2.4456761) q[2];
sx q[2];
rz(-0.35888654) q[2];
rz(0.60838962) q[3];
sx q[3];
rz(-1.5865934) q[3];
sx q[3];
rz(2.2801094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.1818039) q[0];
sx q[0];
rz(-0.53878468) q[0];
sx q[0];
rz(-0.96610075) q[0];
rz(-0.44202647) q[1];
sx q[1];
rz(-1.288488) q[1];
sx q[1];
rz(-2.7332222) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3672542) q[0];
sx q[0];
rz(-1.5440732) q[0];
sx q[0];
rz(-1.7033808) q[0];
x q[1];
rz(2.8595905) q[2];
sx q[2];
rz(-0.15371727) q[2];
sx q[2];
rz(-1.745818) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8345771) q[1];
sx q[1];
rz(-0.6041827) q[1];
sx q[1];
rz(-1.8240364) q[1];
rz(-0.28520198) q[3];
sx q[3];
rz(-2.3894579) q[3];
sx q[3];
rz(-1.3507146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.24892204) q[2];
sx q[2];
rz(-2.9327124) q[2];
sx q[2];
rz(-1.199031) q[2];
rz(-1.8438953) q[3];
sx q[3];
rz(-1.0782995) q[3];
sx q[3];
rz(-0.001999438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2099828) q[0];
sx q[0];
rz(-1.4298121) q[0];
sx q[0];
rz(-2.0344875) q[0];
rz(-0.18868407) q[1];
sx q[1];
rz(-0.37982267) q[1];
sx q[1];
rz(1.3245378) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0640776) q[0];
sx q[0];
rz(-1.0594089) q[0];
sx q[0];
rz(-1.2426503) q[0];
rz(2.1679384) q[2];
sx q[2];
rz(-1.5533548) q[2];
sx q[2];
rz(-0.23250015) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.75809385) q[1];
sx q[1];
rz(-1.0304839) q[1];
sx q[1];
rz(-3.1385149) q[1];
rz(-2.0641493) q[3];
sx q[3];
rz(-2.2098118) q[3];
sx q[3];
rz(0.30917707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.3482427) q[2];
sx q[2];
rz(-2.5620354) q[2];
sx q[2];
rz(-0.14774518) q[2];
rz(-2.31367) q[3];
sx q[3];
rz(-1.4180276) q[3];
sx q[3];
rz(-2.2016321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70656908) q[0];
sx q[0];
rz(-2.7913385) q[0];
sx q[0];
rz(1.6178004) q[0];
rz(1.2500457) q[1];
sx q[1];
rz(-2.3077272) q[1];
sx q[1];
rz(0.42125431) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5180048) q[0];
sx q[0];
rz(-0.54749085) q[0];
sx q[0];
rz(-0.7579375) q[0];
rz(-pi) q[1];
rz(1.4827316) q[2];
sx q[2];
rz(-0.95092623) q[2];
sx q[2];
rz(1.6234835) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.10277412) q[1];
sx q[1];
rz(-1.5802659) q[1];
sx q[1];
rz(-3.0547583) q[1];
rz(-pi) q[2];
rz(-1.8135083) q[3];
sx q[3];
rz(-1.1303217) q[3];
sx q[3];
rz(0.8410359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0201575) q[2];
sx q[2];
rz(-0.73548663) q[2];
sx q[2];
rz(0.45041931) q[2];
rz(1.2823229) q[3];
sx q[3];
rz(-2.4560792) q[3];
sx q[3];
rz(-0.57257563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0192075) q[0];
sx q[0];
rz(-1.6078124) q[0];
sx q[0];
rz(-1.5851198) q[0];
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
rz(-1.6433257) q[3];
sx q[3];
rz(-1.4909496) q[3];
sx q[3];
rz(-0.95076233) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
