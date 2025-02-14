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
rz(-1.5768749) q[0];
sx q[0];
rz(-1.8160507) q[0];
sx q[0];
rz(-2.5817885) q[0];
rz(-1.7436104) q[1];
sx q[1];
rz(-1.3275194) q[1];
sx q[1];
rz(-1.3865857) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5450534) q[0];
sx q[0];
rz(-2.5395416) q[0];
sx q[0];
rz(1.1132147) q[0];
x q[1];
rz(0.3561147) q[2];
sx q[2];
rz(-0.31415161) q[2];
sx q[2];
rz(-0.34700307) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6416642) q[1];
sx q[1];
rz(-1.4870502) q[1];
sx q[1];
rz(0.20831627) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.54744) q[3];
sx q[3];
rz(-0.29072116) q[3];
sx q[3];
rz(0.61222044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6036512) q[2];
sx q[2];
rz(-0.62749046) q[2];
sx q[2];
rz(0.44169912) q[2];
rz(0.80080992) q[3];
sx q[3];
rz(-1.8947314) q[3];
sx q[3];
rz(-2.7504564) q[3];
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
rz(-1.2342728) q[0];
sx q[0];
rz(-1.3751625) q[0];
sx q[0];
rz(0.32670879) q[0];
rz(-1.3455343) q[1];
sx q[1];
rz(-1.6163454) q[1];
sx q[1];
rz(0.84652841) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1382111) q[0];
sx q[0];
rz(-1.6782921) q[0];
sx q[0];
rz(-1.18446) q[0];
x q[1];
rz(-0.2791762) q[2];
sx q[2];
rz(-0.90718944) q[2];
sx q[2];
rz(0.070726591) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4281171) q[1];
sx q[1];
rz(-1.5602116) q[1];
sx q[1];
rz(1.0416563) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3896265) q[3];
sx q[3];
rz(-1.3526254) q[3];
sx q[3];
rz(-0.057540992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9091866) q[2];
sx q[2];
rz(-2.0441983) q[2];
sx q[2];
rz(0.48420134) q[2];
rz(1.3062612) q[3];
sx q[3];
rz(-0.56070119) q[3];
sx q[3];
rz(1.1854712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66121286) q[0];
sx q[0];
rz(-2.4112207) q[0];
sx q[0];
rz(-2.6015306) q[0];
rz(-1.9909319) q[1];
sx q[1];
rz(-1.8692503) q[1];
sx q[1];
rz(0.59404341) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10162486) q[0];
sx q[0];
rz(-1.4007845) q[0];
sx q[0];
rz(3.0307458) q[0];
x q[1];
rz(-0.39806367) q[2];
sx q[2];
rz(-2.4467565) q[2];
sx q[2];
rz(2.8995325) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.84114) q[1];
sx q[1];
rz(-2.3256635) q[1];
sx q[1];
rz(3.0197179) q[1];
rz(-0.066917702) q[3];
sx q[3];
rz(-2.0517931) q[3];
sx q[3];
rz(-1.1851636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.322304) q[2];
sx q[2];
rz(-1.015859) q[2];
sx q[2];
rz(-0.6907531) q[2];
rz(-2.6639719) q[3];
sx q[3];
rz(-0.4126637) q[3];
sx q[3];
rz(2.7121108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2531042) q[0];
sx q[0];
rz(-0.34585837) q[0];
sx q[0];
rz(0.75230569) q[0];
rz(-2.2350156) q[1];
sx q[1];
rz(-0.38175672) q[1];
sx q[1];
rz(-0.22901542) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54628263) q[0];
sx q[0];
rz(-2.7685099) q[0];
sx q[0];
rz(-0.65600106) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53900881) q[2];
sx q[2];
rz(-1.7663562) q[2];
sx q[2];
rz(0.53949088) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.17002925) q[1];
sx q[1];
rz(-1.3336412) q[1];
sx q[1];
rz(-2.8656061) q[1];
x q[2];
rz(-0.58521478) q[3];
sx q[3];
rz(-1.2622331) q[3];
sx q[3];
rz(-1.2206961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0166246) q[2];
sx q[2];
rz(-3.0276073) q[2];
sx q[2];
rz(1.0329186) q[2];
rz(1.0659263) q[3];
sx q[3];
rz(-1.4706127) q[3];
sx q[3];
rz(1.3037995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8835939) q[0];
sx q[0];
rz(-3.0553387) q[0];
sx q[0];
rz(-1.3294543) q[0];
rz(2.6138002) q[1];
sx q[1];
rz(-1.8254435) q[1];
sx q[1];
rz(0.39906183) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2891727) q[0];
sx q[0];
rz(-1.4932844) q[0];
sx q[0];
rz(1.9231251) q[0];
rz(-pi) q[1];
rz(2.6886986) q[2];
sx q[2];
rz(-0.78467227) q[2];
sx q[2];
rz(-1.0396027) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4265824) q[1];
sx q[1];
rz(-1.2599328) q[1];
sx q[1];
rz(-2.3764627) q[1];
rz(2.2656029) q[3];
sx q[3];
rz(-1.4238947) q[3];
sx q[3];
rz(-1.5796652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.993678) q[2];
sx q[2];
rz(-2.3892011) q[2];
sx q[2];
rz(2.1368775) q[2];
rz(2.8435454) q[3];
sx q[3];
rz(-1.7883251) q[3];
sx q[3];
rz(-1.0684439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7816958) q[0];
sx q[0];
rz(-1.5157461) q[0];
sx q[0];
rz(1.9500649) q[0];
rz(-0.53939348) q[1];
sx q[1];
rz(-2.0627779) q[1];
sx q[1];
rz(2.3753812) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.51821) q[0];
sx q[0];
rz(-2.1124387) q[0];
sx q[0];
rz(0.40671273) q[0];
rz(-pi) q[1];
rz(2.4922392) q[2];
sx q[2];
rz(-2.1198065) q[2];
sx q[2];
rz(-1.8446814) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4529383) q[1];
sx q[1];
rz(-0.76111815) q[1];
sx q[1];
rz(1.903694) q[1];
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
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.811565) q[2];
sx q[2];
rz(-1.6946239) q[2];
sx q[2];
rz(-0.88328973) q[2];
rz(2.6106994) q[3];
sx q[3];
rz(-2.5575432) q[3];
sx q[3];
rz(2.365621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74060488) q[0];
sx q[0];
rz(-2.8477836) q[0];
sx q[0];
rz(-0.64939943) q[0];
rz(2.7679288) q[1];
sx q[1];
rz(-0.83761907) q[1];
sx q[1];
rz(1.9992453) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2799046) q[0];
sx q[0];
rz(-1.4122458) q[0];
sx q[0];
rz(1.9025361) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5552005) q[2];
sx q[2];
rz(-1.3609387) q[2];
sx q[2];
rz(-0.02053989) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5264633) q[1];
sx q[1];
rz(-2.3285236) q[1];
sx q[1];
rz(-2.710756) q[1];
rz(-1.4799398) q[3];
sx q[3];
rz(-0.46077076) q[3];
sx q[3];
rz(2.6449869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.00804) q[2];
sx q[2];
rz(-2.4456761) q[2];
sx q[2];
rz(2.7827061) q[2];
rz(-0.60838962) q[3];
sx q[3];
rz(-1.5865934) q[3];
sx q[3];
rz(-2.2801094) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1818039) q[0];
sx q[0];
rz(-0.53878468) q[0];
sx q[0];
rz(-0.96610075) q[0];
rz(-2.6995662) q[1];
sx q[1];
rz(-1.8531046) q[1];
sx q[1];
rz(-2.7332222) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20710564) q[0];
sx q[0];
rz(-1.4382595) q[0];
sx q[0];
rz(0.026959628) q[0];
rz(-pi) q[1];
rz(-0.14773519) q[2];
sx q[2];
rz(-1.5281754) q[2];
sx q[2];
rz(0.10382596) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6679668) q[1];
sx q[1];
rz(-1.7136116) q[1];
sx q[1];
rz(0.98167874) q[1];
rz(-pi) q[2];
rz(2.4099952) q[3];
sx q[3];
rz(-1.3773736) q[3];
sx q[3];
rz(3.132427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8926706) q[2];
sx q[2];
rz(-0.20888027) q[2];
sx q[2];
rz(-1.199031) q[2];
rz(1.8438953) q[3];
sx q[3];
rz(-1.0782995) q[3];
sx q[3];
rz(0.001999438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2099828) q[0];
sx q[0];
rz(-1.7117806) q[0];
sx q[0];
rz(-2.0344875) q[0];
rz(2.9529086) q[1];
sx q[1];
rz(-2.76177) q[1];
sx q[1];
rz(1.8170549) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.077515006) q[0];
sx q[0];
rz(-2.0821838) q[0];
sx q[0];
rz(1.8989424) q[0];
rz(-pi) q[1];
rz(0.97365427) q[2];
sx q[2];
rz(-1.5533548) q[2];
sx q[2];
rz(-2.9090925) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3894818) q[1];
sx q[1];
rz(-2.6012724) q[1];
sx q[1];
rz(-1.5656656) q[1];
x q[2];
rz(-2.5741413) q[3];
sx q[3];
rz(-0.78560116) q[3];
sx q[3];
rz(0.42450617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.79335) q[2];
sx q[2];
rz(-0.57955727) q[2];
sx q[2];
rz(2.9938475) q[2];
rz(-0.82792264) q[3];
sx q[3];
rz(-1.4180276) q[3];
sx q[3];
rz(-0.9399606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70656908) q[0];
sx q[0];
rz(-2.7913385) q[0];
sx q[0];
rz(-1.5237923) q[0];
rz(1.891547) q[1];
sx q[1];
rz(-2.3077272) q[1];
sx q[1];
rz(2.7203383) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73250153) q[0];
sx q[0];
rz(-1.2048462) q[0];
sx q[0];
rz(-2.7247695) q[0];
x q[1];
rz(1.658861) q[2];
sx q[2];
rz(-0.95092623) q[2];
sx q[2];
rz(1.5181092) q[2];
rz(-pi) q[3];
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
x q[2];
rz(-2.6700745) q[3];
sx q[3];
rz(-0.49905159) q[3];
sx q[3];
rz(-1.7744482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0201575) q[2];
sx q[2];
rz(-2.406106) q[2];
sx q[2];
rz(-0.45041931) q[2];
rz(-1.8592698) q[3];
sx q[3];
rz(-2.4560792) q[3];
sx q[3];
rz(2.569017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12238518) q[0];
sx q[0];
rz(-1.6078124) q[0];
sx q[0];
rz(-1.5851198) q[0];
rz(-0.71313329) q[1];
sx q[1];
rz(-1.5819555) q[1];
sx q[1];
rz(-3.1176007) q[1];
rz(2.0959185) q[2];
sx q[2];
rz(-1.1633506) q[2];
sx q[2];
rz(1.9539471) q[2];
rz(-0.73591821) q[3];
sx q[3];
rz(-3.0337742) q[3];
sx q[3];
rz(-0.2119457) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
