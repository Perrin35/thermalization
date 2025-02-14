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
rz(-1.7436104) q[1];
sx q[1];
rz(-1.3275194) q[1];
sx q[1];
rz(-1.3865857) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5450534) q[0];
sx q[0];
rz(-2.5395416) q[0];
sx q[0];
rz(2.028378) q[0];
rz(-pi) q[1];
rz(-1.4580016) q[2];
sx q[2];
rz(-1.8646282) q[2];
sx q[2];
rz(2.4217659) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0885478) q[1];
sx q[1];
rz(-1.3632208) q[1];
sx q[1];
rz(-1.4852086) q[1];
rz(-1.54744) q[3];
sx q[3];
rz(-0.29072116) q[3];
sx q[3];
rz(0.61222044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5379415) q[2];
sx q[2];
rz(-0.62749046) q[2];
sx q[2];
rz(2.6998935) q[2];
rz(0.80080992) q[3];
sx q[3];
rz(-1.8947314) q[3];
sx q[3];
rz(0.39113623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2342728) q[0];
sx q[0];
rz(-1.3751625) q[0];
sx q[0];
rz(2.8148839) q[0];
rz(-1.7960583) q[1];
sx q[1];
rz(-1.6163454) q[1];
sx q[1];
rz(2.2950642) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.617793) q[0];
sx q[0];
rz(-1.1868068) q[0];
sx q[0];
rz(-3.0256173) q[0];
rz(0.88794215) q[2];
sx q[2];
rz(-1.3519716) q[2];
sx q[2];
rz(1.3252979) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0051028) q[1];
sx q[1];
rz(-2.0999036) q[1];
sx q[1];
rz(-3.1293312) q[1];
rz(0.31383456) q[3];
sx q[3];
rz(-2.3646128) q[3];
sx q[3];
rz(1.2858359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9091866) q[2];
sx q[2];
rz(-1.0973944) q[2];
sx q[2];
rz(-0.48420134) q[2];
rz(-1.8353315) q[3];
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
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66121286) q[0];
sx q[0];
rz(-2.4112207) q[0];
sx q[0];
rz(0.54006201) q[0];
rz(-1.1506608) q[1];
sx q[1];
rz(-1.2723424) q[1];
sx q[1];
rz(0.59404341) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0399678) q[0];
sx q[0];
rz(-1.4007845) q[0];
sx q[0];
rz(-0.11084686) q[0];
rz(-1.2582905) q[2];
sx q[2];
rz(-2.2021026) q[2];
sx q[2];
rz(-0.7429276) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.81338367) q[1];
sx q[1];
rz(-1.4821308) q[1];
sx q[1];
rz(-0.81221272) q[1];
rz(1.4433617) q[3];
sx q[3];
rz(-0.48526796) q[3];
sx q[3];
rz(1.8125774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.322304) q[2];
sx q[2];
rz(-1.015859) q[2];
sx q[2];
rz(2.4508396) q[2];
rz(0.47762075) q[3];
sx q[3];
rz(-0.4126637) q[3];
sx q[3];
rz(-0.42948183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-1.2531042) q[0];
sx q[0];
rz(-2.7957343) q[0];
sx q[0];
rz(2.389287) q[0];
rz(0.90657702) q[1];
sx q[1];
rz(-0.38175672) q[1];
sx q[1];
rz(2.9125772) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2370214) q[0];
sx q[0];
rz(-1.2777877) q[0];
sx q[0];
rz(-1.3364393) q[0];
x q[1];
rz(2.7732753) q[2];
sx q[2];
rz(-0.57007705) q[2];
sx q[2];
rz(1.7961479) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4333231) q[1];
sx q[1];
rz(-2.7796942) q[1];
sx q[1];
rz(-0.72558484) q[1];
rz(0.58521478) q[3];
sx q[3];
rz(-1.2622331) q[3];
sx q[3];
rz(-1.9208966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1249681) q[2];
sx q[2];
rz(-3.0276073) q[2];
sx q[2];
rz(1.0329186) q[2];
rz(2.0756663) q[3];
sx q[3];
rz(-1.67098) q[3];
sx q[3];
rz(-1.8377931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2579987) q[0];
sx q[0];
rz(-0.086253919) q[0];
sx q[0];
rz(1.3294543) q[0];
rz(-0.52779245) q[1];
sx q[1];
rz(-1.3161491) q[1];
sx q[1];
rz(2.7425308) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2891727) q[0];
sx q[0];
rz(-1.4932844) q[0];
sx q[0];
rz(-1.2184675) q[0];
rz(-pi) q[1];
x q[1];
rz(2.40995) q[2];
sx q[2];
rz(-1.2564617) q[2];
sx q[2];
rz(0.86282496) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0001155) q[1];
sx q[1];
rz(-2.2908604) q[1];
sx q[1];
rz(1.9898371) q[1];
rz(-pi) q[2];
rz(1.3436704) q[3];
sx q[3];
rz(-0.70763125) q[3];
sx q[3];
rz(-0.16498241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1479147) q[2];
sx q[2];
rz(-2.3892011) q[2];
sx q[2];
rz(1.0047151) q[2];
rz(0.29804722) q[3];
sx q[3];
rz(-1.3532676) q[3];
sx q[3];
rz(-1.0684439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7816958) q[0];
sx q[0];
rz(-1.6258465) q[0];
sx q[0];
rz(1.1915278) q[0];
rz(2.6021992) q[1];
sx q[1];
rz(-1.0788147) q[1];
sx q[1];
rz(0.76621145) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2711002) q[0];
sx q[0];
rz(-1.2249761) q[0];
sx q[0];
rz(-0.99084164) q[0];
rz(-pi) q[1];
rz(-0.91583385) q[2];
sx q[2];
rz(-2.1127491) q[2];
sx q[2];
rz(-0.65108991) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4529383) q[1];
sx q[1];
rz(-2.3804745) q[1];
sx q[1];
rz(-1.2378986) q[1];
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
rz(-pi) q[1];
rz(-1.3300276) q[2];
sx q[2];
rz(-1.6946239) q[2];
sx q[2];
rz(-2.2583029) q[2];
rz(0.53089321) q[3];
sx q[3];
rz(-2.5575432) q[3];
sx q[3];
rz(-2.365621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(0.74060488) q[0];
sx q[0];
rz(-2.8477836) q[0];
sx q[0];
rz(2.4921932) q[0];
rz(-0.3736639) q[1];
sx q[1];
rz(-0.83761907) q[1];
sx q[1];
rz(1.9992453) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4209265) q[0];
sx q[0];
rz(-0.3664137) q[0];
sx q[0];
rz(2.0271676) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0685038) q[2];
sx q[2];
rz(-0.21042779) q[2];
sx q[2];
rz(-3.0874012) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5264633) q[1];
sx q[1];
rz(-0.81306902) q[1];
sx q[1];
rz(-0.43083664) q[1];
rz(-pi) q[2];
rz(2.029923) q[3];
sx q[3];
rz(-1.6111501) q[3];
sx q[3];
rz(-0.99276517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1335527) q[2];
sx q[2];
rz(-2.4456761) q[2];
sx q[2];
rz(0.35888654) q[2];
rz(0.60838962) q[3];
sx q[3];
rz(-1.5549992) q[3];
sx q[3];
rz(0.86148328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95978874) q[0];
sx q[0];
rz(-2.602808) q[0];
sx q[0];
rz(-2.1754919) q[0];
rz(0.44202647) q[1];
sx q[1];
rz(-1.8531046) q[1];
sx q[1];
rz(-2.7332222) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7743384) q[0];
sx q[0];
rz(-1.5975195) q[0];
sx q[0];
rz(-1.7033808) q[0];
rz(-1.5277065) q[2];
sx q[2];
rz(-1.7183964) q[2];
sx q[2];
rz(-1.6809631) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8345771) q[1];
sx q[1];
rz(-0.6041827) q[1];
sx q[1];
rz(1.3175563) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8281874) q[3];
sx q[3];
rz(-0.85581778) q[3];
sx q[3];
rz(1.7325213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.24892204) q[2];
sx q[2];
rz(-2.9327124) q[2];
sx q[2];
rz(1.199031) q[2];
rz(1.8438953) q[3];
sx q[3];
rz(-2.0632931) q[3];
sx q[3];
rz(-0.001999438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2099828) q[0];
sx q[0];
rz(-1.4298121) q[0];
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
rz(1.3281836) q[0];
sx q[0];
rz(-1.285901) q[0];
sx q[0];
rz(2.6064742) q[0];
rz(-pi) q[1];
rz(-0.021090551) q[2];
sx q[2];
rz(-2.167835) q[2];
sx q[2];
rz(-1.326438) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3894818) q[1];
sx q[1];
rz(-0.5403203) q[1];
sx q[1];
rz(-1.575927) q[1];
rz(-0.56745133) q[3];
sx q[3];
rz(-0.78560116) q[3];
sx q[3];
rz(-0.42450617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.79335) q[2];
sx q[2];
rz(-0.57955727) q[2];
sx q[2];
rz(-0.14774518) q[2];
rz(-2.31367) q[3];
sx q[3];
rz(-1.7235651) q[3];
sx q[3];
rz(2.2016321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4350236) q[0];
sx q[0];
rz(-0.35025418) q[0];
sx q[0];
rz(1.6178004) q[0];
rz(1.891547) q[1];
sx q[1];
rz(-0.83386546) q[1];
sx q[1];
rz(0.42125431) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68116248) q[0];
sx q[0];
rz(-1.1831076) q[0];
sx q[0];
rz(-1.1739385) q[0];
rz(-1.4827316) q[2];
sx q[2];
rz(-2.1906664) q[2];
sx q[2];
rz(1.6234835) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.359671) q[1];
sx q[1];
rz(-3.0542448) q[1];
sx q[1];
rz(-3.0328301) q[1];
rz(-0.45205595) q[3];
sx q[3];
rz(-1.789942) q[3];
sx q[3];
rz(-0.62458405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1214352) q[2];
sx q[2];
rz(-0.73548663) q[2];
sx q[2];
rz(0.45041931) q[2];
rz(-1.8592698) q[3];
sx q[3];
rz(-2.4560792) q[3];
sx q[3];
rz(-0.57257563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[2];
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
rz(-2.0959185) q[2];
sx q[2];
rz(-1.9782421) q[2];
sx q[2];
rz(-1.1876456) q[2];
rz(-2.4056744) q[3];
sx q[3];
rz(-0.10781846) q[3];
sx q[3];
rz(2.929647) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
