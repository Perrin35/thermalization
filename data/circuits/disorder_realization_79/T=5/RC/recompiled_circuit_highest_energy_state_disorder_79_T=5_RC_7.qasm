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
rz(1.8140732) q[1];
sx q[1];
rz(10.811364) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0836285) q[0];
sx q[0];
rz(-2.1037408) q[0];
sx q[0];
rz(-2.8468639) q[0];
x q[1];
rz(2.8459889) q[2];
sx q[2];
rz(-1.6787375) q[2];
sx q[2];
rz(-2.2578277) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6416642) q[1];
sx q[1];
rz(-1.4870502) q[1];
sx q[1];
rz(2.9332764) q[1];
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
rz(-pi) q[1];
x q[1];
rz(1.6036512) q[2];
sx q[2];
rz(-0.62749046) q[2];
sx q[2];
rz(-2.6998935) q[2];
rz(0.80080992) q[3];
sx q[3];
rz(-1.8947314) q[3];
sx q[3];
rz(-2.7504564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2342728) q[0];
sx q[0];
rz(-1.7664302) q[0];
sx q[0];
rz(-2.8148839) q[0];
rz(-1.3455343) q[1];
sx q[1];
rz(-1.5252472) q[1];
sx q[1];
rz(-0.84652841) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1382111) q[0];
sx q[0];
rz(-1.4633006) q[0];
sx q[0];
rz(1.9571327) q[0];
rz(2.2536505) q[2];
sx q[2];
rz(-1.789621) q[2];
sx q[2];
rz(-1.8162948) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4281171) q[1];
sx q[1];
rz(-1.5813811) q[1];
sx q[1];
rz(1.0416563) q[1];
x q[2];
rz(0.31383456) q[3];
sx q[3];
rz(-2.3646128) q[3];
sx q[3];
rz(-1.8557567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9091866) q[2];
sx q[2];
rz(-2.0441983) q[2];
sx q[2];
rz(-0.48420134) q[2];
rz(1.3062612) q[3];
sx q[3];
rz(-0.56070119) q[3];
sx q[3];
rz(1.1854712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66121286) q[0];
sx q[0];
rz(-2.4112207) q[0];
sx q[0];
rz(-2.6015306) q[0];
rz(1.9909319) q[1];
sx q[1];
rz(-1.2723424) q[1];
sx q[1];
rz(0.59404341) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68349379) q[0];
sx q[0];
rz(-2.9389295) q[0];
sx q[0];
rz(-2.1432102) q[0];
rz(-pi) q[1];
rz(2.743529) q[2];
sx q[2];
rz(-2.4467565) q[2];
sx q[2];
rz(-0.24206012) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.6642528) q[1];
sx q[1];
rz(-2.3788733) q[1];
sx q[1];
rz(1.4422756) q[1];
rz(-pi) q[2];
rz(2.0527127) q[3];
sx q[3];
rz(-1.6301117) q[3];
sx q[3];
rz(-2.7869567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.322304) q[2];
sx q[2];
rz(-1.015859) q[2];
sx q[2];
rz(-0.6907531) q[2];
rz(0.47762075) q[3];
sx q[3];
rz(-2.728929) q[3];
sx q[3];
rz(-2.7121108) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2531042) q[0];
sx q[0];
rz(-0.34585837) q[0];
sx q[0];
rz(-2.389287) q[0];
rz(0.90657702) q[1];
sx q[1];
rz(-0.38175672) q[1];
sx q[1];
rz(2.9125772) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54628263) q[0];
sx q[0];
rz(-0.37308274) q[0];
sx q[0];
rz(0.65600106) q[0];
rz(-0.53900881) q[2];
sx q[2];
rz(-1.7663562) q[2];
sx q[2];
rz(-0.53949088) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6743857) q[1];
sx q[1];
rz(-1.3027281) q[1];
sx q[1];
rz(1.3246791) q[1];
rz(-0.52336542) q[3];
sx q[3];
rz(-0.65306758) q[3];
sx q[3];
rz(3.0619597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0166246) q[2];
sx q[2];
rz(-3.0276073) q[2];
sx q[2];
rz(1.0329186) q[2];
rz(-2.0756663) q[3];
sx q[3];
rz(-1.67098) q[3];
sx q[3];
rz(-1.3037995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(1.2579987) q[0];
sx q[0];
rz(-0.086253919) q[0];
sx q[0];
rz(1.3294543) q[0];
rz(2.6138002) q[1];
sx q[1];
rz(-1.3161491) q[1];
sx q[1];
rz(-0.39906183) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3947537) q[0];
sx q[0];
rz(-1.2195713) q[0];
sx q[0];
rz(-0.082562692) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45289401) q[2];
sx q[2];
rz(-0.78467227) q[2];
sx q[2];
rz(2.10199) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0001155) q[1];
sx q[1];
rz(-0.85073227) q[1];
sx q[1];
rz(-1.9898371) q[1];
rz(-pi) q[2];
rz(0.19029096) q[3];
sx q[3];
rz(-2.2566593) q[3];
sx q[3];
rz(-0.13026419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1479147) q[2];
sx q[2];
rz(-0.75239158) q[2];
sx q[2];
rz(1.0047151) q[2];
rz(-0.29804722) q[3];
sx q[3];
rz(-1.3532676) q[3];
sx q[3];
rz(-2.0731488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7816958) q[0];
sx q[0];
rz(-1.6258465) q[0];
sx q[0];
rz(-1.1915278) q[0];
rz(-2.6021992) q[1];
sx q[1];
rz(-1.0788147) q[1];
sx q[1];
rz(2.3753812) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6233826) q[0];
sx q[0];
rz(-1.0291539) q[0];
sx q[0];
rz(-0.40671273) q[0];
rz(-0.79121235) q[2];
sx q[2];
rz(-0.82381604) q[2];
sx q[2];
rz(-0.3283161) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.12745652) q[1];
sx q[1];
rz(-1.3434504) q[1];
sx q[1];
rz(0.83782398) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8155072) q[3];
sx q[3];
rz(-0.38077106) q[3];
sx q[3];
rz(2.0920193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.811565) q[2];
sx q[2];
rz(-1.4469688) q[2];
sx q[2];
rz(2.2583029) q[2];
rz(-0.53089321) q[3];
sx q[3];
rz(-0.58404946) q[3];
sx q[3];
rz(-2.365621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4009878) q[0];
sx q[0];
rz(-2.8477836) q[0];
sx q[0];
rz(-0.64939943) q[0];
rz(-0.3736639) q[1];
sx q[1];
rz(-0.83761907) q[1];
sx q[1];
rz(1.9992453) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72066618) q[0];
sx q[0];
rz(-0.3664137) q[0];
sx q[0];
rz(-1.1144251) q[0];
rz(-pi) q[1];
rz(3.0685038) q[2];
sx q[2];
rz(-0.21042779) q[2];
sx q[2];
rz(-0.054191438) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.879952) q[1];
sx q[1];
rz(-1.8790207) q[1];
sx q[1];
rz(0.76517268) q[1];
x q[2];
rz(2.029923) q[3];
sx q[3];
rz(-1.5304426) q[3];
sx q[3];
rz(-2.1488275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1335527) q[2];
sx q[2];
rz(-0.69591659) q[2];
sx q[2];
rz(-2.7827061) q[2];
rz(2.533203) q[3];
sx q[3];
rz(-1.5865934) q[3];
sx q[3];
rz(0.86148328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.95978874) q[0];
sx q[0];
rz(-2.602808) q[0];
sx q[0];
rz(-0.96610075) q[0];
rz(-2.6995662) q[1];
sx q[1];
rz(-1.8531046) q[1];
sx q[1];
rz(-2.7332222) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3672542) q[0];
sx q[0];
rz(-1.5975195) q[0];
sx q[0];
rz(1.7033808) q[0];
x q[1];
rz(2.9938575) q[2];
sx q[2];
rz(-1.5281754) q[2];
sx q[2];
rz(-3.0377667) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4736259) q[1];
sx q[1];
rz(-1.7136116) q[1];
sx q[1];
rz(-0.98167874) q[1];
rz(-2.4099952) q[3];
sx q[3];
rz(-1.3773736) q[3];
sx q[3];
rz(-3.132427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.24892204) q[2];
sx q[2];
rz(-0.20888027) q[2];
sx q[2];
rz(1.199031) q[2];
rz(-1.8438953) q[3];
sx q[3];
rz(-1.0782995) q[3];
sx q[3];
rz(-0.001999438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93160981) q[0];
sx q[0];
rz(-1.4298121) q[0];
sx q[0];
rz(2.0344875) q[0];
rz(-2.9529086) q[1];
sx q[1];
rz(-0.37982267) q[1];
sx q[1];
rz(-1.3245378) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.077515006) q[0];
sx q[0];
rz(-1.0594089) q[0];
sx q[0];
rz(1.8989424) q[0];
rz(-0.021090551) q[2];
sx q[2];
rz(-2.167835) q[2];
sx q[2];
rz(1.8151547) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.81428567) q[1];
sx q[1];
rz(-1.568157) q[1];
sx q[1];
rz(-2.1111108) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0641493) q[3];
sx q[3];
rz(-0.93178082) q[3];
sx q[3];
rz(-2.8324156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.79335) q[2];
sx q[2];
rz(-2.5620354) q[2];
sx q[2];
rz(2.9938475) q[2];
rz(0.82792264) q[3];
sx q[3];
rz(-1.4180276) q[3];
sx q[3];
rz(-2.2016321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70656908) q[0];
sx q[0];
rz(-2.7913385) q[0];
sx q[0];
rz(1.5237923) q[0];
rz(1.2500457) q[1];
sx q[1];
rz(-2.3077272) q[1];
sx q[1];
rz(0.42125431) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5180048) q[0];
sx q[0];
rz(-2.5941018) q[0];
sx q[0];
rz(2.3836552) q[0];
x q[1];
rz(-0.6217072) q[2];
sx q[2];
rz(-1.499147) q[2];
sx q[2];
rz(-0.10393427) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.359671) q[1];
sx q[1];
rz(-3.0542448) q[1];
sx q[1];
rz(-3.0328301) q[1];
x q[2];
rz(2.6700745) q[3];
sx q[3];
rz(-0.49905159) q[3];
sx q[3];
rz(1.7744482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1214352) q[2];
sx q[2];
rz(-0.73548663) q[2];
sx q[2];
rz(2.6911733) q[2];
rz(-1.2823229) q[3];
sx q[3];
rz(-2.4560792) q[3];
sx q[3];
rz(0.57257563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.12238518) q[0];
sx q[0];
rz(-1.5337802) q[0];
sx q[0];
rz(1.5564729) q[0];
rz(-0.71313329) q[1];
sx q[1];
rz(-1.5819555) q[1];
sx q[1];
rz(-3.1176007) q[1];
rz(2.6789011) q[2];
sx q[2];
rz(-1.0924871) q[2];
sx q[2];
rz(0.15747216) q[2];
rz(1.4982669) q[3];
sx q[3];
rz(-1.4909496) q[3];
sx q[3];
rz(-0.95076233) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
